"""
Take a VCF and create a gemini compatible database
"""
import sys

import itertools as it
import re
import zlib
import cPickle
import time
from collections import defaultdict

import numpy as np
import sqlalchemy as sql
from pedagree import Ped
import geneimpacts
from sqlalchemy import String, Float, Integer, Boolean
import cyvcf2

import cProfile
import StringIO
import pstats
import contextlib

import multiprocessing


ALT_MAX = 30
REF_MAX = 25

def grouper(n, iterable):
    iterable = iter(iterable)
    piece = list(it.islice(iterable, n))
    while piece:
        yield piece
        piece = list(it.islice(iterable, n))

@contextlib.contextmanager
def profiled():
    pr = cProfile.Profile()
    pr.enable()
    yield
    pr.disable()
    s = StringIO.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('time')
    ps.print_stats(100)
    # uncomment this to see who's calling what
    # ps.print_callers()
    print s.getvalue()

def set_column_length(column, length):
    #if t.engine
    table = column.table
    e = table.metadata.bind
    e.dispose()
    with e.connect() as e:
        c = table.columns[column.name]
        if c.type.length >= length:
            return
        c.type.length = length
        if e.dialect.name == "postgres":
            e.execute('ALTER TABLE %s ALTER COLUMN %s TYPE VARCHAR(%d)' %
                                (table.name, c.name, length))
        elif e.dialect.name == "mysql":
            e.execute('ALTER TABLE %s MODIFY %s VARCHAR(%d)' %
                                (table.name, c.name, length))

def unpack_blob(blob):
    return cPickle.loads(zlib.decompress(blob))

def xpack_blob(obj):
    if obj is None: return ''
    return buffer(blosc.compress(obj.tostring(), obj.dtype.itemsize, clevel=5, shuffle=True))

def pack_blob(obj, _none=zlib.compress(cPickle.dumps(None, cPickle.HIGHEST_PROTOCOL))):
    if obj is None: return _none
    return zlib.compress(cPickle.dumps(obj, cPickle.HIGHEST_PROTOCOL), 1)

def clean(name):
    """
    turn a vcf id into a db name
    """
    return name.replace("-", "_").replace(" ", "_").strip('"').strip("'")

def info_parse(line,
        _patt=re.compile("(\w+)=(\"[^\"]+\"|[^,]+)")):
    """
    >>> ret = info_parse('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">')
    >>> assert ret == {"ID": "AC", "Number": "A", "Type": "Integer","Description": '"Allele count in genotypes, for each ALT allele, in the same order as listed"'}, ret
    """
    assert line.startswith("##INFO=")
    stub = line.split("=<")[1].rstrip(">")
    return dict(_patt.findall(stub))


type_lookups = {
        "Integer": sql.Integer(),
        "Float": sql.Float(),
        "Flag": sql.Boolean(),
        "Character": sql.String(1),
        "String": sql.String(5),
        }

class VCFDB(object):
    gt_cols = ("gts", "gt_types", "gt_phases", "gt_depths", "gt_ref_depths",
               "gt_alt_depths", "gt_quals")

    effect_list = ["CSQ", "ANN", "EFF"]
    _black_list = []

    def __init__(self, vcf_path, db_path, ped_path=None, blobber=pack_blob,
            black_list=None):
        self.vcf_path = vcf_path
        if not db_path.startswith(("sqlite:", "mysql", "postgres")):
                db_path = "sqlite:///" + db_path
        self.db_path = db_path
        self.engine = sql.create_engine(db_path)
        self.impacts_headers = {}
        self.metadata = sql.MetaData(bind=self.engine)
        self.blobber = blobber
        self.ped_path = ped_path
        self.black_list = list(VCFDB._black_list) + list(VCFDB.effect_list) + (black_list or [])

        self.vcf = cyvcf2.VCF(vcf_path)
        self.pool = multiprocessing.Pool(2)
        # we use the cache to infer the lengths of string fields.
        self.cache = list(it.islice(self.vcf, 10000))
        self.create_columns()
        self.samples = self.create_samples()
        self.load()
        self.index()

    def _set_variant_properties(self, v, d):
        d['type'] = v.var_type
        d['sub_type'] = v.var_subtype
        d['call_rate'] = v.call_rate
        d['num_hom_ref'] = v.num_hom_ref
        d['num_het'] = v.num_het
        d['num_hom_alt'] = v.num_hom_alt
        d['aaf'] = v.aaf
        # TODO:
        #d['hwe'] = v.hwe
        #d['inbreeding_coef'] = ??
        #d['pi'] = ??

    def _load(self, iterable, create, start):

        variants = []
        keys = set()
        last = False
        i = 0
        for i, v in enumerate(iterable, start=start):
            #if i > 1000: break
            d = dict(v.INFO)
            for c in self.gt_cols:
                # named gt_bases in cyvcf2 and gts in db
                arr = v.gt_bases if c == "gts" else getattr(v, c, None)
                if arr is not None:
                    arr = arr[self.sample_idxs]
                #d[c] = self.blobber(arr)
                d[c] = arr

            d['chrom'], d['start'], d['end'] = v.CHROM, v.start, v.end
            d['ref'], d['alt'] = v.REF, ",".join(v.ALT)

            d['qual'], d['filter'], d['vcf_id'] = v.QUAL, v.FILTER, v.ID
            d['variant_id'] = i
            self._set_variant_properties(v, d)

            # TODO: just save required keys outside.
            keys.update(d.keys())

            variants.append(d)
            # http://docs.sqlalchemy.org/en/latest/faq/performance.html
            if not create and (i % 10000) == 0:
                self.insert(variants, keys, i)
                variants = variants[:0]

        if len(variants) != 0:
            self.insert(variants, keys, i, create=create)

        return i

    def load(self):
        self.t0 = self.t = time.time()

        i = self._load(self.cache, create=True, start=1)
        self.cache = []
        #with profiled():
        self._load(self.vcf, create=False, start=i+1)
        self.proc.join()

    @property
    def variants_inserter(self):
        try:
            return self._variants_inserter
        except AttributeError:
            self._variants_inserter = self.variants.insert()
            return self._variants_inserter

    @property
    def variant_impacts_inserter(self):
        try:
            return self._variant_impacts_inserter
        except AttributeError:
            self._variant_impacts_inserter = self.variant_impacts.insert()
            return self._variant_impacts_inserter

    def check_column_lengths(self, dicts, cols):
        change_cols = defaultdict(int)
        for name, c in cols.items():
            l = c.type.length
            for d in dicts:
                if len(d.get(name) or '') > l:
                    change_cols[c.name] = max(change_cols[c.name], len(d.get(name)))
                    #set_column_length(c, len(d.get(name)))
        change_cols = dict(change_cols)
        for c, l in change_cols.items():
            sys.stderr.write("changing varchar field '%s' to length %d\n" %
                                     (c,  l))
        return change_cols



    def insert(self, variants, keys, i, create=False):
        variant_impacts = []
        ivariants = []
        """
        for variant, impacts in self.pool.imap(gene_info, ((v,
            self.impacts_headers, self.blobber, self.gt_cols) for
                                                   v in variants), 2000):
            variant_impacts.extend(impacts)
            ivariants.append(variant)
        """
        for v in variants:
            variant, impacts = gene_info((v, self.impacts_headers, self.blobber,
                self.gt_cols))
            variant_impacts.extend(impacts)
            ivariants.append(variant)

        variants = ivariants
        vlengths = vilengths = {}
        if create:
            self.create(variants, variant_impacts)

        elif self.engine.dialect.name != "sqlite":
            with self.refresh():
                vlengths = self.check_column_lengths(variants, {c.name: c for c in self.variants_columns if
                    c.type.__class__.__name__ == "String"})

                vilengths = self.check_column_lengths(variant_impacts, {c.name: c for c in
                    self.variant_impacts_columns if c.type.__class__.__name__ ==
                    "String"})

        def update(d, keys):
            # fill empty columns
            u = dict.fromkeys(keys)
            u.update(d)
            # set aaf to default to -1.
            for k in u:
                if k.startswith(("af_", "aaf_")) or k.endswith(("_af", "_aaf")):
                    if u[k] is None:
                        u[k] = -1.0
            return u

        print vlengths
        self._insert(vlengths, [update(v, keys) for v in variants],
                     vilengths, variant_impacts)
        #if hasattr(self, "proc2"):
        #    self.proc2.join()
        #self.proc2 = self._insert(self.variant_impacts_inserter, variant_impacts)
        vt = vit = 0
        fmt = "%d variant_impacts:%d\tchunk time:%.1f\tv:%.2f\tvi:%.2f\t%.2f variants/second\n"
        vps = i / float(time.time() - self.t0)
        sys.stderr.write(fmt % (i, len(variant_impacts), time.time() - self.t,
            vt, vit, vps))
        self.t = time.time()

    def _insert(self, vlengths, v_objs, vilengths, vi_objs):
        if hasattr(self, "proc"):
            self.proc.join()

        #self.refresh()

        for name, clen in vlengths.items():
            col = self.variants.columns[name]
            set_column_length(col, clen)

        for name, clen in vilengths.items():
            col = self.variant_impacts.columns[name]
            set_column_length(col, clen)

        self.proc = multiprocessing.Process(target=_insert, args=((v_objs, vi_objs, self.engine),))
        self.proc.start()
        return
        tx = time.time()
        try:
            # (2006, 'MySQL server has gone away'
            # if you see this, need to increase max_allowed_packet and/or other
            # params in my.cnf (or we should detect and reduce the chunk size)
            if len(objs) > 6000:
                for group in grouper(5000, objs):
                    self.engine.execute(stmt, list(group))
            else:
                self.engine.execute(stmt, objs)
            return time.time() - tx
        except:
            for o in objs:
                self.engine.execute(stmt, o)
            raise

    def create_columns(self):
        self.variants_columns = self.get_variants_columns()
        self.variant_impacts_columns = self.get_variant_impacts_columns()

    def create(self, dvariants, dvariant_impacts):
        # update the lengths of the string columns based on the variants that
        # we've seen so far
        v_cols = {c.name: c for c in self.variants_columns if c.type.__class__.__name__ == "String"}
        self._create(dvariants, v_cols)
        for c in self.variants_columns:
            if c.name.endswith("_af") or c.name.startswith("af_"):
                print c.name, c.type, c.default

        vi_cols = {c.name: c for c in self.variant_impacts_columns if c.type.__class__.__name__ == "String"}
        self._create(dvariant_impacts, vi_cols)

        self._create_tables()

    def _create(self, dicts, cols):
        exclude_cols = set()
        for name, col in cols.items():
            if name in exclude_cols: continue
            for d in dicts:
                try:
                    if col.type.length < len(d.get(name) or ''):
                        #col.type.length = int(1.618 * len(d[name]) + 0.5)
                        col.type.length = int(1.2 * len(d[name]) + 0.5)
                except:
                    print name, col.type.length
                    raise
                if col.type.length > 200:
                    col.type = sql.TEXT()
                    exclude_cols.add(name)
                    break

    def _create_tables(self):
        # sets the column lengths for strings
        self.variants = sql.Table("variants", self.metadata, *self.variants_columns)
        for c in self.variants.columns:
            if c.name.endswith(("_af", "_aaf")) or c.name.startswith(("af_", "aaf_")):
                c.default, c.type = sql.ColumnDefault(-1.0), Float()

        self.variant_impacts = sql.Table("variant_impacts", self.metadata, *self.variant_impacts_columns)
        self.variant_impacts.drop(self.engine, checkfirst=True)
        self.variants.drop(self.engine, checkfirst=True)

        self.variants.create()
        self.variant_impacts.create()
        self.create_vcf_header_table()

    def create_vcf_header_table(self):
        h = self.vcf.raw_header
        t = sql.Table("vcf_header", self.metadata,
                      #sql.Column("vcf_header", sql.TEXT(len(h)))
                      sql.Column("vcf_header", sql.TEXT)
                      )
        t.drop(self.engine, checkfirst=True)
        t.create(checkfirst=True)
        self.engine.execute(t.insert(), [dict(vcf_header=h)])

    def get_variant_impacts_columns(self):
        return [sql.Column("variant_id", Integer,
                           sql.ForeignKey("variants.variant_id"), nullable=False),
                ] + self.variants_gene_columns()

    def refresh(self, reflect=True):
        self.engine.dispose()
        conn = self.engine.connect()
        self.metadata = sql.MetaData(bind=conn)
        if reflect:
            self.metadata.reflect()
        return conn

    def index(self):
        sys.stderr.write("indexing ... ")
        conn = self.refresh()
        t0 = time.time()
        sql.Index("idx_variants_chrom_start", self.variants.c.chrom, self.variants.c.start).create()
        sql.Index("idx_variants_exonic", self.variants.c.is_exonic).create()
        sql.Index("idx_variants_coding", self.variants.c.is_coding).create()
        sql.Index("idx_variants_impact", self.variants.c.impact).create()
        sql.Index("idx_variants_impact_severity", self.variants.c.impact_severity).create()
        sys.stderr.write("finished in %.1f seconds...\n" % (time.time() - t0))
        sys.stderr.write("total time: in %.1f seconds...\n" % (time.time() - self.t0))

    def create_samples(self):
        p = Ped(self.ped_path)
        samples = self.vcf.samples
        cols = ['sample_id', 'family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
        idxs = []
        rows = []
        for i, s in enumerate(p.samples(), start=1):
            idxs.append(samples.index(s.sample_id))
            assert s.sample_id in samples
            if i == 0:
                cols.extend(s.attrs)
            rows.append([i, s.family_id, s.sample_id, str(s.paternal_id), str(s.maternal_id),
                '1' if s.sex == 'male' else '2' if s.sex == 'female' else '-9',
                '2' if s.affected is True else '1' if s.affected is False else '-9',
                ] + s.attrs)

        scols = [sql.Column('sample_id', Integer, primary_key=True)]
        for i, col in enumerate(cols[1:], start=1):
            vals = [r[i] for r in rows]
            l = max(len(v) for v in vals)
            scols.append(sql.Column(col, String(l)))

        t = sql.Table('samples', self.metadata, *scols)
        t.drop(self.engine, checkfirst=True)
        t.create(checkfirst=True)

        self.engine.execute(t.insert(), [dict(zip(cols, r)) for r in rows])

        # track the order to pull from the genotype fields.
        self.sample_idxs = np.array(idxs)

    def get_variants_columns(self):
        columns = self.variants_default_columns()
        columns.extend(self.variants_calculated_columns())
        columns.extend(self.variants_gene_columns())
        columns.extend(self.variants_sv_columns())
        columns.extend(self.variants_info_columns(self.vcf.raw_header))
        columns.extend(self.variants_genotype_columns())
        return columns

    def variants_default_columns(self):
        return [
            sql.Column("variant_id", Integer(), primary_key=True),
            sql.Column("chrom", String(10)),
            sql.Column("start", Integer()),
            sql.Column("end", Integer()),
            sql.Column("vcf_id", String(12)),
            #sql.Column("anno_id", Integer()),
            sql.Column("ref", String(REF_MAX)),
            sql.Column("alt", String(ALT_MAX)),
            sql.Column("qual", Float()),
            sql.Column("filter", String(10)),
           ]

    def variants_gene_columns(self):
        # all of these are also stored in the variant_impacts table.
        return [
            sql.Column("gene", String(20)),
            sql.Column("transcript", String(20)),
            sql.Column("is_exonic", Boolean()),
            sql.Column("is_coding", Boolean()),
            sql.Column("is_lof", Boolean()),
            sql.Column("is_splicing", Boolean()),
            sql.Column("exon", String(8)),
            sql.Column("codon_change", String(8)),
            sql.Column("aa_change", String(8)),
            sql.Column("aa_length", String(8)),
            sql.Column("biotype", String(50)),
            sql.Column("impact", String(20)),
            sql.Column("impact_so", String(20)),
            sql.Column("impact_severity", String(4)),
            sql.Column("polyphen_pred", String(20)),
            sql.Column("polyphen_score", Float()),
            sql.Column("sift_pred", String(20)),
            sql.Column("sift_score", Float()),
            ]


    def variants_calculated_columns(self):
        return [
            sql.Column("type", String(8)),
            sql.Column("sub_type", String(20)),
            sql.Column("call_rate", Float()),
            sql.Column("num_hom_ref", Integer()),
            sql.Column("num_het", Integer()),
            sql.Column("num_hom_alt", Integer()),
            sql.Column("aaf", Float()),
            sql.Column("hwe", Float()),
            sql.Column("inbreeding_coef", Float()),
            sql.Column("pi", Float()),
            ]

    def variants_sv_columns(self):
        return []

    def variants_genotype_columns(self):

        #return [sql.Column(name, sql.BLOB()) for name in self.gt_cols]
        return [sql.Column(name, sql.LargeBinary()) for name in self.gt_cols]

    def update_impacts_headers(self, hdr_dict):
        "keep the description so we know how to parse the CSQ/ANN fields"

        desc = hdr_dict["Description"]
        if hdr_dict["ID"] == "ANN":
            parts = [x.strip("\"'") for x in re.split("\s*\|\s*", desc.split(":", 1)[1].strip('" '))]
        elif hdr_dict["ID"] == "EFF":
            parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
        elif hdr_dict["ID"] == "CSQ":
            parts = [x.strip(" [])'(\"") for x in re.split("\||\(", desc.split(":", 1)[1].strip())]
        else:
            raise Exeception("don't know how to use %s as annotation" % hdr_dict["ID"])
        self.impacts_headers[hdr_dict["ID"]] = parts

    def variants_info_columns(self, raw_header):
        "create Column() objects for each entry in the info field"
        for l in (x.strip() for x in raw_header.split("\n")):
            if not l.startswith("##INFO"):
                continue

            d = info_parse(l)
            if d["ID"] in self.effect_list:
                self.update_impacts_headers(d)
                continue

            id = clean(d["ID"].lower())
            if d["ID"] in self.black_list or id in self.black_list:
                continue

            if id == "id":
                id = "idx"
            c = sql.Column(id, type_lookups[d["Type"]], primary_key=False)
            if id.endswith(("_af", "_aaf")) or id.startswith(("af_", "aaf_", "an_")):
                #NOTE: this isn't getting filled with -1.0 in the database.
                c = sql.Column(id, Float(), default=-1.0)
            yield c


def _insert(args):
    v_objs, vi_objs, engine = args
    #print len(objs), len(vi_objs)
    engine.dispose()
    with engine.connect() as engine:
        m = sql.MetaData(bind=engine)
        m.reflect()
        v_stmt = m.tables["variants"].insert()
        vi_stmt = m.tables["variant_impacts"].insert()
        for objs, stmt in ((v_objs, v_stmt), (vi_objs, vi_stmt)):
            try:
                # (2006, 'MySQL server has gone away'
                # if you see this, need to increase max_allowed_packet and/or other
                # params in my.cnf (or we should detect and reduce the chunk size)
                if len(objs) > 6000:
                    for group in grouper(5000, objs):
                        engine.execute(stmt, group)
                else:
                    engine.execute(stmt, objs)
            except:
                for o in objs:
                    engine.execute(stmt, o)
                raise

def gene_info(d_and_impacts_headers):
    # this is parallelized as it's only simple objects and the gene impacts
    # stuff is slow.
    d, impacts_headers, blobber, gt_cols = d_and_impacts_headers
    impacts = []
    for k in (eff for eff in ("CSQ", "ANN", "EFF") if eff in d):
        if k == "CSQ":
            impacts.extend(geneimpacts.VEP(e, impacts_headers[k], checks=False) for e in d[k].split(","))
        elif k == "ANN":
            impacts.extend(geneimpacts.SnpEff(e, impacts_headers[k]) for e in d[k].split(","))
        elif k == "EFF":
            impacts.extend(geneimpacts.OldSnpEff(e, impacts_headers[k]) for e in d[k].split(","))
        #del d[k] # save some memory

    top = geneimpacts.Effect.top_severity(impacts)
    if isinstance(top, list):
        top = top[0]
    keys = ('gene', 'transcript', 'is_exonic', 'is_coding', 'is_splicing',
            'is_lof', 'exon', 'codon_change', 'aa_change', 'aa_length',
            'biotype', 'top_consequence', 'so', 'effect_severity',
            'polyphen_pred', 'polyphen_score', 'sift_pred', 'sift_score')

    for k in keys:
        d[k] = getattr(top, k)

    d['impact'] = top.top_consequence
    d['impact_so'] = top.so
    d['impact_severity'] = top.effect_severity

    for c in gt_cols:
        d[c] = blobber(d[c])


    # TODO: check "exonic" vs is_exonic"
    # TODO: check top_consequence
    gimpacts = []
    for impact in impacts:
        #gimpacts.append({k: getattr(impact, k) for k in keys})
        gimpacts.append(dict(variant_id=d['variant_id'],
                             gene=impact.gene, transcript=impact.transcript,
                             is_exonic=impact.is_exonic, is_coding=impact.is_coding,
                             is_splicing=impact.is_splicing, is_lof=impact.is_lof,
                             exon=impact.exon, codon_change=impact.codon_change,
                             aa_change=impact.aa_change, aa_length=impact.aa_length,
                             biotype=impact.biotype, top_consequence=impact.top_consequence,
                             so=impact.so, effect_severity=impact.effect_severity,
                             polyphen_pred=impact.polyphen_pred,
                             polyphen_score=impact.polyphen_score,
                             sift_pred=impact.sift_pred,
                             sift_score=impact.sift_score))
    return d, gimpacts

if __name__ == "__main__":

    import doctest
    doctest.testmod()

    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument("VCF")
    p.add_argument("db")
    p.add_argument("ped")
    a = p.parse_args()

    v = VCFDB(a.VCF, a.db, a.ped)
