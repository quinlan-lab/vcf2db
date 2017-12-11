from __future__ import print_function
import os
from vcf2db import VCFDB, get_dburl, clean
import atexit
import sqlalchemy as sql
import sys

vcf = "tests/test.vcf"
ped = "tests/test.ped"
db = "tests/xx.db"


n_variants = sum(1 for x in open(vcf) if not x[0] == "#")

def rm(f):
    try:
        os.unlink(f)
    except OSError:
        pass

#atexit.register(rm, db)

def test_extra_impacts():
    db = "tests/xx-impacts.db"
    impacts_extras=["MQ", "QD", "HRun", "M-CAP_pred"]
    v = VCFDB(vcf, db, ped, impacts_extras=impacts_extras)
    eng = sql.create_engine(get_dburl(db))
    metadata = sql.MetaData(bind=eng)
    metadata.reflect()

    yield check_header, metadata
    yield check_samples, metadata
    yield check_variants, metadata
    yield check_aaf, metadata
    yield check_cap, metadata
    yield check_sample_genotype_counts, metadata

    vi = metadata.tables["variant_impacts"]
    ixs = list(map(clean, impacts_extras))

    # check that all of the impacts_extras are in the variant_impacts table.
    cs = {c.name: c for c in vi.columns}
    for c in ixs:
        assert c in cs

def test_load():

    v = VCFDB(vcf, db, ped)

    eng = sql.create_engine(get_dburl(db))
    metadata = sql.MetaData(bind=eng)
    metadata.reflect()

    res, = next(iter(eng.execute("select count(*) from variants")))
    assert res == n_variants, res

    assert "variants" in metadata.tables
    assert "variant_impacts" in metadata.tables
    assert "samples" in metadata.tables
    assert "vcf_header" in metadata.tables

    yield check_header, metadata
    yield check_samples, metadata
    yield check_variants, metadata
    yield check_aaf, metadata
    yield check_cap, metadata
    yield check_sample_genotype_counts, metadata

def check_cap(metadata):
    tbl = metadata.tables["variants"]
    vs = [x[0] for x in sql.select([tbl.c.m_cap_pred]).execute()]
    assert vs[0] == "testing"

def check_aaf(metadata):

    tbl = metadata.tables["variants"]
    afs = [x[0] for x in sql.select([tbl.c.max_aaf_all]).execute()]
    assert not [x for x in afs if x is None], afs


def check_header(metadata):

    vh = metadata.tables["vcf_header"]
    assert len(vh.columns) == 1


def check_samples(metadata):
    samples = metadata.tables["samples"]
    expected = [x.strip('# ') for x in next(open(ped)).split()]

    cols = [c.name for c in samples.columns]
    assert not set(expected) - set(cols),  set(expected) - set(cols)

    expected = [x.rstrip().split()[-1] for i, x in enumerate(open(ped)) if i > 0]
    obs = [x[0] for x in sql.select([samples.c.extra]).execute()]
    assert obs == expected, obs

def check_variants(metadata):

    tbl = metadata.tables["variants"]

    dels = [x[0] for x in sql.select([tbl.c.dels]).execute()]
    assert all(d == 0 for d in dels), dels


def test_load_expand():
    expand = ['gt_types', 'gt_ref_depths', 'gt_alt_depths']
    v = VCFDB(vcf, db, ped, expand=expand)

    eng = sql.create_engine(get_dburl(db))
    metadata = sql.MetaData(bind=eng)
    metadata.reflect()

    for e in expand:
        yield check_expanded_length, metadata, e

        yield check_expanded_columns, metadata, e

def check_expanded_length(metadata, field):

    assert "sample_" + field in metadata.tables

    tbl = metadata.tables['sample_' + field]
    qry = sql.select([sql.func.count(tbl.c.variant_id)])
    a_res, = metadata.bind.execute(qry).fetchone()

    vtbl = metadata.tables['variants']
    qry = sql.select([sql.func.count(vtbl.c.variant_id)])
    b_res, = metadata.bind.execute(qry).fetchone()

    assert a_res == b_res

def check_sample_genotype_counts(metadata):
    cs = metadata.tables['sample_genotype_counts']
    stmt = sql.select([cs])

    sums = []
    for s in metadata.bind.execute(stmt):
        sums.append(sum(s[1:]))

    assert len(set(sums)) == 1



def check_expanded_columns(metadata, field):

    tsamples = metadata.tables['samples']

    stmt = sql.select([tsamples.c.name])
    samples = ['%s' % v for v, in metadata.bind.execute(stmt)]

    texp = metadata.tables['sample_' + field]

    for s in samples:
        texp.columns["sample_%s" % (s, )]
