vcf2db
======

create a gemini-compatible database from a VCF:

```
python vcf2db.py diseaseX.anno.vcf.gz disease_x.ped x.db
python vcf2db.py diseaseX.anno.vcf.gz disease_x.ped "postgres://brentp:password@localhost/gemini"
python vcf2db.py diseaseX.anno.vcf.gz disease_x.ped "mysql://brentp:password@localhost/gemini"
```

With sqlite3. This inserts at about 1200 variants / second including time to index.

**NOTE** while this allows loading into `mysql` and `postgres`, you will need gemini version
from github to use the database once it is loaded into `mysql` and `postgres`. Due to some [idiosyncrasies](http://docs.aws.amazon.com/efs/latest/ug/nfs4-unsupported-features.html), Amazon's Elastic File Storage (EFS) is not supported for the creation of sqlite3 databases. Elastic Block Storage (EBS) is suitable for this step.

Installation
============

```
git clone https://github.com/quinlan-lab/vcf2db
cd vcf2db
conda install -y gcc snappy # install the C library for snappy
conda install -c nlesc python-snappy=0.5
conda install -c bioconda cyvcf2 peddy
pip install -r requirements.txt
```

Annotation
==========

vcf2db now supports using bcftools csq so you can annotate with bcftools like:

```
wget ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.chr.gff3.gz
./bcftools csq --local-csq -p R -g Homo_sapiens.GRCh37.82.chr.gff3.gz -f $fasta $vcf
```

How It Works
============

Previously (and currently), gemini kept a bunch of vetted annotations along with the gemini install
and annotated an incoming VCF with those annotations as it was loaded into gemini. This is nice for
users but by *de-coupling* the annotation from the loading, we have more flexiblility.

This script pulls annotations that are defined in the INFO field, using the types defined in the header,
to create a database schema. It expects a `CSQ` tag from VEP or a `ANN` tag from snpEff in order to
determine the associated gene and consequence.

This means that the user is responsible for annotating their own VCF--though we will provide a simple
means to do this with [vcfanno](https://github.com/brentp/vcfanno). 

At this point, the script works and creates a gemini-compatible database. **It is therefore possible to use
gemini with GRCh38 or other organisms**. But the utility will depend on the resources that are available
for the given genome build and organism.

Workflow
========

Optional
--------

If there are available resources that indicate common variants, for example, 1KG or ExAC for GRCh37, then
it is useful to annotate with the allele frequencies in those populations.

To annotate a VCF with [vcfanno](https://github.com/brentp/vcfanno), follow (or use) [this example configuration](https://github.com/brentp/vcfanno/blob/4e6fd8e7f58e7e561520a6c988c91c6360a7bc42/example/gem.conf)

Functional Annotation
---------------------

Use VEP or snpEff to annotate the VCF by consequence.

Load
----

Gather a PED file and load with the script:

```
python vcf2db.py some.annotated.vcf.gz some.ped my.gemini.db
```

To have the sample fields expanded into separate tables so that they can be used INFO
SQL queries directly, use:
```
python vcf2db.py some.annotated.vcf.gz some.ped my.gemini.db --expand gt_types --expand gt_ref_depths --expand gt_alt_depths
```
