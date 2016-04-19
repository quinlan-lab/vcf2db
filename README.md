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
from github to use the database once it is loaded into `mysql` and `postgres`.

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

Installation
============

```
pip install -r requirements.txt
```

Then use as a script.

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
python vcf2db.py some.annotated.vcf.gz my.gemini.db some.ped
```

To have the sample fields expanded into separate tables, use:
```
python vcf2db.py some.annotated.vcf.gz my.gemini.db some.ped --expand gt_types --expand gt_ref_depths --expand gt_alt_depths
```
