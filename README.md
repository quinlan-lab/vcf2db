vcf2db
======

create a gemini-compatible database from a VCF:

```
python vcf2db.py diseaseX.anno.vcf.gz x.db disease_x.ped
python vcf2db.py diseaseX.anno.vcf.gz "postgres://brentp:password@localhost/gemini" disease_x.ped
python vcf2db.py diseaseX.anno.vcf.gz "mysql://brentp:password@localhost/gemini" disease_x.ped
```

With sqlite3. This inserts at about 1200 variants / second including time to index.
