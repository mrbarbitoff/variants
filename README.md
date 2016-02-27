# variants
This repository contains scripts for variant analysis and interpretation enhancement

As of now, the only function supplied is making a CSV-formatted file from VCF with intelligent parsing of dbNSFP functional prediction annotation (https://sites.google.com/site/jpopgen/dbNSFP) and SnpEff effects. 

Example usage of the corresponding script:

```
csvKnit.py in.vcf fathmm-mkl.txt out.csv

    in.vcf - a fully annotated VCF-file
    fathmm-mkl.txt - fathmm-MKL output file with annotations for all SNPs - if not present, leave blank
    out.csv - filename for output CSV
```
