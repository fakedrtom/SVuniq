SVuniq
=============================================================================

Overview
========
Just trying to make a VCF parser that removes SVs (or pieces of SVs that overlap)
if they are present in public SV data (like GNOMAD, CCDG).

Work in progress...

Installation
============

Documentation
================

Citation
================

Acknowledgements
================

This is how I pulled data from gnomAD:

```
zcat gnomad_v2_sv.sites.bed.gz | awk '$5!="BND" && $7=="PASS"' | cut -f 1-3,5,29-31 | grep -v "^#" | awk 'BEGIN {print "#CHROM\tSTART\tEND\tSVTYPE\tAN\tAC\tAF"} {printf "%s\t" "%d\t" "%d\t" "%s\t" "%d\t" "%d\t" "%0.6f\n", $1, $2, $3, $4, $5, $6, $7}' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"$6"\t"$5}' | bgzip -c > gnomad.bed.gz
```

And how I combined the data:

```
cat <(zcat ccdg.bed.gz | grep -v "^#" | awk '$4 != "BND"' | awk '{print $0"\tCCDG"}') <(zcat gnomad.bed.gz | grep -v "^#" | awk '{print $0"\tgnomAD"}') | gsort /dev/stdin ~/genome_ref/genome.hg19.file | awk 'BEGIN {print "#CHROM\tSTART\tEND\tSVTYPE\tAF\tAC\tAN\tSOURCE"} {print $0}' | bgzip -c > all_svs.bed.gz
```
