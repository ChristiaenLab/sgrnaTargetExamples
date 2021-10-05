# sgrnaTargetExamples
An example script of how to seach for sgRNA targets from a set of sequences. 
Requires: 
Java 8 (for FlashFry, check with `java -version`)
bedtools 2.27+
R 4.0+
R packages: GenomicRanges,rtracklayer, Biostrings, GenomicFeatures


To download requisite files:
```
git clone https://github.com/ChristiaenLab/sgrnaTargetExamples
cd sgrnaTargetExamples
make
```
If this doesn't work, the files URLs can be seen in `setup.sh`.
Usage:
```
Rscript targets.R
```
