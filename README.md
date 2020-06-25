# Genc2020_RNAseq
Analysis codes for Genc et al. 2020. Please find the raw files in the GEO database, the accession ID is `GSE153225`

## Expression quantification using salmon

```
# index the reference genome
salmon  index -k 23 -t Drosophila_melanogaster.BDGP6.cdna.all.fa.gz -i BDGP6.salmon 

# salmon quant
salmon quant -i BDGP6.salmon -l A -r $fastq_file -p 2 --gcBias --fldMean 350
```

## Differential gene expression analysis using DESeq2

Please run the R code `run_deseq2.R` from this repository.
