This code is associated with the paper from Genc et al., "Homeostatic plasticity fails at the intersection of autism-gene mutations and a novel class of common genetic modifiers". eLife, 2020. http://doi.org/10.7554/eLife.55775

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
