#### eQTL exploring
The reference is ajaffe's brainseq2 https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_full/run_eqtls_dlpfc.R

#### libraries
```r
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
```

load expression data
```r
load("/home/data1/R/ED/expr_cutoff/rse_gene_modified.Rdata")
# extract pd and rpkms
pd = colData(rse_gene)
geneRpkm = assays(rse_gene)$rpkm
```

load SNPs genotypes data
```r

```
