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
1. retrieve AN GWAS positive SNPs
```r
library(data.table)
dat<-fread("/home/data1/R/prs/pgcAN2/pgcAN2.2019-07.vcf.tsv.gz")
length(which(dat$PVAL<1e-8))
# [1] 285  # only export AN GWAS positive SNPs
pos=dat[dat$PVAL<1e-8,]
rm(dat)
```
2. query genotypes for those SNPs
```r
library(SNPlocs.Hsapiens.dbSNP149.GRCh38) # tried 155 version
snps <- SNPlocs.Hsapiens.dbSNP149.GRCh38

# two SNPs don't have rs#
which(pos$ID=='.')
#[1]  81 167, mannually update the rs# for those two SNPs
pos$ID[81]='rs3923475'
pos$ID[167]='rs4858799'
# query SNPs hg38 version position
pos.hg38= snpsById(snps,pos$ID,ifnotfound="drop")

summary(pos$ID %in% pos.hg38.155$RefSNP_id)
#   Mode   FALSE    TRUE 
# logical      32     253 
summary(pos$ID %in% pos.hg38$RefSNP_id)
#   Mode   FALSE    TRUE 
#logical      32     253

## update hg38 verison position
pos$hg38=data.frame(pos.hg38)$pos[match(pos$ID,pos.hg38$RefSNP_id)]
```
