library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)

# load rse_gene data in AMD
load('expr_cutoff/rse_gene.Rdata')

# read in DSM4 Eating disorder subtype information from Dropbox
subtype=read.table("ED_subtype.txt",header=TRUE)
colData(rse_gene)$ED.subtype=subtype$subtype[match(rse_gene$BrNum,subtype$BrNum)]
