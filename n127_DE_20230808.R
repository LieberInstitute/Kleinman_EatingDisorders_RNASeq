setwd("/home/data1/R/ED")
load("expr_cutoff/rse_gene.Rdata")
rse_gene$Group = factor(rse_gene$Group, levels(rse_gene$Group)[c(1,3,2)] )
with(colData(rse_gene), table(Group))

mod = model.matrix(~Group + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
                   data = colData(rse_gene))


load("count_data/degradation_rse_EatingDisorders_caudate.Rdata")
cov_rse = cov_rse_caudate

qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
k = num.sv(log2(assays(cov_rse)$counts+1), mod)	# 12
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
# [1] 44.600 11.100  5.630  3.050  2.330  1.830  1.430  1.310  1.180  0.999  0.869  0.823

modQsva = cbind(mod, qSVs)
library(rtracklayer)
library(recount)

library(recount.bwtool) 
library(BiocParallel)

pd = colData(rse_gene)
geneRpkm = recount::getRPKM(rse_gene, length="Length")
pca = prcomp(t(log2(geneRpkm+1)))
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]
# [1] 25.90 15.60  6.12  4.56  3.35
summary(lm(qsvBonf$x[,1] ~ pd$Group + pd$AgeDeath + pd$mitoRate + pd$rRNA_rate 
           + pd$totalAssignedGene + pd$RIN + pd$overallMapRate))

Call:
lm(formula = qsvBonf$x[, 1] ~ pd$Group + pd$AgeDeath + pd$mitoRate + 
    pd$rRNA_rate + pd$totalAssignedGene + pd$RIN + pd$overallMapRate)

Residuals:
     Min       1Q   Median       3Q      Max 
-10.1524  -2.4184   0.1897   2.3378  11.5867 

Coefficients:
                       Estimate Std. Error t value Pr(>|t|)    
(Intercept)           1.404e+02  1.060e+01  13.247  < 2e-16 ***
pd$GroupMDD           3.498e-01  9.041e-01   0.387 0.699548    
pd$GroupED            9.481e-01  9.212e-01   1.029 0.305506    
pd$AgeDeath          -4.512e-03  2.081e-02  -0.217 0.828713    
pd$mitoRate           1.210e+02  2.897e+01   4.177 5.68e-05 ***
pd$rRNA_rate         -6.727e+04  1.879e+04  -3.579 0.000501 ***
pd$totalAssignedGene -2.166e+02  8.314e+00 -26.055  < 2e-16 ***
pd$RIN               -2.232e+00  5.472e-01  -4.079 8.25e-05 ***
pd$overallMapRate    -2.990e+01  1.179e+01  -2.536 0.012518 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 3.812 on 118 degrees of freedom
Multiple R-squared:  0.9347,	Adjusted R-squared:  0.9302 
F-statistic:   211 on 8 and 118 DF,  p-value: < 2.2e-16

library(limma)
library(edgeR)
