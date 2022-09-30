####################################
# from /dcl02/lieber/ptsd/RNAseq/VA_PTSD_RNAseq/R_scripts/DEG_scripts/qSV_model_analysis_DLPFC_PTSD.R
#####################################

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)

setwd('/dcl02/lieber/ajaffe/EatingDisorder_Data/differential_expression/')

#load rse_gene object
load('/dcl02/lieber/ajaffe/EatingDisorder_Data/expr_cutoff/rse_gene.Rdata')
rse_gene$Group = factor(rse_gene$Group, levels(rse_gene$Group)[c(2,3,1)] )	# reorder

## model
with(colData(rse_gene), table(Group))
# Group
     # ED     MDD Control
     # 42      43      42


# no Sex, Race, Region (all CAUC, F, Caudate)
mod = model.matrix(~Group + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
		data = colData(rse_gene))

#load cov_rse object
load("/dcl02/lieber/ajaffe/EatingDisorder_Data/count_data/degradation_rse_EatingDisorders_caudate.Rdata", 
		verbose = TRUE)
cov_rse = cov_rse_caudate

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)	# 12
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]
# [1] 44.600 11.100  5.630  3.050  2.330  1.830  1.430  1.310  1.180  0.999  0.869  0.823
# > sum(getPcaVars(qsvBonf)[1:k])
# [1] 75.151

modQsva = cbind(mod, qSVs)

save(qsvBonf, qSVs, mod, modQsva, cov_rse, file = 'rdas/ED_qsvs.Rdata')
dir.create("pdf", showWarnings = FALSE)
pdf(file='pdf/qsvs_var_explained.pdf', h=5, w=5, useDingbats = FALSE)
	plot(getPcaVars(qsvBonf)[1:k], pch=20)
dev.off()


#####################################
#Explore qSVs (adapted from https://github.com/LieberInstitute/qsva_brain/blob/master/brainseq_phase2_qsv/explore_qsvs.R)
#####################################

library(jaffelab)
library(rtracklayer)
library(recount)
library(recount.bwtool)
library(BiocParallel)
library(SummarizedExperiment)
library(readxl)
library(devtools)

setwd('/dcl02/lieber/ajaffe/EatingDisorder_Data/differential_expression/')

# gene-level expression for degradation data
load('/dcl02/lieber/ajaffe/EatingDisorder_Data/expr_cutoff/rse_gene.Rdata')
rse_gene$Group = factor(rse_gene$Group, levels(rse_gene$Group)[c(1,3,2)] )	# reorder

geneRpkm = recount::getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))

# percent of variance explained by pcas
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]
# > getPcaVars(pca)[1:5]
#[1] 25.90 15.60  6.12  4.56  3.35


## pca of degradation plots
pdf('pdf/qSV_pca.pdf', useDingbats = FALSE)
plot(pd$totalAssignedGene ~ qsvBonf$x[,1],
    ylab = "Gene Assignment Rate", pch=20,
    xlab=paste0("sv1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$RIN ~ qsvBonf$x[,1],
    ylab = "RIN", pch=20,
    xlab=paste0("sv1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$mitoRate ~ qsvBonf$x[,1],
    ylab = "mitoRate", pch=20,
    xlab=paste0("sv1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$rRNA_rate ~ qsvBonf$x[,1],
    ylab = "rRNA Rate", pch=20,
	xlab=paste0("sv1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$overallMapRate ~ qsvBonf$x[,1],
    ylab = "Overall Map Rate", pch=20,
	xlab=paste0("sv1: ",getPcaVars(qsvBonf)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
dev.off()

## pca plots
pdf('pdf/log2rpkm_pca_qual.pdf', useDingbats = FALSE)
plot(pd$totalAssignedGene ~ pca$x[,1],
    ylab = "Gene Assignment Rate", pch=20,
    xlab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$RIN ~ pca$x[,1],
    ylab = "RIN", pch=20,
    xlab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$mitoRate ~ pca$x[,1],
    ylab = "mitoRate", pch=20,
    xlab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$rRNA_rate ~ pca$x[,1],
    ylab = "rRNA Rate", pch=20,
	xlab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
plot(pd$overallMapRate ~ pca$x[,1],
    ylab = "Overall Map Rate", pch=20,
	xlab=paste0("pca1: ",getPcaVars(pca)[1],"% Var Expl"), 
	col = c("orange", "cadetblue2","magenta","black")[factor(pd$Region)])
	legend('bottomleft', sort(unique(pd$Region)), lwd=2, 
	col = c("orange", "cadetblue2","magenta","black"), bty='n')
dev.off()
