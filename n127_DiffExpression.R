setwd("/home/data1/R/ED")
load("expr_cutoff/rse_gene.Rdata")
rse_gene$Dx = factor(rse_geneDx, levels(rse_gene$Dx)[c(1,3,2)] )
with(colData(rse_gene), table(Dx))

mod = model.matrix(~Dx + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
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
summary(lm(qsvBonf$x[,1] ~ pd$Dx + pd$AgeDeath + pd$mitoRate + pd$rRNA_rate 
           + pd$totalAssignedGene + pd$RIN + pd$overallMapRate))



#load BMI information
bmi=read.csv('ED_BMI.csv')
bmi=bmi[match(colData(rse_gene)$BrNum,bmi$BrNum),]
identical(paste(bmi$BrNum),colData(rse_gene)$BrNum)
# [1] TRUE
colnames(bmi)[2:5]=c('AgeOnsetMdd','height','BMI','weight')
colData(rse_gene)=cbind(colData(rse_gene),bmi[,2:5])

######################

library(limma)
library(edgeR)
#########  add snpPCs #########
load("/home/data1/R/ED/snpPCs/n127_ED_MDS_indep-pairwise_200_100_0.2.rda.rda")
identical(rownames(mds),rse_gene$BrNum)
# [1] TRUE
colData(rse_gene)=cbind(colData(rse_gene),mds)

#Get rid of confounders
deleteVars <- c('mitoRate','rRNA_rate','totalAssignedGene','RIN')

rse_gene$Dx=factor(rse_gene$Dx,c("Control","MDD","ED"))
modQsva = model.matrix(~0 + Dx + AgeDeath + overallMapRate +BMI+ qSVs[,1:2] + snpPC1+snpPC2+snpPC3, 
		data = colData(rse_gene))
colnames(modQsva)[1:3] = levels(rse_gene$Dx)

colnames(modQsva)[7:8]=c('PC1','PC2')

cmtx <- makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva)
##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=FALSE)
fitGene = lmFit(vGene)
ebGene = eBayes(contrasts.fit(fitGene,cmtx))

sigGeneOV = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:2,		# for overall (F")
	p.value = 1, number=nrow(rse_gene))
	
sigGeneCNT = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:3,		# for mdd->ED and cnt->ed ("ED.Control")
	p.value = 1, number=nrow(rse_gene))

## significance levels
pvalMat = as.matrix(ebGene$p.value)
qvalMat = pvalMat
qvalMat[,1:3] = p.adjust(pvalMat[,1:3],method="fdr") 
colnames(pvalMat) = paste0("p_",colnames(pvalMat))
colnames(qvalMat) = paste0("q_",colnames(qvalMat))

# sigGeneCNT = cbind(sigGeneCNT,pvalMat,qvalMat)  
# rownames are not identical

sigGeneCNT=cbind(sigGeneCNT,pvalMat[rownames(sigGeneCNT),],qvalMat[rownames(sigGeneCNT),])

sigGeneCNT = sigGeneCNT[,-c(10,11)]

sigGeneOV = sigGeneOV[,-c(10,11)]

# for check results briefly
sigGeneCNT$Symbol[which(sigGeneCNT$'q_ED-MDD'<0.05)]
