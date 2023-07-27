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

# degradation data
load("count_data/degradation_rse_EatingDisorders_caudate.Rdata")
cov_rse = cov_rse_caudate

## model
rse_gene$Group = factor(rse_gene$Group, levels(rse_gene$Group)[c(1,3,2)] )
with(colData(rse_gene), table(Group))

## use ED.subtype == 'AN' only
table(rse_gene$Group,rse_gene$ED.subtype)
  #         AN Bulimia NOS
 # Control  0       0   0
 # MDD      0       0   0
 # ED      17       9  14
id=rse_gene$ED.subtype %in% c('Bulimia','NOS')

rse_gene.an=rse_gene[,!(id)]

mod = model.matrix(~Group + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
		data = colData(rse_gene.an))

## get qSVs for top bonferroni
identical(colnames(cov_rse),rse_gene$RNum)
# [1] TRUE
cov_rse.an=cov_rse[,!(id)]
identical(colnames(cov_rse.an),rse_gene.an$RNum)
# [1] TRUE
qsvBonf = prcomp(t(log2(assays(cov_rse.an)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse.an)$counts+1), mod)	# 10
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

modQsva = cbind(mod, qSVs)

library(rtracklayer)
library(recount)
library(recount.bwtool)  
library(BiocParallel)

geneRpkm = recount::getRPKM(rse_gene.an, length="Length")
pd = colData(rse_gene.an)

summary(lm(qsvBonf$x[,1] ~ pd$Group + pd$AgeDeath + pd$mitoRate + pd$rRNA_rate 
	+ pd$totalAssignedGene + pd$RIN + pd$overallMapRate))
#Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
#(Intercept)           1.491e+02  1.249e+01  11.930  < 2e-16 ***
#pd$GroupMDD           2.665e-01  9.330e-01   0.286 0.775803    
#pd$GroupED            1.491e+00  1.185e+00   1.258 0.211638    
#pd$AgeDeath          -1.043e-02  2.277e-02  -0.458 0.647850    
#pd$mitoRate           9.911e+01  3.118e+01   3.179 0.001995 ** 
#pd$rRNA_rate         -6.230e+04  2.228e+04  -2.797 0.006252 ** 
#pd$totalAssignedGene -2.149e+02  9.708e+00 -22.138  < 2e-16 ***
#pd$RIN               -2.289e+00  6.422e-01  -3.564 0.000574 ***
#pd$overallMapRate    -3.905e+01  1.436e+01  -2.719 0.007784 ** 


#####################################
# Case-Control 
#####################################
library(limma)
library(edgeR)

#########  add snpPCs #########
snppcs=read.table(file='./LIBD_Brain.topmed_102920.maf_0.005.info_0.8.call_0.8.uniq_brnum.sex_imputed.chr1_22.geno_0.01.maf_0.05.hwe_0.001.noambig.clump_200_100_0.2.twice.pca.txt',
sep=" ", header=TRUE)

summary(rse_gene.an$BrNum %in% snppcs$FID)
#    Mode   FALSE    TRUE 
# logical       3     101
rse_gene.an$BrNum[which(!(rse_gene.an$BrNum %in% snppcs$FID))]
# [1] "Br1060" "Br5004" "Br6254"
colData(rse_gene.an)$snpPC1=snppcs$C1[match(rse_gene.an$BrNum,snppcs$FID)]
colData(rse_gene.an)$snpPC2=snppcs$C2[match(rse_gene.an$BrNum,snppcs$FID)]
colData(rse_gene.an)$snpPC3=snppcs$C3[match(rse_gene.an$BrNum,snppcs$FID)]

#Get rid of confounders
deleteVars <- c('mitoRate','rRNA_rate','totalAssignedGene','RIN')

modQsva = model.matrix(~0 + Group + AgeDeath + overallMapRate + qSVs[,1:3] + snpPC1+snpPC2+snpPC3, 
		data = colData(rse_gene.an))
colnames(modQsva)[1:3] = levels(rse_gene.an$Group)
colnames(modQsva)[6:8]=c('PC1','PC2','PC3')

cmtx <- makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva)

##### GENE ######
dge = DGEList(counts = assays(rse_gene.an)$counts,
	genes = rowData(rse_gene.an))
dge = calcNormFactors(dge)

ids=which(!is.na(rse_gene.an$snpPC1))
dgem=dge[,ids]
vGene = voom(dgem,modQsva, plot=FALSE)

fitGene = lmFit(vGene)
ebGene = eBayes(contrasts.fit(fitGene,cmtx))

sigGeneOV = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:2,		# for overall (F")
	p.value = 1, number=nrow(rse_gene.an))
sigGeneCNT = topTable(eBayes(contrasts.fit(fitGene,cmtx)), coef=1:3,		# for mdd->ED and cnt->ed ("ED.Control")
	p.value = 1, number=nrow(rse_gene.an))

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


outGeneOV = sigGeneOV[rownames(rse_gene.an),]
outGeneCNT = sigGeneCNT[rownames(rse_gene.an),]

sum(outGeneOV$P.Value < 0.005)
# 872
sum(outGeneCNT$"p_ED-Control" < 0.005)
# 413
save.image(file='./DE/AN_contrast_gene_20230727.Rdata')

write.csv(sigGeneCNT,file='./DE/AN_3qSVsPCs_3snpPCs_noBMI_contrast_pvalues.csv')
