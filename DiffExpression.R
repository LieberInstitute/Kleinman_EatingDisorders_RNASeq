library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)

setwd('/users/cli/ED')  ## only work on cluster R enviorment, not Rstudio of cluster

#load rse_gene object
load('/dcl02/lieber/ajaffe/EatingDisorder_Data/expr_cutoff/rse_gene.Rdata')

# write.csv(colData(rse_gene),file='/users/cli/ED/n127_rse_gene_colData.csv')

####################  not use #############
#load BMI information
bmi=read.csv('ED_BMI.csv')
bmi=bmi[match(colData(rse_gene)$BrNum,bmi$BrNum),]
identical(paste(bmi$BrNum),colData(rse_gene)$BrNum)
# [1] TRUE
colnames(bmi)[2:5]=c('AgeOnsetMdd','height','BMI','weight')
colData(rse_gene)=cbind(colData(rse_gene),bmi[,2:5])
####################  not use #############


## model
rse_gene$Group = factor(rse_gene$Group, levels(rse_gene$Group)[c(1,3,2)] )
with(colData(rse_gene), table(Group))

mod = model.matrix(~Group + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
		data = colData(rse_gene))
#Br5525 doen's have BMI data, set as 28

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

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("recount")
BiocManager::install('LieberInstitute/recount.bwtool')

library(rtracklayer)
library(recount)
library(recount.bwtool)  ## can't install
library(BiocParallel)

geneRpkm = recount::getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))

# percent of variance explained by pcas
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]
#[1] 25.90 15.60  6.12  4.56  3.35

# top 1000 expressed regions associated with degradation
load("/dcl02/lieber/ajaffe/EatingDisorder_Data/differential_expression/rdas/ED_qsvs.Rdata", verbose = TRUE)

# ## get qSVs for top bonferroni
# qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))
getPcaVars(qsvBonf)[1:5]
#[1] 44.60 11.10  5.63  3.05  2.33

summary(lm(qsvBonf$x[,1] ~ pd$Group + pd$AgeDeath + pd$mitoRate + pd$rRNA_rate 
	+ pd$totalAssignedGene + pd$RIN + pd$overallMapRate))
	
	# Coefficients:
	#                        Estimate Std. Error t value Pr(>|t|)    
	# (Intercept)           1.392e+02  1.071e+01  12.995  < 2e-16 ***
	# pd$GroupED            1.050e+00  9.299e-01   1.129  0.26113    
	# pd$GroupMDD           3.449e-01  9.051e-01   0.381  0.70388    
	# pd$AgeDeath          -4.216e-03  2.083e-02  -0.202  0.83999    
	# pd$mitoRate           1.203e+02  2.901e+01   4.147 6.42e-05 ***
	# pd$rRNA_rate         -6.787e+04  1.883e+04  -3.605  0.00046 ***
	# pd$totalAssignedGene -2.156e+02  8.407e+00 -25.645  < 2e-16 ***
	# pd$RIN               -2.189e+00  5.501e-01  -3.979  0.00012 ***
	# pd$overallMapRate    -3.066e+01  1.183e+01  -2.590  0.01081 *  
	# pd$BMI                3.903e-02  4.543e-02   0.859  0.39203


#####################################
# Case-Control 
#####################################
library(limma)
library(edgeR)

#########  add snpPCs #########
snppcs=read.table(file='/users/cli/ED/LIBD_Brain.topmed_102920.maf_0.005.info_0.8.call_0.8.uniq_brnum.sex_imputed.chr1_22.geno_0.01.maf_0.05.hwe_0.001.noambig.clump_200_100_0.2.twice.pca.txt',
sep=" ", header=TRUE)

summary(rse_gene$BrNum %in% snppcs$FID)
#    Mode   FALSE    TRUE 
# logical       5     122
rse_gene$BrNum[which(!(rse_gene$BrNum %in% snppcs$FID))]
# [1] "Br1060" "Br5519" "Br5004" "Br6019" "Br6254"
colData(rse_gene)$snpPC1=snppcs$C1[match(rse_gene$BrNum,snppcs$FID)]
colData(rse_gene)$snpPC2=snppcs$C2[match(rse_gene$BrNum,snppcs$FID)]
colData(rse_gene)$snpPC3=snppcs$C3[match(rse_gene$BrNum,snppcs$FID)]


#Get rid of confounders
deleteVars <- c('mitoRate','rRNA_rate','totalAssignedGene','RIN')


modQsva = model.matrix(~0 + Group + AgeDeath + overallMapRate + qSVs[,1:3] + snpPC1+snpPC2+snpPC3, 
		data = colData(rse_gene))
colnames(modQsva)[1:3] = levels(rse_gene$Group)
# modQsva = model.matrix(~0 + Group + AgeDeath + overallMapRate + BMI + qSVs[,1:3], data = colData(rse_gene))
colnames(modQsva)[6:8]=c('PC1','PC2','PC3')

cmtx <- makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva)

##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
	genes = rowData(rse_gene))
dge = calcNormFactors(dge)


ids=which(!is.na(rse_gene$snpPC1))
dgem=dge[,ids]

# pdf('pdf/voom_qsva.pdf', useDingbats = FALSE)
	vGene = voom(dge,modQsva, plot=FALSE)
	vGene = voom(dgem,modQsva, plot=FALSE)
# dev.off()
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


outGeneOV = sigGeneOV[rownames(rse_gene),]
outGeneCNT = sigGeneCNT[rownames(rse_gene),]

sum(outGeneOV$P.Value < 0.005)
# 214
sum(outGeneCNT$"p_ED-Control" < 0.005)
# 311
save.image(file='ED_contrast_gene.Rdata')

write.csv(sigGeneCNT,file='3qSVsPCs_noBMI_contrast_pvalues.csv')


#### Examining the number of DE genes
# https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
summary(decideTests(ebGene))
#       MDD-Control ED-MDD ED-Control
# Down           361      0         78
# NotSig       24456  25148      25033
# Up             331      0         37

dt <- decideTests(ebGene)
## change the threshold as fdr q<0.05
dtsig = sigGeneCNT[rownames(dt),]


sigid=list()
sigid[[1]]=which(dtsig$`q_MDD-Control`<0.05)
sigid[[2]]=which(dtsig$`q_ED-MDD`<0.05)
sigid[[3]]=which(dtsig$`q_ED-Control`<0.05)
names(sigid)=c('MDD-Control','ED-MDD','ED-Control')
colnames(dtsig)[10:12]=c('MDD-Control','ED-MDD','ED-Control')

dt[,1]=0
dt[,2]=0
dt[,3]=0

# up regulation
up=which(dtsig[sigid[[1]],'MDD-Control']>0)
dt[sigid[[1]],1][up]=1

up=which(dtsig[sigid[[2]],'ED-MDD']>0)
dt[sigid[[2]],2][up]=1

up=which(dtsig[sigid[[3]],'ED-Control']>0)
dt[sigid[[3]],3][up]=1

# down regulation
down=which(dtsig[sigid[[1]],'MDD-Control']<0)
dt[sigid[[1]],1][down]=-1

down=which(dtsig[sigid[[2]],'ED-MDD']<0)
dt[sigid[[2]],2][down]=-1

down=which(dtsig[sigid[[3]],'ED-Control']<0)
dt[sigid[[3]],3][down]=-1

summary(dt)
#        MDD-Control ED-MDD ED-Control
# Down           163     18         87
# NotSig       24827  25125      25022
# Up             158      5         39

jpeg(file='vennDiagram.jpg',height=5,width=6,unit='in',res=300)
vennDiagram(dt, circle.col=c("turquoise", "salmon","purple"))
dev.off()

save.image(file='ED_contrast_gene_n122.Rdata')
