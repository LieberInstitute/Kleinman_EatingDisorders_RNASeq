# to solve Dr. Walter's question about :
# Our preliminary data suggests that aripiprazole best reduces anxiety and improves eating behavior in AN 
# but would like to understand the molecular biology mechanisms. 
# Would be terrific if it was possible for your group to look at ratio of DAD2 and 5HT1A autoreceptors and 
# post-synaptic receptors, and ratio of DA D2R to 5HT1A R as well as DAD1R compared to controls.  

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(rtracklayer)
library(recount)
library(recount.bwtool)  ## can't install
library(BiocParallel)

library(limma)
library(edgeR)

load("expr_cutoff/rse_gene_modified.Rdata")
mod = model.matrix(~Dx + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
                   data = colData(rse_gene))

#load cov_rse object
load("count_data/degradation_rse_EatingDisorders_caudate.Rdata", 
     verbose = TRUE)
cov_rse = cov_rse_caudate

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

##qsva
k = num.sv(log2(assays(cov_rse)$counts+1), mod)	# 12
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]


geneRpkm = recount::getRPKM(rse_gene, length="Length")
pd = colData(rse_gene)

# do pca on genes
pca = prcomp(t(log2(geneRpkm+1)))

# percent of variance explained by pcas
pcaVars = getPcaVars(pca)
getPcaVars(pca)[1:5]


### Differential Expression ###
table(rse_gene$Dx)
rse_gene$Dx=factor(rse_gene$Dx,c("Control","MDD","ED"))

modQsva = model.matrix(~0 + Dx + AgeDeath + overallMapRate +BMI+ qSVs[,1:2] + snpPC1+snpPC2+snpPC3, 
                       data = colData(rse_gene))
colnames(modQsva)[1:3] = levels(rse_gene$Dx)
# modQsva = model.matrix(~0 + Group + AgeDeath + overallMapRate + BMI + qSVs[,1:3], data = colData(rse_gene))
colnames(modQsva)[7:8]=c('PC1','PC2')

cmtx <- makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva)
##### GENE ######
dge = DGEList(counts = assays(rse_gene)$counts,
              genes = rowData(rse_gene))
dge = calcNormFactors(dge)
vGene = voom(dge,modQsva, plot=FALSE)

## locate interested genes by Dr. Walter
drd=which(startsWith(rowData(rse_gene)$Symbol,"DRD"))
ht=which(startsWith(rowData(rse_gene)$Symbol,"HT"))

# add the interested genes normalized expression vGene$E as the colData(rse_gene)
identical(colnames(vGene$E),rse_gene$RNum)
# [1] TRUE
colData(rse_gene)=cbind(colData(rse_gene),t(vGene$E[c(drd,ht),]))
colnames(colData(rse_gene))[95:121]=rowData(rse_gene)$Symbol[c(drd,ht)]

# vGene.d = voom(dge,cbind(modQsva,rse_gene$DRD3), plot=FALSE) # no difference with vGene
modQsva.d = model.matrix(~0 + Dx + DRD2+AgeDeath + overallMapRate +BMI+ qSVs[,1:2] + snpPC1+snpPC2+snpPC3, 
                       data = colData(rse_gene))
colnames(modQsva.d)[1:3] = levels(rse_gene$Dx)

colnames(modQsva.d)[8:9]=c('PC1','PC2')

cmtx.d <- makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva.d)
vGene.d = voom(dge,modQsva.d,plot=FALSE)
fitGene = lmFit(vGene.d)
ebGene = eBayes(contrasts.fit(fitGene,cmtx.d))
sigGeneCNT = topTable(eBayes(contrasts.fit(fitGene,cmtx.d)), coef=1:3,		# for mdd->ED and cnt->ed ("ED.Control")
                      p.value = 1, number=nrow(rse_gene))
sigGeneOV = topTable(eBayes(contrasts.fit(fitGene,cmtx.d)), coef=1:2,		# for overall (F")
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


### make in batch
gene.int=rowData(rse_gene)$Symbol[c(drd,ht)]
modQsva.list<-lapply(gene.int,function(x){
  mod=model.matrix(~0 + Dx +AgeDeath + overallMapRate +BMI+ qSVs[,1:2] + snpPC1+snpPC2+snpPC3+colData(rse_gene)[,x], 
                              data = colData(rse_gene))
  colnames(mod)[1:3] = levels(rse_gene$Dx)
  colnames(mod)[7:8]=c('PC1','PC2')
  colnames(mod)[12]=x
  return(mod)
})
names(modQsva.list)=gene.int

## manually change the "-" & "_" to "." for compitable with makeContrasts function
colnames(modQsva.list[["HTT-AS"]])[12]="HTT.AS"
colnames(modQsva.list[["HTT-AS1_1"]])[12]="HTT.AS1.1"
colnames(modQsva.list[["HTR5A-AS1"]])[12]="HTR5A.AS1"

cmtx.list <- lapply(1:27,function(x) {
  test=makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva.list[[x]])
  print(x)
  return(test)
}) 

## batch calculate the sigGeneCNT
outGeneCNT.list <- lapply(10:27,function(x){
  vGene = voom(dge,modQsva.list[[x]],plot=FALSE)
  fitGene = lmFit(vGene)
  ebGene = eBayes(contrasts.fit(fitGene,cmtx.list[[x]]))
  sigGeneCNT = topTable(eBayes(contrasts.fit(fitGene,cmtx.list[[x]])), coef=1:3,		# for mdd->ED and cnt->ed ("ED.Control")
                        p.value = 1, number=nrow(rse_gene))
  sigGeneOV = topTable(eBayes(contrasts.fit(fitGene,cmtx.list[[x]])), coef=1:2,		# for overall (F")
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
  print(x)
  print(sigGeneCNT$Symbol[which(sigGeneCNT$'q_ED-MDD'<0.05)])
  return(sigGeneCNT)
  
})

names(outGeneCNT.list)=gene.int
names(cmtx.list)=gene.int
save(gene.int,cmtx.list,modQsva.list,outGeneCNT.list,file="DE/20231012-results/DRD_HT_ratios_DE_results_20240301.rda")

## write results into xlsx file
library(writexl)
write_xlsx(outGeneCNT.list,"DE/20231012-results/DRD_HT_ratios_DE_results_20240301.xlsx")
