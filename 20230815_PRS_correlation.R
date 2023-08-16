library(enrichplot)
library(DOSE)
library(data.table)
library(magrittr)

## /home/data1/R/prs/ed_n127
files=list.files(pattern=".profile$")
files=files[8:14]

prs=lapply(files,function(x){
pr=read.table(file=x,header=TRUE)
 return(pr)
})

names(prs)=files

pr.score=sapply(files,function(x) {
	pr=prs[[x]][,6]
	return(pr)
 })

pr.score=data.frame(pr.score)
rownames(pr.score)=prs[[1]][,2]
colnames(pr.score)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')

## annotate to Brain Number
load("/home/data1/R/idcheck/refgeno_n2630_20230522.rda")
rm(refgeno)

pr.score$BrNum=brains$BrNumCor[match(rownames(pr.score),brains$ID)]
## read in PCs
pcs <- fread("n127.eigenvec", header=F) %>%
    setnames(., colnames(.), c("FID", "IID", paste0("PC",1:20)) )
identical(rownames(pr.score),pcs$IID)
#[1] TRUE
pr.score=cbind(pr.score,pcs[,3:22]
save(pr.score,file='n125_prs.rda')


### load expression data
load("/home/data1/R/ED/DE/n127_ED_bmi_2qsvsPC_3snpPCs.rda")

## align prs with rse_gene, exclude 2 controls
which(!rse_gene$BrNum %in% pr.score$BrNum)
# [1]  5 57
prs=pr.score[match(rse_gene$BrNum,pr.score$BrNum),]
exp=vGene$E
prs$Group=rse_gene$Group
prs$group=sapply(prs$Group,unclass)


library(psych)
library(writexl)
data=data.frame(gene=exp[1,],prs1=prs[,1],age=rse_gene$AgeDeath,group=prs$group)
result=corr.test(x=data,y=NULL,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
# x could be a matrix/dataframe of expression
# y could be a matrix/dataframe of prs
# result[[1]] is the correlation matrix
# result[[4]] is the unadjusted probability values, result$p
# result[[5]] is the adjusted probability values, or result$p.adj


rownames(prs)=rse_gene$BrNum
prs.res=sapply(prs[,1:7],function(x){
  fit=lm(x~rse_gene$AgeDeath+group+PC1+PC2+PC3+PC4+PC5,data=prs)
  return(residuals(fit))
})

dim(prs.res)
# [1] 125   7

exp=exp[,-c(5,57)]
# COR=corr.test(x=t(exp),y=prs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
# take long time

library(parallel)
COR.list=mclapply(rownames(exp),function(x) {
	result=corr.test(x=exp[x,],y=prs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
	return(result)
},mc.cores=8)
COR=sapply(COR.list,function(x){
 return(x$p)
 })
COR=t(COR)
colnames(COR)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(COR)=rownames(exp)

## choose the p range/cutoff of PRS has most correlation with expression
lapply(colnames(COR),function(x){
	length(which(COR[,x]<0.05))
})
# [[1]] 844 [[2]] 3552 [[3]] 3617 [[4]] 1608 [[5]] 658 [[6]] 626 [[7]] 573
lapply(colnames(COR.adj),function(x){
	length(which(COR.adj[,x]<0.05))
})
# 85 652 556 180 93 88 87
## choose COR.adj
save(COR,COR.adj,exp,prs,prs.res,file='/home/data1/R/prs/ed_n127/n125_pearsonCor_result.rda')


############################## 
## run enrichment analysis ###
##############################
####
library(org.Hs.eg.db)
library(clusterProfiler)
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

id=which(COR.adj[,2]<0.05)
sig=rowData(rse_gene)[id,]
sig=unique(as.character(sig$EntrezID[!is.na(sig$EntrezID)]))
length(sig)
# [1] 506

goBP=enrichGO(gene          = sig,
               universe      = geneUniverse,
               OrgDb         = org.Hs.eg.db,
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.1,
               qvalueCutoff  = 0.5,
               readable      = TRUE)
goMF=enrichGO(gene          = sig,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
# no sig
goCC=enrichGO(gene          = sig,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE) 
write_xlsx(list("BP"=goBP@result,"MF"=goMF@result,"CC"=goCC@result), 
				'./20230816_p.adjless0.05_PRS_GO.xlsx')


pdf('PRS exp correlation BP_GO.pdf',height=10,width=8)
dotplot(goBP,showCategory = 20,title="Expression PRS correlation p.adj<0.05 Biological Process GO")
dev.off()

pdf('PRS exp correlation MF_GO.pdf',height=10,width=8)
dotplot(goMF,showCategory = 20,title="Expression PRS correlation p.adj<0.05 Molecular Functions GO")
dev.off()

pdf('PRS exp correlation CC_GO.pdf',height=10,width=8)
dotplot(goCC,showCategory = 20,title="Expression PRS correlation p.adj<0.05 Cellular Components GO")
dev.off()

save(COR,COR.adj,exp,prs,prs.res,goBP,goMF,goCC,file='/home/data1/R/prs/ed_n127/n125_pearsonCor_result.rda')
