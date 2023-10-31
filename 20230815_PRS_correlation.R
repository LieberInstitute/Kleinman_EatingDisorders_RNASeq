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
#load("/home/data1/R/ED/DE/n127_ED_bmi_2qsvsPC_3snpPCs.rda")
load("/home/data1/R/ED/expr_cutoff/rse_gene_modified.Rdata") # Dx changed 3 people 20231031
load("n127_ED_bmi_2qsvs_3snpPCs.rda")

## align prs with rse_gene, exclude 2 controls
which(!rse_gene$BrNum %in% pr.score$BrNum)
# [1]  5 57
prs=pr.score[match(rse_gene$BrNum,pr.score$BrNum),]
exp=vGene$E
#prs$Group=rse_gene$Group
#prs$group=sapply(prs$Group,unclass)
prs$Group=rse_gene$Dx  # Dx changed 3 people 20231031
prs$group=as.numeric(factor(prs$Group))


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

## correlation p values original	       
COR.p=sapply(COR.list,function(x){
 return(x$p)
 })
COR.p=t(COR.p)
colnames(COR.p)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(COR.p)=rownames(exp)

## choose the p range/cutoff of PRS has most correlation with expression
lapply(colnames(COR.p),function(x){
	length(which(COR.p[,x]<0.05))
})
# [[1]] 844 [[2]] 3552 [[3]] 3617 [[4]] 1608 [[5]] 658 [[6]] 626 [[7]] 573
# 535, 2105, 1782, 857, 452, 406, 392  # Dx changed 3 people 20231031

## correlation p values adjusted	       
COR.padj=sapply(COR.list,function(x){
 return(x$p.adj)
 })
COR.padj=t(COR.padj)
colnames(COR.padj)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(COR.padj)=rownames(exp)
lapply(colnames(COR.padj),function(x){
	length(which(COR.padj[,x]<0.05))
})
# 85 652 556 180 93 88 87
# 51, 257, 198, 100, 68, 63, 62    # Dx changed 3 people 20231031
	       
## choose COR.adj
#save(COR,COR.adj,exp,prs,prs.res,file='/home/data1/R/prs/ed_n127/n125_pearsonCor_result.rda')
save(COR.list,COR.p,COR.padj,exp,prs,prs.res,file='/home/data1/R/prs/ed_n127/n125_pearsonCor_result_20231031.rda')


############
# check the correlation dirrection
COR.r=sapply(COR.list,function(x){
 return(x$r)
 })
COR.r=t(COR.r)
colnames(COR.r)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(COR.r)=rownames(exp)	



	       
############################## 
## run enrichment analysis ###
##############################
####
library(org.Hs.eg.db)
library(clusterProfiler)
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

id=which(COR[,2]<0.05 & COR.r[,2]<0) # select negative correlation genes
sig=rowData(rse_gene)[id,]
sig=unique(as.character(sig$EntrezID[!is.na(sig$EntrezID)]))
length(sig)
# [1] 506 for total, 84 for negative and positive correlation 422 genes

###### top 3 correlated genes
COR.p=sapply(COR.list,function(x){
 return(x$p)
 })
COR.p=t(COR.p)
colnames(COR.p)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(COR.p)=rownames(exp)
top.cor=data.frame(ensg=rownames(exp),Symbol=rowData(rse_gene)$Symbol,EntrezID=rowData(rse_gene)$EntrezID,cor.r=COR.r[,2],cor.p=COR.p[,2],cor.p.adj=COR[,2])
top.cor=top.cor[order(top.cor$cor.p),]
write_xlsx(list("PRS_exp_COR"=top.cor),'PRS_exp_COR_prange0.05.xlsx')
	       
	       
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
genes=data.frame(cbind(rowData(rse_gene)[id,"Symbol"],COR[id,2],COR.r[id,2]))
colnames(genes)=c('Symbol','p.adj','r')
write_xlsx(list("gene"=genes,"BP"=goBP@result,"MF"=goMF@result,"CC"=goCC@result), 
				'./20230824_p.adjless0.05_negCor_PRS_GO.xlsx')


pdf('PRS exp correlation BP_GO.pdf',height=10,width=8)
dotplot(goBP,showCategory = 20,title="Expression PRS correlation p.adj<0.05 Biological Process GO")
dev.off()

pdf('PRS exp correlation MF_GO.pdf',height=10,width=8)
dotplot(goMF,showCategory = 20,title="Expression PRS correlation p.adj<0.05 Molecular Functions GO")
dev.off()

pdf('PRS exp correlation CC_GO.pdf',height=10,width=8)
dotplot(goCC,showCategory = 20,title="Expression PRS correlation p.adj<0.05 Cellular Components GO")
dev.off()

save(COR,COR.p,COR.r,exp,prs,prs.res,goBP,goMF,goCC,file='/home/data1/R/prs/ed_n127/n125_pearsonCor_result.rda')


######## ED only ###########
prs$AgeDeath=rse_gene$AgeDeath
edexp=vGene$E[,colData(rse_gene)$Group=='ED']

edonlyprs.res=sapply(prs[,1:7],function(x){
  fit=lm(x~AgeDeath+PC1+PC2+PC3+PC4+PC5,data=prs,subset=(prs$Group=='ED'))
  return(residuals(fit))
})
# CORed=corr.test(x=edexp,y=edonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
CORed.list=mclapply(rownames(edexp),function(x) {
	result=corr.test(x=edexp[x,],y=edonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
	return(result)
},mc.cores=8)
CORed=sapply(CORed.list,function(x){
 return(x$p.adj)
 })
CORed=t(CORed)
colnames(CORed)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(CORed)=rownames(edexp)
# choose the most genes column
lapply(colnames(CORed),function(x){
	length(which(CORed[,x]<0.05))
})
# 71 177 61 30 22 18 20

######## MDD only ###########

mddexp=vGene$E[,colData(rse_gene)$Group=='MDD']

mddonlyprs.res=sapply(prs[,1:7],function(x){
  fit=lm(x~AgeDeath+PC1+PC2+PC3+PC4+PC5,data=prs,subset=(prs$Group=='MDD'))
  return(residuals(fit))
})
# CORed=corr.test(x=edexp,y=edonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
CORmdd.list=mclapply(rownames(mddexp),function(x) {
	result=corr.test(x=mddexp[x,],y=mddonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
	return(result)
},mc.cores=8)
CORmdd=sapply(CORmdd.list,function(x){
 return(x$p.adj)
 })
CORmdd=t(CORmdd)
colnames(CORmdd)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(CORmdd)=rownames(mddexp)
# choose the most genes column
lapply(colnames(CORmdd),function(x){
	length(which(CORmdd[,x]<0.05))
})
# 41 222 164 160 152 152 150
