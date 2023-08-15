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
data=data.frame(gene=exp[1,],prs1=prs[,1],age=rse_gene$AgeDeath,group=prs$group)
result=corr.test(x=data,y=NULL,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
# x could be a matrix/dataframe of expression
# y could be a matrix/dataframe of prs
# result[[1]] is the correlation matrix
# result[[4]] is the unadjusted probability values, result$p
# result[[5]] is the adjusted probability values, or result$p.adj


rownames(prs)=rse_gene$BrNum
prs.res=sapply(prs[,1:7],function(x){
  fit=lm(x~rse_gene$AgeDeath+group+mds$snpPC1+mds$snpPC2+mds$snpPC3+mds$snpPC4+mds$snpPC5,data=prs)
  return(residuals(fit))
})

dim(prs.res)
# [1] 125   7

exp=exp[,-c(5,57)]
COR=corr.test(x=t(exp),y=prs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
# take long time

library(parallel)
COR.list=mclapply(rownames(exp),function(x) {
	result=corr.test(x=exp[x,],y=prs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
	return(result)
},mc.cores=8)
COR.3snpPC=sapply(COR.list,function(x){
 return(x$p)
 })
COR.3snpPC=t(COR.3snpPC)
colnames(COR.3snpPC)=c('0.001','0.05','0.1','0.2','0.3','0.4','0.5')
rownames(COR.3snpPC)=rownames(exp)

which(COR.3snpPC==min(COR.3snpPC))%%25148
COR.3snpPC[22693,] # is the mininal exp

save(COR,exp,prs,prs.res,file='/home/data1/R/prs/ed_n127/n125_pearsonCor_result.rda')

x=COR$p
min(x)
# [1] 3.113981e-05
head(x[order(x[,3]),])
#                        0.001         0.05          0.1          0.2
# ENSG00000230673.3  0.27920174 7.689147e-04 3.113981e-05 6.315763e-05
# ENSG00000074582.12 0.04300055 4.193266e-05 5.586354e-05 1.112479e-03
length(which(x[,3]<0.01))
# [1] 391


############################## 
## run enrichment analysis ###
##############################
####
library(org.Hs.eg.db)
library(clusterProfiler)
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

id=which(COR$p[,3]<0.01)
sig=rowData(rse_gene)[id,]
sig=unique(as.character(sig$EntrezID[!is.na(sig$EntrezID)]))
length(sig)
# [1] 266

moduleGeneList =list(sig=sig)
goBP <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff  = 1, qvalueCutoff  = 1,
                       readable= TRUE)
# no significant term

############# Group Separately ##################
## eating disorder only  ###
