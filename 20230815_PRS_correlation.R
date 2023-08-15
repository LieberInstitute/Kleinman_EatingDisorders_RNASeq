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
