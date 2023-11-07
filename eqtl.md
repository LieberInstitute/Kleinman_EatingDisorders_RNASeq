#### eQTL exploring
The reference is ajaffe's brainseq2 https://github.com/LieberInstitute/brainseq_phase2/blob/master/eQTL_full/run_eqtls_dlpfc.R

#### libraries
```r
library(SummarizedExperiment)
library(jaffelab)
library(MatrixEQTL)
library(sva)
```

load expression data
```r
load("/home/data1/R/ED/expr_cutoff/rse_gene_modified.Rdata")
# extract pd and rpkms
pd = colData(rse_gene)
geneRpkm = assays(rse_gene)$rpkm
```

load SNPs genotypes data
1. retrieve AN GWAS positive SNPs
```r
library(data.table)
dat<-fread("/home/data1/R/prs/pgcAN2/pgcAN2.2019-07.vcf.tsv.gz")
length(which(dat$PVAL<1e-8))
# [1] 285  # only export AN GWAS positive SNPs
pos=dat[dat$PVAL<1e-8,]
rm(dat)
```
2. get hg38 locations for those SNPs
```r
library(SNPlocs.Hsapiens.dbSNP149.GRCh38) # tried 155 version
snps <- SNPlocs.Hsapiens.dbSNP149.GRCh38

# two SNPs don't have rs#
which(pos$ID=='.')
#[1]  81 167, mannually update the rs# for those two SNPs
pos$ID[81]='rs3923475'
pos$ID[167]='rs4858799'
# query SNPs hg38 version position
pos.hg38= snpsById(snps,pos$ID,ifnotfound="drop")

summary(pos$ID %in% pos.hg38.155$RefSNP_id)
#   Mode   FALSE    TRUE 
# logical      32     253 
summary(pos$ID %in% pos.hg38$RefSNP_id)
#   Mode   FALSE    TRUE 
#logical      32     253

## update hg38 verison position
pos$hg38=data.frame(pos.hg38)$pos[match(pos$ID,pos.hg38$RefSNP_id)]

save(pos,file="AN_GWAS_1e-8.rda")
```
3. get genotypes for those SNPs
```r
library(snpStats)
load('/home/data1/R/genotype/n2630/n2630_imputed.rda') # very big, may have conflict libraries with previous loaded, so start a new session with R

queryid=pos[!is.na(pos$hg38),]
dim(queryid)
# [1] 253  15

queryid$id=sapply(1:nrow(queryid),function(x){
  subid=which(sample[[3]]$chromosome==queryid[x,]$CHROM & sample[[3]]$position==queryid[x,]$hg38)
  return(subid)
})
# clean up unfound SNPs
queryid$n2630ID=paste(queryid$id)
which(queryid$n2630ID=="integer(0)")
# [1]  10  19  32  74  91 147 164 172 186 249

queryid$snp.name=sapply(queryid$id,function(x){
  sample[[3]][x,]$snp.name
})

queryid$allele.1=sapply(queryid$id,function(x){
  sample[[3]][x,]$allele.1
})

queryid$allele.2=sapply(queryid$id,function(x){
  sample[[3]][x,]$allele.2
})

# query imputed genotypes
snplist=list()
queryid=queryid[order(queryid$n2630ID),] # sort by queryid n2630ID
idlist=as.list(queryid$id[1:243]) # remove those unfound SNPs

snplist=lapply(idlist,function(x) {
  data.frame(as(sample[[1]][,x],'numeric'))
})
names(snplist)=queryid$ID[1:243]
genotypes=do.call('cbind',snplist)
identical(colnames(genotypes),gsub("\\:","\\.",queryid$snp.name[1:243]))
# [1] TRUE
colnames(genotypes)=queryid$ID[1:243]
```

4. covert to BrNum
```r
load('/home/data1/R/idcheck/refgeno_n2630_20230522.rda')
identical(rownames(genotypes), brains$ID)
#[1] TRUE
genotypes$BrNum=brains$BrNumCor

summary(genotypes$BrNum %in% rse_gene$BrNum)
#   Mode   FALSE    TRUE 
# logical    2503     127

n127=genotypes[match(rse_gene$BrNum,genotypes$BrNum),]
rownames(n127)=n127$BrNum
n127=t(n127[,1:243])
save(queryid,pos,n127,file="/home/data1/R/ED/DE/20231012-results/AN_GWAS_1e-8_n127_genotypes.rda")
```

### Statistic Model
```r
pd$Dx=factor(pd$Dx,levels=c('Control','ED','MDD'))
summary(pd$Dx)
# Control      ED     MDD 
#     42      41      44
mod = model.matrix(~Dx + AgeDeath + snpPC1 + snpPC2 + snpPC3,	data = pd)
```

### create SNP objects
```r
library(MatrixEQTL)
theSnps = SlicedData$new(as.matrix(n127))
theSnps$ResliceCombined(sliceSize = 50000)

snpspos = queryid[,c("ID","CHROM","hg38")][1:243,]
colnames(snpspos) = c("name","chr","pos")
```

### do PCA
```r
library(sva)
library(recount)
geneRpkm = recount::getRPKM(rse_gene, length="Length")
pcaGene = prcomp(t(log2(geneRpkm+1)))
kGene = num.sv(log2(geneRpkm+1), mod) # kGene is 11
genePCs = pcaGene$x[,1:kGene]
pcaVars = getPcaVars(pcaGene)
getPcaVars(pcaGene)[1:5]
# [1] 25.90 15.60  6.12  4.56  3.35

covsGene = SlicedData$new(t(cbind(mod[,-1],genePCs)))
```

### feature annotation
```r
posGene = as.data.frame(rowRanges(rse_gene))[,1:3]
posGene$name = rownames(posGene)
posGene = posGene[,c(4,1:3)]
```

### sliced expression data
```r
geneSlice = SlicedData$new(log2(geneRpkm+1))
geneSlice$ResliceCombined(sliceSize = 5000)

```

### Run EQTLs
```r
print("Starting eQTLs")

meGene = Matrix_eQTL_main(snps=theSnps, gene = geneSlice, 
	cvrt = covsGene, output_file_name.cis =  ".ctxt" ,
	pvOutputThreshold.cis = 1,  pvOutputThreshold=0,
	snpspos = snpspos, genepos = posGene, 
	useModel = modelLINEAR,	cisDist=5e5,
	pvalue.hist = 100,min.pv.by.genesnp = TRUE)

# extract
geneEqtl = meGene$cis$eqtls
geneEqtl$gene = as.character(geneEqtl$gene)
geneEqtl$snps = as.character(geneEqtl$snps)
# no cis sig SNP was found

## trans SNPs eqtl
length(which(meGene$trans$eqtls$FDR<0.05))
# [1] 173

transeqtl=meGene$trans$eqtls[meGene$trans$eqtls$FDR<0.05,]

# annotate
transeqtl$Symbol=rowData(rse_gene)$Symbol[match(transeqtl$gene,rowData(rse_gene)$gencodeID)]
table(transeqtl$Symbol)
# WDR6 
# 173 
```

### check genotypes frequency
```r
library(reshape2)
df=data.frame(t(n127))

# control group
dfct=df[rse_gene$Dx=='Control',]
genocountsCT=lapply(colnames(dfct),function(x){
reshape2::melt(table(dfct[,x]))
})
names(genocountsCT)=colnames(dfct)
CT.counts=do.call("rbind",genocountsCT)

# ED group
dfed=df[rse_gene$Dx=='ED',]
genocountsED=lapply(colnames(dfed),function(x){
reshape2::melt(table(dfed[,x]))
})
names(genocountsED)=colnames(dfed)
ED.counts=do.call("rbind",genocountsED)

# MDD group
dfmdd=df[rse_gene$Dx=='MDD',]
genocountsMDD=lapply(colnames(dfmdd),function(x){
reshape2::melt(table(dfmdd[,x]))
})
names(genocountsMDD)=colnames(dfmdd)
MDD.counts=do.call("rbind",genocountsMDD)
```
