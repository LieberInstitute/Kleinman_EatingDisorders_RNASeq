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

