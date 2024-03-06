load("/home/data1/R/idcheck/refgeno_n2630_20230522.rda")
load("/home/data1/R/ED/expr_cutoff/rse_gene.Rdata")

library(SummarizedExperiment)
library(readr)
library(stringr)
library(GenomicRanges)
library(jaffelab)

pd=colData(rse_gene)

fam=read.table("/home/data1/R/genotype/n2630/ex.fam")
colnames(fam) = c("FID", "IID", "MID", "PID", "SEX","PHENO")

identical(fam$IID, brains$ID)
# [1] TRUE
fam$BrNum=brains$BrNumCor
fam = fam[!duplicated(fam$BrNum),] # remove 650 if they have 1M
# no duplication

famOut = fam[which(fam$BrNum %in% pd$BrNum),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)


bfile = "/home/data1/R/genotype/n2630/ex"
newbfile = "/home/data1/R/ED/snpPCs/n127"

## extract
system(paste("/home/data1/R/genotype/n2222/plink --bfile", bfile, 
	"--keep samples_to_extract.txt --geno 0.1 --maf 0.05 --hwe 0.000001 --make-bed --out", newbfile))

# ## independent and cluster
system(paste("/home/data1/R/genotype/n2222/plink --bfile", newbfile, "--indep 100 10 1.25 --out", newbfile))

## MDS components	
system(paste0("/home/data1/R/genotype/n2222/plink --bfile ", newbfile, 
	" --cluster --mds-plot 10 --extract ",newbfile, ".prune.in --out ", newbfile))

# ## A transpose
system(paste("/home/data1/R/genotype/n2222/plink --bfile", newbfile,
	"--recode A-transpose --out", newbfile))


################
## read in #####

## read in genotypes
genotypes  = read_delim(paste0(newbfile, ".traw"), delim="\t")

snp = as.data.frame(genotypes[,-(1:6)])
identical(colnames(snp), paste0("0_",famOut$IID))
# [1] TRUE
colnames(snp)=famOut$BrNum
snp = as.matrix(snp[,pd$BrNum])

#### read in MDS
mds = read.table(paste0(newbfile, ".mds"), 
	header=TRUE,as.is=TRUE)
identical(mds$IID,famOut$IID)
# [1] TRUE
rownames(mds)=famOut$BrNum
mds = mds[,-(1:3)]
colnames(mds) = paste0("snpPC",1:ncol(mds))

mds=mds[match(pd$BrNum,rownames(mds)),]

#############
## save #####
save(mds, snp,  compress=TRUE,
	file = "n127_ED_genotypes.rda")
save(mds, 
	file = "n127_ED_MDS.rda")



##################
## calculate caucsian only female only mds
demo=read.table("n2630_demo.txt",header=TRUE,sep='\t')
# data came from LIMS
id=match(brains$BrNumCor,paste0("Br",demo$Br))
brains$AgeDeath=demo$AgeDeath[id]
brains$Race=demo$Race[id]
brains$Sex=demo$Sex[id]

identical(fam$IID,brains$ID)
#[1] TRUE
famOut = fam[which((brains$Race %in% c('CAUC','Multi-Racial')) & brains$Sex=='F'),]

bfile = "/home/data1/R/genotype/n2630/ex"
newbfile = "/home/data1/R/ED/snpPCs/n640_cauc_f_"


#####################
# try different plink settings
famOut = fam[which(fam$BrNum %in% pd$BrNum),]
write.table(famOut[,1:2], "samples_to_extract.txt",
	col.names=FALSE, row.names=FALSE, quote=FALSE)

bfile = "/home/data1/R/genotype/n2630/ex"
newbfile = "/home/data1/R/ED/snpPCs/n127_test_"

## extract
system(paste("/home/data1/R/genotype/n2222/plink --bfile", bfile, 
	"--keep samples_to_extract.txt --geno 0.01 --maf 0.05 --hwe 0.001 --make-bed --out", newbfile))

# ## independent and cluster
system(paste("/home/data1/R/genotype/n2222/plink --bfile", newbfile, "--indep 200 100 1.25 --out", newbfile))
# change settings based on n1908 version txt file name
system(paste("/home/data1/R/genotype/n2222/plink --bfile", newbfile, "--indep-pairwise 200 100 0.2 --out", newbfile))

## MDS components	
system(paste0("/home/data1/R/genotype/n2222/plink --bfile ", newbfile, 
	" --cluster --mds-plot 10 --extract ",newbfile, ".prune.in --out ", newbfile))
