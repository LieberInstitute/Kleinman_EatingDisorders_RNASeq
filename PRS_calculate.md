### Reading the AN data file
```r
ls()library(data.table)
dat<-fread("pgcAN2.2019-07.vcf.tsv.gz")
## ED data downloaded from https://pgc.unc.edu/for-researchers/download-results/
## 14671980.zip
# Filter out SNPs
result <- dat[IMPINFO > 0.8]z
# impinfo is imputation quality, and MAF already filtered by >=0.01, ad the readme.pdf
# Output the gz file
fwrite(result, "pgcAN2.gz", sep="\t")
head(dat)
   CHROM       POS         ID REF ALT         BETA     SE    PVAL NGT IMPINFO
1:     8 101592213 rs62513865   T   C  0.009197572 0.0265 0.72760   0   0.955
2:     8 106973048 rs79643588   A   G  0.000199980 0.0230 0.99390   0   0.981
3:     8 108690829 rs17396518   G   T  0.010900374 0.0137 0.42550   0   0.958
4:     8 108681675   rs983166   C   A  0.027595712 0.0137 0.04361   0   0.966
5:     8 103044620 rs28842593   C   T -0.024302939 0.0193 0.20800   0   0.874
6:     8 109712249 rs35107696   I   D  0.018399683 0.0159 0.24680   0   0.989
   NEFFDIV2  NCAS  NCON                              DIRE
1: 23160.95 16992 55525 -+-+----+------+--+-++++-++-++++-
2: 23160.95 16992 55525 ------+--+++-+-+++---++-+-+++--++
3: 23160.95 16992 55525 +---+++++-+++-+--+--+++--++--++--
4: 23160.95 16992 55525 +---+++++++-+-++++--+-++--+--++++
5: 23160.95 16992 55525 +-+------+-+++-++-+----++----+++-
6: 23160.95 16992 55525 ++-+-+-+++---+-+---+++-+-++-++-+-

```

### remove duplicate SNPs
```bash 
gunzip -c pgcAN2.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > pgcAN2.nodup.gz
```

### remove Ambiguous SNPs
```bash
gunzip -c pgcAN2.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
gzip > pgcAN2.QC.gz

cat pgcAN2.QC.gz | wc -l
1042227 # [1] 6663431 in total
```




## Extract n127 ED cohort from n2630 data set
```bash
/home/data1/R/genotype/n2222/plink \
--bfile /home/data1/R/genotype/n2630/ex \
--keep /home/data1/R/ED/snpPCs/samples_to_extract.txt \
--make-bed --out \
/home/data1/R/prs/ed_n127/n127
```

### standard GWAS QC
```bash
/home/data1/R/genotype/n2222/plink \
    --bfile n127 \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out n127.QC

cat n127.QC.snplist | wc -l
8986425
```

### perform pruning to remove highly correlated SNPs:
```bash
/home/data1/R/genotype/n2222/plink \
    --bfile n127 \
    --keep n127.QC.fam \
    --extract n127.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out n127.QC

cat n127.QC.prune.out | wc -l
8170778
```

### Heterozygosity rates computed by prune.in
```bash
/home/data1/R/genotype/n2222/plink \
    --bfile n127 \
    --extract n127.QC.prune.in \
    --keep n127.QC.fam \
    --het \
    --out n127.QC

--------------
815647 variants and 127 people pass filters and QC.
Note: No phenotypes present.
--het: 792717 variants scanned, report written to n127.QC.het .
```
This will generate the n127.QC.het file, which contains F coefficient estimates for assessing heterozygosity. We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean, which can be performed using the following R command (assuming that you have R downloaded, then you can open an R session by typing R in your terminal):
```r
library(data.table)
# Read in file
dat <- fread("n127.QC.het")
# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], "n127.valid.sample", sep="\t") 
id=which(!dat$IID %in% valid$IID)
load("/home/data1/R/idcheck/refgeno_n2630_20230522.rda")
brains[brains$ID %in% dat$IID[id],]
#  BrNumCor
#  Br1060 control
#  Br5228 control 
# 2 samples were excluded
q() # exit R
```

### remove Mismatching SNPs
1. Load the bim file, the summary statistic and the QC SNP list into R
```r
# magrittr allow us to do piping, which help to reduce the 
# amount of intermediate data types
library(data.table)
library(magrittr)
# Read in bim file 
bim <- fread("n127.bim") %>%
    # Note: . represents the output from previous step
    # The syntax here means, setnames of the data read from
    # the bim file, and replace the original column names by 
    # the new names
    setnames(., colnames(.), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")) %>%
    # And immediately change the alleles to upper cases
    .[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
# Read in summary statistic data (require data.table v1.12.0+)
AN <- fread("../pgcAN2/pgcAN2.QC.gz") %>%
    # And immediately change the alleles to upper cases
    .[,c("REF","ALT"):=list(toupper(REF), toupper(ALT))]
colnames(AN)[1:5]=c("CHR","BP","SNP","A2","A1") # change colnames for match n127
## REF is A2, ALT is A1 according to an2019.readme.pdf

# Read in QCed SNPs
qc <- fread("n127.QC.snplist", header=F)
```

2. Identify SNPs that require strand flipping
```r
# Function for calculating the complementary allele
complement <- function(x){
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
} 

# # Merge summary statistic with target
# colnames(AN)[3]='SNPrs'
# colnames(AN)[15]='SNP' # convert snp name to chr:bp:alt:ref
# info <- merge(df19, AN, by="SNP")
# dim(info)
# # [1] 4458383      21
# summary(info$SNP %in% qc$V1)
#   Mode   FALSE    TRUE 
# logical 4420415   37968 --------------------too less, maybe use annotation methods
##############################################
## use annotation methods get rs#
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
snpcount(snps)
seqinfo(snps)

## check all the AN snps's hg38 position
## AN[132,]has no rs#, remove this line
an.hg38=snpsById(snps,AN$SNP[-132],ifnotfound="drop")
# 5861522 SNPs were found
save(an.hg38,file="../pgcAN2/AN.hg38.rda")
an.hg38.df=data.frame(an.hg38)

bim$snpname=paste(bim$CHR,bim$BP,sep=":")
an.hg38.df$snpname=paste(an.hg38.df$seqnames,an.hg38.df$pos,sep=":")
summary(an.hg38.df$snpname %in% bim$snpname)
#    Mode   FALSE    TRUE 
# logical  202757 5658765 

## merge bim and AN by chr:pos
info = bim[bim$snpname %in% an.hg38.df$snpname,]
dim(info)
# [1] 5660463       7
info$rs=an.hg38.df$RefSNP_id[match(info$snpname,an.hg38.df$snpname)]
id=match(info$rs,AN$SNP)
info$A1=AN$A1[id]
info$A2=AN$A2[id]


# ## use chr+pos to annotation n127 snp rs#, which could be find overlaps with AN data
# for(i in 1:nrow(an.hg38.df)){
#     id=which(bim$CHR==an.hg38.df$seqnames[i] & bim$BP==an.hg38.df$pos[i])
#     bim$RS[id]=an.hg38.df$RefSNP_id[i]
# } ## this takes too long time

summary(info$SNP %in% qc$V1)
#    Mode   FALSE    TRUE 
# logical  105282 5555181 

# filter QCed SNPs
info <- info[info$SNP %in% qc$V1,]
save(info, file="filtered_QCed_SNPs.rda")

# Get SNPs that have the same alleles across base and target
info.match <- info[A1 == B.A1 & A2 == B.A2, SNP]
length(info.match)
# [1] 1,204,450
length(info[A1 == B.A2 & A2 == B.A1, SNP])
# [1] 4,346,266

# Identify SNPs that are complementary between base and target
com.snps <- info[sapply(B.A1, complement) == A1 &
                    sapply(B.A2, complement) == A2, SNP]
length(com.snps)
#[1] 417        

# Now update the bim file
bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
            sapply(B.A2, complement))]
```

3. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)
```r
# identify SNPs that need recoding
recode.snps <- info[B.A1==A2 & B.A2==A1, SNP]
length(recode.snps)
# [1] 4346266
# Update the bim file
bim[SNP %in% recode.snps, c("B.A1", "B.A2") :=
        list(B.A2, B.A1)]

# identify SNPs that need recoding & complement
com.recode <- info[sapply(B.A1, complement) == A2 &
                    sapply(B.A2, complement) == A1, SNP]
length(com.recode) 
# [1] 417                   
# Now update the bim file
bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
            sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim[,c("SNP", "B.A1")], "n127.a1", col.names=F, sep="\t") # 18326354 snps remain
```

4. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)
```r
mismatch <- bim[!(SNP %in% info.match |
                    SNP %in% com.snps |
                    SNP %in% recode.snps |
                    SNP %in% com.recode), SNP]
length(mismatch)
# [1] 12773486
write.table(mismatch, "n127.mismatch", quote=F, row.names=F, col.names=F)
```

### Duplicate SNPs, make sure to remove any duplicate SNPs
```r
summary(!duplicated(bim$SNP))
#     Mode     TRUE 
#  logical 18326354
```

### Sex chromosomes
A sex check can be performed in PLINK, in which individuals are called as females if their X chromosome homozygosity estimate (F statistic) is < 0.2 and as males if the estimate is > 0.8.
```bash
/home/data1/R/genotype/n2222/plink \
--bfile n127 \
--extract n127.QC.prune.in \
--keep n127.valid.sample \
--check-sex \
--out n127.QC
```
n127.QC.sexcheck containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.
15 brains have problem!
```r
library(data.table)
# Read in file
valid <- fread("n127.valid.sample")
dat <- fread("n127.QC.sexcheck")[FID%in%valid$FID]
fwrite(dat[,c("FID","IID")], "n127.QC.valid", sep="\t")  # still use all samples

```

### Sample overlap/relatedness
```bash
/home/data1/R/genotype/n2222/plink \
--bfile n127 \
--extract n127.QC.prune.in \
--keep n127.QC.valid \
--rel-cutoff 0.125 \
--out n127.QC
```


## Generate final QC'ed target data file
```bash
/home/data1/R/genotype/n2222/plink \
--bfile n127 \
--make-bed \
--keep n127.QC.rel.id \
--out n127.QC \
--extract n127.QC.snplist \
--exclude n127.mismatch \
--a1-allele n127.a1
```
--a1-allele: 5552868 assignments made.
5552868 variants and 125 people pass filters and QC.




### Update Effect Size
When the effect size relates to disease risk and is thus given as an odds ratio (OR), rather than BETA (for continuous traits), then the PRS is computed as a product of ORs. To simplify this calculation, we take the natural logarithm of the OR so that the PRS can be computed using summation instead (which can be back-transformed afterwards). We can obtain the transformed summary statistics with R:
```r
AN <- fread("../pgcAN2/pgcAN2.QC.gz") %>%
    # And immediately change the alleles to upper cases
    .[,c("REF","ALT"):=list(toupper(REF), toupper(ALT))]
colnames(AN)[1:5]=c("CHR","BP","SNP","A2","A1") # change colnames for match n127

## update AN snps names as our data
AN$rs=AN$SNP
id=match(AN$rs,info$rs)
AN$SNP=info$SNP[id]

na.id=which(is.na(AN$SNP))
AN$SNP[na.id]=AN$rs[na.id]


AN$BETA <- log(AN$OR)  ## AN data already use BETA value

colnames(AN)[8]="P"
write.table(AN, "AN.QC.Transformed", quote=F, row.names=F)
```

### Clumping
```bash
/home/data1/R/genotype/n2222/plink \
    --bfile n127.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump AN.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out n127
```
1110562 more top variant IDs missing; see log file.
--clump: 175702 clumps formed from 5552866 top variants.
Results written to n127.clumped .

n127.clumped, containing the index SNPs after clumping is performed. We can extract the index SNP ID by performing the following command:
```bash
awk 'NR!=1{print $3}' n127.clumped >  n127.valid.snp
cat n127.valid.snp | wc -l
# 175704
```

## Generate PRS
plink provides a convenient function --score and --q-score-range for calculating polygenic scores.

We will need three files:

1. The base data file: AN.QC.Transformed
2. A file containing SNP IDs and their corresponding P-values ($3 because SNP ID is located in the third column; $8 because the P-value is located in the eighth column)
```bash
awk '{print $3,$8}' AN.QC.Transformed > SNP.pvalue
```
3. A file containing the different P-value thresholds for inclusion of SNPs in the PRS. Here calculate PRS corresponding to a few thresholds for illustration purposes:
```bash
echo "0.001 0 0.001" > range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
```

calculate PRS
```bash
/home/data1/R/genotype/n2222/plink \
    --bfile n127.QC \
    --score AN.QC.Transformed 3 5 6 header \
    --q-score-range range_list SNP.pvalue \
    --extract n127.valid.snp \
    --out n127
```
$3 is snpname, $5 is A1, $6 is BETA value, effect size

175702 variants and 125 people pass filters and QC.
Note: No phenotypes present.
Warning: 6487823 lines skipped in --score file (6487729 due to variant ID
mismatch, 94 due to allele code mismatch); see n127.nopred for details.
--score: 175608 valid predictors loaded.
Warning: 6487824 lines skipped in --q-score-range data file.
--score: 7 ranges processed.
Results written to n127.*.profile.

### Accounting for Population Stratification
calculate PCs
```bash
# First, we need to perform prunning
/home/data1/R/genotype/n2222/plink  \
    --bfile n127.QC \
    --indep-pairwise 200 50 0.25 \
    --out n127
# Then we calculate the first 20 PCs
/home/data1/R/genotype/n2222/plink  \
    --bfile n127.QC \
    --extract n127.prune.in \
    --pca 20 \
    --out n127
```

