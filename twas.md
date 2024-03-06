### TWAS Analysis Explore
References:
http://gusevlab.org/projects/fusion/#compute-your-own-predictive-models
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4767558/

### Installation
```bash
# fusion software
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
cd fusion_twas-master

# 1000 Genomes
wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
tar xjvf LDREF.tar.bz2

# R packages
install.packages(c('optparse','RColorBrewer'))
library(devtools)
install_github("gabraham/plink2R/plink2R",repos = NULL)

# for computing own weights (AN GWAS), additional packages
# 1. Add the bundled GCTA binary gcta_nr_robust to path (coded by Po-Ru Loh for robust non-linear optimization)
# 2. Download and install PLINK2, add plink to path
# 3. Launch R and install the following required libraries:
install.packages(c('glmnet','methods'))
```

#### Input: GWAS summary statistics
The primary input is genome-wide summary statistics in LD-score format. At minimum, this is a flat file with a header row containing the following fields:

SNP – SNP identifier (rsID)
A1 – first allele (effect allele)
A2 – second allele (other allele)
Z – Z-scores, sign with respect to A1.

use AN paper GWAS results, the SNP id should matching to the LDREF/*.bim (1000G Ref, rs#), only the matched SNPs will be used for prediction.
AN paper GWAS result: /home/data1/R/prs/pgcAN2/pgcAN2.QC.gz
<strong> Calculate Z-score from beta and SE
https://bios25328.hakyimlab.org/post/2022/04/13/calculate-z-score-p-value-chi2-stat-from-gwas/
Z=beta/SE
Verify P value, only got absolute value of P, and small difference with AN original P values:
Z=pnorm(-abs(Z)) * 2
```r
library(data.table)
AN=fread("/home/data1/R/prs/ed_n127/AN.QC.Transformed")

AN$Z=AN$BETA/AN$SE
AN.twas=AN[,c(3,5,4,16,15)]
write.table(AN.twas, "/home/data1/R/twas/an/AN.twas", quote=F, row.names=F, col.names=T)
```

## 1. Build gene-level PLINK .bed files (build_bims_gene.R)
```r
library("SummarizedExperiment")
library("jaffelab")
library("data.table")
library("sessioninfo")
library("getopt")
library("BiocParallel")
library("tidyr")
library("here")
library(sva)


## load rse
load("/home/data1/R/ED/expr_cutoff/rse_gene_modified.Rdata")
rse <- rse_gene
assays(rse)$raw_expr <- assays(rse_gene)$RPKM # no rpkm assay in the rse
geneRpkm = recount::getRPKM(rse_gene, length="Length")
assays(rse)$raw_expr <- geneRpkm

## model
rse$Dx = factor(rse$Dx, c("Control","MDD","ED") )
with(colData(rse), table(Dx))

mod = model.matrix(~Dx + AgeDeath + mitoRate + rRNA_rate + totalAssignedGene + RIN + overallMapRate, 
		data = colData(rse))

## load degradation data to calculate expression PCs
load("/home/data1/R/ED/count_data/degradation_rse_EatingDisorders_caudate.Rdata")
cov_rse = cov_rse_caudate

## get qSVs for top bonferroni
qsvBonf = prcomp(t(log2(assays(cov_rse)$counts+1)))

## qsva, gene expression PCs
k = num.sv(log2(assays(cov_rse)$counts+1), mod)	# 12
qSVs = qsvBonf$x[,1:k]
getPcaVars(qsvBonf)[1:k]

pcs <- qSVs

pd = colData(rse)
colnames(pd)
colData(rse) <- cbind(colData(rse), pcs)

mod = model.matrix(~Dx + AgeDeath + overallMapRate +BMI+ PC1+PC2 + snpPC1+snpPC2+snpPC3, 
		data = colData(rse))
## Regress out effects. 
assays(rse)$clean_expr = cleaningY(log2(assays(rse)$raw_expr + 1), mod, P = 2)

colnames(rse) <- colData(rse)$BrNum

## load bim file
bim_file="/home/data1/R/prs/ed_n127/n127.QC"
bim <- fread(
    paste0(bim_file, ".bim"),
    col.names = c("chr", "snp", "position", "basepair", "allele1", "allele2")
)

# convert 23 to X, as is std in plink
bim$chr <- as.character(bim$chr)
# bim[chr == "23", ]$chr <- "X"  # I don't have 23

bim_gr <- GRanges(
    paste0("chr", bim$chr),
    IRanges(bim$basepair, width = 1)
)
mcols(bim_gr) <- bim[, -c("chr", "basepair")]

## Based on http://gusevlab.org/projects/fusion/#computing-your-own-functional-weights
## they restrict to 500kb on each side
rse_window <- resize(rowRanges(rse), width(rowRanges(rse)) + 500000 * 2, fix = "center")
mcols(rse_window) <- NULL

## Keep only those feature windows with some SNPs nearby
keep_feat <- which(countOverlaps(rse_window, bim_gr) > 0)
rse <- rse[keep_feat, ]
#  number of genes per chr
rse_window <- rse_window[keep_feat, ]
table(seqnames(rse))
#  chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9 chr10 chr11 chr12 chr13 
#  2344  1691  1443   957  1262  1270  1218   931   943  1014  1290  1300   506 
# chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22  chrX  chrY  chrM 
#   930   980  1168  1435   431  1493   629   263   584     0     0     0 
save(rse, file = "/home/data1/R/twas/ed/subsetted_rse.RData")

## Simplify RSE file to reduce mem
assays(rse) <- list("clean_expr" = assays(rse)$clean_expr)
mcols(rowRanges(rse)) <- NULL

## Subset to features
setwd("/home/data1/R/twas/ed")
dir.create("snp_files", showWarnings = FALSE)
dir.create("bim_files", showWarnings = FALSE)

## For checking if the triplet of bim files exists
check_bim <- function(x) {
    all(file.exists(paste0(x, c(".fam", ".bed", ".bim"))))
}

## For matching, due to complicated BrNums
## this function was used in filter_data/filter_snps.R
brnumerical <- function(x) {
    as.integer(gsub("Br|_.*", "", x))
}

## Create the bim files
i <- seq_len(nrow(rse))
i.names <- rownames(rse)
save(i, i.names, file = "i_info.RData")

write.table(data.frame(i, i.names), file = "input_ids.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

bim_files <- bpmapply(function(i, feat_id, clean = TRUE) {
    if (i == 1 || i %% 1000 == 0) {
        message("*******************************************************************************")
        message(paste(Sys.time(), "building bim file for i =", i, "corresponding to feature", feat_id))
    }

    base <- paste0("gene_", i)
    filt_snp <- paste0("filtered_snps_", base, ".txt")

    # change this file
    filt_bim <- gsub(".txt", "", filt_snp)
    filt_snp <- file.path("snp_files", filt_snp)
    dir.create(file.path("bim_files", base), showWarnings = FALSE)
    filt_bim <- file.path("bim_files", base, filt_bim)

    ## Re-use the bim file if it exists already
    if (check_bim(filt_bim)) {
        return(TRUE)
    }

    j <- subjectHits(findOverlaps(rse_window[i], bim_gr))

    fwrite(
        bim[j, "snp"], # %>% unique()
        file = filt_snp,
        sep = "\t", col.names = FALSE
    )

    system(paste(
        "/home/data1/R/genotype/n2222/plink --bfile", bim_file, "--extract", filt_snp,
        "--make-bed --out", filt_bim, "--memory 5000 --threads 1 --silent"
    ))

    ## Edit the "phenotype" column of the fam file
    filt_fam <- fread(paste0(filt_bim, ".fam"),
        col.names = c("famid", "w_famid", "w_famid_fa", "w_famid_mo", "sex_code", "phenotype")
    )

    ## Note BrNums might be duplicated, hence the use of match()
    m <- match(brnumerical(filt_fam$w_famid), brnumerical(colData(rse)$BrNum))
    stopifnot(all(!is.na(m)))

    ## Use cleaned expression for now. Could be an argument for the code.
    filt_fam$phenotype <- assays(rse)$clean_expr[i, m]

    ## Ovewrite fam file (for the phenotype info)
    fwrite(filt_fam, file = paste0(filt_bim, ".fam"), sep = " ", col.names = FALSE)
    ## Clean up
    #if (clean) unlink(filt_snp)

    return(check_bim(filt_bim))
}, i, i.names,  SIMPLIFY = FALSE)
```

## 2. Compute weights according to our dataset
```bash
FILELIST=$(echo "/home/data1/R/twas/ed/input_id10.txt")

while read -r NUM _; do
  parallel Rscript /home/data1/R/twas/scripts/FUSION.compute_weights.R \
   --bfile /home/data1/R/twas/ed/bim_files/gene_{}/filtered_snps_gene_{} \
   --tmp tmp_files/gene_{} \
   --out out_files/gene_{} \
   --PATH_plink /home/data1/R/genotype/n2222/plink \
   --PATH_gcta /home/ran/anaconda3/pkgs/gcta-1.94.1-h9ee0642_0/bin/gcta \
   --PATH_gemma /home/ran/anaconda3/pkgs/gemma-0.98.3-hb4ccc14_0/bin/gemma \
   --models top1,lasso --hsq_p 1.0001 --verbose 1 --save_hsq ::: $NUM
done < $FILELIST

# multiple cores
parallel -j 4
parallel -P 2
```

## 3. Merge Individual TWAS gene weights
```r
library("SummarizedExperiment")
library("sessioninfo")
library("getopt")
library("BiocParallel")

## Use the rse file from build_bims.R
rse_file <- file.path("subsetted_rse.RData") #

load("i_info.RData", verbose = TRUE)

## Compute TWAS weights in parallel
output_status <- mapply(function(i, feat_id, clean_bim = FALSE) {

    ## Check that the output file exists
    out_file <- file.path("out_files", paste0("gene_", i))

    ## Clean up
    if (clean_bim) unlink(dirname(filt_bim), recursive = TRUE)

    if(i%%1000==0){print(i)}
    return(file.exists(paste0(out_file, ".wgt.RDat")))
}, i, i.names,  SIMPLIFY = FALSE) 
output_status <- unlist(output_status)
names(output_status) <- paste0("gene_", i)


rdat_files <- dir("out_files", ".wgt.RDat", full.names = TRUE)
wglist <- paste0("gene",".list")

write.table(rdat_files, file = "gene.list", row.names = FALSE, col.names = FALSE, quote = FALSE)

system(paste0("Rscript /home/data1/R/twas/scz/fusion_twas-master/utils/FUSION.profile_wgt.R ", wglist, " > ","gene", ".profile.err 2>&1"))
# output file is gene.profile.err
# rename the gene.profile.err to wglist_summary.txt


# creating the .pos file
pos_match <- match(gsub("out_files/|\\.wgt\\.RDat", "", rdat_files),names(output_status))

pos_info <- data.frame(
    "WGT" = rdat_files,
    ## use the feature id
    "ID" = names(rowRanges(rse))[pos_match],
    "CHR" = gsub("chr", "", seqnames(rowRanges(rse)[pos_match])),
    "P0" = start(rowRanges(rse)[pos_match]),
    "P1" = end(rowRanges(rse)[pos_match]),
    "geneID" = mcols(rowRanges(rse))[, "gencodeID"][pos_match],
    stringsAsFactors = FALSE
)

sapply(pos_info, function(x) sum(is.na(x)))
pos_info$ID[pos_info$ID == ""] <- NA
pos_info$geneID[pos_info$geneID == ""] <- NA
write.table(pos_info[, -which(colnames(pos_info) == "geneID")],
    file = paste0("gene", ".pos"),
    row.names = FALSE, col.names = TRUE, quote = FALSE
)
save(pos_info, file = "pos_info.Rdata")
```

## 4. Apply Weights
```r
# change reference bim file rs# to what we have
library(data.table)
AN=fread('/home/data1/R/twas/an/an.twas.backup')

bim=read.table('/home/data1/brainseq_phase2/reference_hg38/LDREF_hg38/1000G.EUR.22.bim',sep='\t')

# get rs# only
bim$rs=sapply(strsplit(bim$V2, "\\:"), `[`, 1)
summary(bim$rs %in% AN$rs)
#    Mode   FALSE    TRUE 
# logical    1023   13255 
bim$SNP=AN$SNP[match(bim$rs,AN$rs)]
id=which(is.na(bim$SNP))
bim$SNP[id]=bim$rs[id]
bim=bim[,c(1,8,3:6)]
write.table(bim,'/home/data1/brainseq_phase2/reference_hg38/LDREF_hg38_libd/1000G.EUR.22.bim',sep='\t',
    quote=FALSE,col.names=FALSE,row.names=FALSE)
```
```bash
Rscript /home/data1/R/twas/scz/fusion_twas-master/FUSION.assoc_test.R \
    --sumstats /home/data1/R/twas/an/an.twas.backup \
    --weights /home/data1/R/twas/ed/gene.pos \
    --weights_dir /home/data1/R/twas/ed/ \
    --ref_ld_chr /home/data1/brainseq_phase2/reference_hg38/LDREF_hg38_libd/1000G.EUR. \
    --chr 22 \
    --out AN_gwas/chr22.dat
# same script for other chromosomes
```


## 5. combine twas into a single Rdata file
```r
library("SummarizedExperiment")
library("ggplot2")
library("gplots")
library("VennDiagram")
library("RColorBrewer")
library("readr")
library("here")
library("getopt")
library("GenomicRanges")
library("dplyr")
library("data.table")

datfiles=list.files("AN_gwas",pattern=".dat$")
dats=lapply(datfiles,function(x){
    dat=read.table(paste0("AN_gwas/",x),sep='\t',header=TRUE)
    return(dat)
})
names(dats)=gsub("\\.dat","",datfiles)

twas=do.call('rbind',dats)

load("subsetted_rse.RData")
# get rse ranges
rse_ranges <- rowRanges(rse) %>%
    ranges() %>%
    as.data.table()

# Set std. col names, specifically ID
colnames(rse_ranges) <- c("start", "end", "width", "ID")
rse_ranges=cbind(rse_ranges,rowData(rse)$Symbol,rowData(rse)$gene_type)

# convert to dt
twas_exp_dt <- as.data.table(twas)

# get ranges from rse
twas_exp_ranges <- merge(twas_exp_dt, rse_ranges, by = "ID")

# find midpoint for each gene
twas_mean_dist <- round(twas_exp_ranges[, c(width / 2) + start])

# add midpoint col
twas_exp_fin <- cbind(twas_exp_ranges, twas_mean_dist)

dir.create("rda", showWarnings = FALSE)
save(twas_exp_fin, file = "rda/twas_exp_ranges.Rdata")
```

## 6. TWAS plot generation (generate_twas_plot.R)
```r
library(ggplot2)
library(ggrepel)
library(dplyr)
library(data.table)
library(plotly)
library(htmlwidgets)
library(sessioninfo)
library(SummarizedExperiment)
library(ggpubr)
library(tools)

data.table::setDTthreads(threads = 1)
# Sourcing Data/Inst. Vars. ####
load("rda/twas_exp_ranges.Rdata")

dir.create(file.path("analysis/plots"),
           showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path("analysis/tables"),
           showWarnings = FALSE,
           recursive = TRUE)

# Filter N/A Z scores
twas_exp_fin$TWAS.Z=as.numeric(twas_exp_fin$TWAS.Z)
colnames(twas_exp_fin)[24:25]=c('genesymbol','genetype')
twas_z <- twas_exp_fin %>% filter(!is.na(TWAS.Z))

don <- list()

axisdf <- list()

don_key <- list()

p <- list()

intctv_plot <- list()

fin_plot <- list()

don[[1]] <-
        twas_z %>%
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len = max(end)) %>%
        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(twas_z,
                  .,
                  by = c("CHR" = "CHR")) %>%
        
        # Add a cumulative position of each SNP
        arrange(CHR, twas_mean_dist) %>%
        mutate(BPcum = twas_mean_dist + tot)

axisdf[[1]] = don[[1]] %>% group_by(CHR) %>% summarise(center = (max(BPcum) + min(BPcum)) / 2)

# Prepare text description for each SNP:
    don[[1]]$text <-
        paste0(
            "Gene Symbol: ",
            don[[1]]$genesymbol,
           # "\nENSEMBL Gene ID: ",
           # don[[1]]$geneid,
           # "\nBrain Subregion: ",
           # don[[1]]$region,
            "\nChromosome: ",
            don[[1]]$CHR,
            "\nStart Position: ",
            don[[1]]$start,
            "\nEnd Position: ",
            don[[1]]$end,
            "\nZ score: ",
            don[[1]]$TWAS.Z %>% round(2)
        )
    
    don_key[[1]] <-
        highlight_key(don[[1]], ~ genesymbol, group = "Gene Symbol")


# TWAS Z Manhattan Plot ####
pdf(file = "analysis/plots/AN_TWAS_ManhattanPlot.pdf")
# storing ggplot as an object3
twas_exp_fin$region="Caudate"
sig <- qnorm(1 - 0.025 / table(twas_exp_fin$region))

# Bonferroni Correction
    sig_bonf <- sig[[1]]
ggplot(don_key[[1]], aes(x = BPcum, y = TWAS.Z, text = text)) +
        
        ggtitle(paste0("Gene Windows of ", "Caudate" , " TWAS")) +
        # Show all points
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
        scale_color_manual(values = rep(c("#861657", "#D56AA0"), 22)) +
        geom_hline(
            yintercept = c(sig_bonf,-sig_bonf),
            color = "grey40",
            linetype = "dashed"
        ) +
        
        # custom X axis:
        scale_x_continuous(labels = axisdf[[1]]$CHR, breaks = axisdf[[1]]$center) +
        scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
        
        # Custom the theme:
        theme_bw() +
        theme(
            legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        )

dev.off()

# Scatter plots ####
pdf(
    'analysis/plots/AN_TWAS_ScatterPlots.pdf',
    useDingbats = FALSE,
    width = 10,
    height = 10
)

twas_z$region="Caudate"
twas_z_select <-
    select(twas_z, ID, genesymbol, TWAS.Z, TWAS.P, region) %>%
    as.data.table()

# FDR calculation
twas_z_select$fdr.p <-
    p.adjust(twas_z_select$TWAS.P, 'fdr')

twas_z$fdr.p<- p.adjust(twas_z$TWAS.P,'fdr')
# bonferroni correction
twas_z$bonf.p<- p.adjust(twas_z$TWAS.P,'bonferroni')

write.csv(twas_z,file='analysis/tables/AN_TWAS_20240112.csv')

save.image("rda/generate_plots_data.RData")
```

## 7. Enrichment tests on TWAS FDR <5% genes
```r
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(devtools)
library(xlsx)

data.table::setDTthreads(threads = 1)

## Define background universe of genes ####

load("rda/generate_plots_data.RData")

## Find Entrez IDs ####
# create a new EntrezID column in twas_z which contains the Entrez IDs  in the rse
twas_z[, EntrezID := rowData(rse)[match(twas_z$ID, rowData(rse)$gencodeID),]$EntrezID]
twas_z_fdr <- twas_z[fdr.p < 0.05,]

go_caud <- enrichGO(gene=twas_z_fdr$EntrezID, OrgDb="org.Hs.eg.db",
        keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 1,
         pAdjustMethod = "fdr", universe = twas_z$EntrezID,
         qvalueCutoff = 1, readable = TRUE)

## Save as XLSX ####
go_dt <- as.data.table(go_caud)
go_dt <- go_dt[order(qvalue),]
go_dt$GeneRatio = paste0(" ", go_dt$GeneRatio)
write.csv(go_dt, "analysis/tables/AN_EnrichmentTest.csv", 
        col.names=TRUE, row.names=TRUE)
```
