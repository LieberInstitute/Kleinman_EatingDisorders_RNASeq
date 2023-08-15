#volcano plots
#based on C:\Users\licha\Documents\R\LIBD_VA_PTSD_RNAseq_4Region-main\R_scripts\Plotting\Volcano\volcano_plots.R

library(SummarizedExperiment)
library(readxl)
library(jaffelab)
library(edgeR)
library(limma)
library(recount)
library(TeachingDemos) # for shadow text
library(clusterProfiler)
library(org.Hs.eg.db)
library(VariantAnnotation)
library(RColorBrewer)
library(biomaRt)
library(readxl)

setwd("/home/data1/R/ED/DE")

load("ED_contrast_gene_n122.Rdata")
# https://support.bioconductor.org/p/53177/

# contrast fit is
# cmtx <- makeContrasts( "MDD-Control", "ED-MDD","ED-Control", levels= modQsva)

# coef=1 is topTable for MDD-Control
mddct = topTable(eBayes(contrasts.fit(fitGene,cmtx)),number=nrow(rse_gene),coef=1)

# coef=2 is topTable for ED-mdd
edmdd = topTable(eBayes(contrasts.fit(fitGene,cmtx)),number=nrow(rse_gene),coef=2)

# coef=3 is topTable for ED-Control
edct = topTable(eBayes(contrasts.fit(fitGene,cmtx)),number=nrow(rse_gene),coef=3)


library(ggrepel)
# plot overall pvalues, got them from coef only equal to 2
ggplot(data=edmdd, aes(x=logFC, y=-log10(P.Value),  label=Symbol)) +
        geom_point(pch=21) +
        theme_minimal() +
        # geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-0.6, 0.6), col="red") +
        geom_hline(yintercept=-log10(0.05), col="red")


# https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

  BiocManager::install('EnhancedVolcano')

# library(airway)
library(magrittr)
library(EnhancedVolcano)
# The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-4.
p1=EnhancedVolcano(mddct,lab = mddct$Symbol,
  x = 'logFC',y = 'P.Value',title = 'MDD vs Control',pCutoff=1e-4)

p2=EnhancedVolcano(edmdd,lab = edmdd$Symbol,
  selectLab = c('AC027612.1','FXYD3','LMOD1','TRAF1','BEST1','TRNT1','RP1-95L4.4','ALDH1L1','STX11','SLCO4A1','RNF144B','SMIM3','SNORD115-14','NEDD9','SRGN','RP11-338N10.3','AC027612.1','FXYD3','LMOD1','TRAF1','BEST1','TRNT1','RP1-95L4.4','ALDH1L1','STX11','SLCO4A1','RNF144B','SMIM3','SNORD115-14','NEDD9','SRGN','RP11-338N10.3'),
  x = 'logFC',y = 'P.Value',title = 'ED vs MDD',pCutoff=1e-4)

p3=EnhancedVolcano(edct,lab = edct$Symbol,x = 'logFC',y = 'P.Value',title = 'ED vs Control',pCutoff=1e-4)


library(gridExtra)
jpeg(file='volcano.jpeg',res=300,width=12,height=8,unit='in')
grid.arrange(p1,p2,p3,ncol=3)
dev.off()
