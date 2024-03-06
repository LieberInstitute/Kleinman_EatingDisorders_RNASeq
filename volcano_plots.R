#volcano plots

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
library(writexl)

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

write_xlsx(list("ed_mdd"=edmdd,"ed_ct"=edct,"mdd_ct"=mddct),
          'n127_ED_cmtx_LFC.xlsx')

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
                   selectLab = c('MAFF','SERPINE1','LMOD1'),
                   x = 'logFC',y = 'P.Value',title = 'ED vs MDD',pCutoff=1e-4, 
                   labSize=3.5,legendPosition = 'bottom',
    legendLabSize = 10,
    legendIconSize = 4.0)
p3=EnhancedVolcano(edct,lab = edct$Symbol,x = 'logFC',y = 'P.Value',title = 'ED vs Control',pCutoff=1e-4)


library(gridExtra)
jpeg(file='volcano.jpeg',res=300,width=12,height=8,unit='in')
grid.arrange(p1,p2,p3,ncol=3)
dev.off()

#### 20230822 #### /home/data1/R/ED/DE/plots
png(file='ED-MDD volcano.png',height=10,width=6,res=300,unit="in")
EnhancedVolcano(edmdd,lab = edmdd$Symbol,
                selectLab = c('MAFF','SERPINE1','LMOD1'),
                x = 'logFC',y = 'P.Value',title = 'ED vs MDD',pCutoff=1e-4, 
                labSize=3.5,legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 4.0)
dev.off()

png(file='MDD Control volcano.png',height=10,width=6,res=300,unit="in")
EnhancedVolcano(mddct,lab = mddct$Symbol,
              #  selectLab = c('MAFF','SERPINE1','LMOD1'),
                x = 'logFC',y = 'P.Value',title = 'MDD vs Control',pCutoff=1e-4, 
                labSize=3.5,legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 4.0)
dev.off()

png(file='ED Control volcano.png',height=10,width=6,res=300,unit="in")
EnhancedVolcano(edct,lab = edct$Symbol,
              selectLab = c('DIO2','SERPINE1','HEMK1','SERPINA3','FAM110C'),
                x = 'logFC',y = 'P.Value',title = 'ED vs Control',pCutoff=1e-4, 
                labSize=3.5,legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 4.0)
dev.off()
