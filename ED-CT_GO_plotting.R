## ED-CT GO plotting

load("DE/ED_contrast_gene_n122.Rdata")
library(httpgd)
library(org.Hs.eg.db)
library(clusterProfiler)

geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

# 52 & 73 genes
id1=which(sigGeneCNT$`p_ED-Control`<0.005 & sigGeneCNT$`p_MDD-Control`<0.005)
id2=which(sigGeneCNT$`p_ED-Control`<0.005 & sigGeneCNT$`p_MDD-Control`>0.005)

sigem1=as.character(sigGeneCNT$EntrezID[id1])
sigem2=as.character(sigGeneCNT$EntrezID[id2])
sigem1=sigem1[!is.na(sigem1)]
sigem2=sigem2[!is.na(sigem2)]

goBP1=enrichGO(gene          = sigem1,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)

goBP2=enrichGO(gene          = sigem2,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
# q value <0.05 no significant, switch to p value 0.005

goMF1=enrichGO(gene          = sigem1,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)

goMF2=enrichGO(gene          = sigem2,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
# no significant in MF
goCC1=enrichGO(gene          = sigem1,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)

goCC2=enrichGO(gene          = sigem2,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE) 
# BP and CC has significant
save(goBP1,goBP2,goMF1,goMF2,goCC1,goCC2,id1,id2,sigem1,sigem2,file='/home/data1/R/ED/DE/GO/GO_ED-CT-MDD_20220913.rda')

# plot
dbp=dotplot(goBP2,showCategory = 20,title="ED-CT only p<0.005 Biological Process GO")
dcc=dotplot(goCC2,showCategory = 20, title="ED-CT only p<0.005 Cellular Component GO")
cowplot::plot_grid(dbp,dcc,ncol=2)

jpeg(file='/home/data1/R/ED/DE/GO/GO_ED-CT-MDD_BP&CC.jpeg',height=8,width=16,res=300,unit='in')


# write result
write.csv(goBP2@result,file='/home/data1/R/ED/DE/GO/GO_ED-CT-MDD_BP.csv')
write.csv(goCC2@result,file='/home/data1/R/ED/DE/GO/GO_ED-CT-MDD_CC.csv')
