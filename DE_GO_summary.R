# summary DE GO for eating disorder
# 
# load("DE/ED_contrast_gene_n122.Rdata")
library(httpgd)
library(org.Hs.eg.db)
library(clusterProfiler)
library(writexl)

id=which(sigGeneCNT$`p_ED-Control`<0.01)
sigem=as.character(sigGeneCNT$EntrezID[id])
sigem=sigem[!is.na(sigem)]
length(sigem)

goBP=enrichGO(gene          = sigem,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
			 

goMF=enrichGO(gene          = sigem,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)

goCC=enrichGO(gene          = sigem,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)


write_xlsx(list("BP"=goBP@result,"MF"=goMF@result,"CC"=goCC@result), 
				'DE/GO/01_ED-CT.xlsx')
