id=which(sigGeneCNT$`q_ED-MDD`<0.05)
id=sigGeneCNT[id,]$gencodeID
id=match(id,rownames(rowData(rse_gene)))

CORed$p[id,5]
#  ENSG00000232531.3 ENSG00000089356.17 ENSG00000163431.12 ENSG00000056558.10 TRAF1
#        0.094605982        0.043727376        0.127994406        0.003125535 
# ENSG00000167995.15 ENSG00000072756.16 TRNT1  ENSG00000217648.1 ENSG00000144908.13 
#        0.105303554        0.007855430        0.061725201        0.089103620 
#  ENSG00000135604.9 ENSG00000101187.15  ENSG00000137393.9  ENSG00000256235.1 
#        0.021541229        0.094409311        0.276700651        0.265060903 
#  ENSG00000199960.1 ENSG00000111859.16  ENSG00000122862.4  ENSG00000269978.1 
#        0.737881385        0.064811903        0.094068670        0.720695267 
#  ENSG00000257243.1  ENSG00000183888.4  ENSG00000090376.9 ENSG00000117151.12 
#        0.027910376        0.148946545        0.251302033        0.958205557 
#  ENSG00000262919.7 ENSG00000096070.19 ENSG00000076067.11 
#        0.245466868        0.286369492        0.023484883 

df=data.frame(edexp[,id],PRS=edonlyprs.res[,5])
colnames(df)=c(rowData(rse_gene)$Symbol[id],"PRS")



library(ggpubr)
library(enrichplot)
sp1= ggscatter(df, x = "TRAF1", y = "PRS",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson", label.x = 3, label.y = 5)
   
sp2= ggscatter(df, x = "TRNT1", y = "PRS",
   add = "reg.line",  # Add regressin line
   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson", label.x = 4, label.y = 5)
   
   
   
#############
# plot DE genes GO
goBP=enrichGO(gene          = sigem,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
d1=dotplot(goBP,showCategory = 20,title="ED-MDD p<0.005 Biological Process GO")
c1=cnetplot(goBP,foldChange=sigem,circular=TRUE,colorEdge=TRUE)
goBP2<- pairwise_termsim(goBP)
t1=treeplot(goBP2)
bt1=emapplot(goBP2)
pdf('EDvsCT less005 GO BP.pdf',height=20,width=20)
cowplot::plot_grid(d1,c1,t1,bt1,ncol=2)
dev.off()

goMF=enrichGO(gene          = sigec,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "MF",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)

d3=dotplot(goBP,showCategory = 20,title="MDD-CT p<0.005 Biological Process GO")
c3=cnetplot(goBP,foldChange=sigmc,circular=TRUE,colorEdge=TRUE)
goBP2<- pairwise_termsim(goBP)
t3=treeplot(goBP2)
bt3=emapplot(goBP2)
pdf('MDDvsCT less005 GO BP.pdf',height=20,width=20)
cowplot::plot_grid(d3,c3,t3,bt3,ncol=2)
dev.off()


goCC=enrichGO(gene          = sigmc,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.2,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
d4=dotplot(goCC,showCategory = 20,title="MDD-CT p<0.005 Cellular Component GO")			  
c4=cnetplot(goCC,foldChange=sigmc,circular=TRUE,colorEdge=TRUE)
goCC2<- pairwise_termsim(goCC)
t4=treeplot(goCC2)
bt4=emapplot(goCC2)
pdf('MDDvsCT less005 GO CC.pdf',height=20,width=20)
cowplot::plot_grid(d4,c4,t4,bt4,ncol=2)
dev.off()


