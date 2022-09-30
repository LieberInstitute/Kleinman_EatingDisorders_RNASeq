load("prs_an.pgc_2019.rda")
load("/home/data1/R/ED/DE/ED_contrast_gene_n122.Rdata")

which(!(colData(rse_gene)$BrNum %in% score$FID))
#[1] 27
colData(rse_gene)$BrNum[27]
#[1] "Br1732"
# R17347, Br1732, MDD

#Br1732 in n2108 version

##### find vGene$E and prs correlations

# align vGene$E with prs score
exp=vGene$E
which(colnames(exp)=="R17347")
#[1] 27
exp=exp[,-27]

id=match(colData(rse_gene)$BrNum, score$FID)
identical(score$FID[id][-27],colData(rse_gene)$BrNum[-27])
#[1] TRUE
edprs=score[id,][-27,]

identical(edprs$FID,colData(rse_gene)$BrNum[-27])
# [1] TRUE
identical(colnames(exp),colData(rse_gene)$RNum[-27])
# [1] TRUE
save(edprs,exp,file='/home/data1/R/ED/prs_score/n121_Exp_prs.rda')

## pearson correlation
## methods in ER & AJ's 2018 paper
# https://github.com/LieberInstitute/Code_for_WGCNA_MS_ER
library(psych)
edprs$Group=colData(rse_gene)$Group[-27]
edprs$group=sapply(edprs$Group,unclass)

# test code
data=data.frame(gene=exp[1,],prs1=edprs[,3],age=edprs[,15],group=edprs$group)
result=corr.test(x=data,y=NULL,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
# x could be a matrix/dataframe of expression
# y could be a matrix/dataframe of prs
# result[[1]] is the correlation matrix
# result[[4]] is the unadjusted probability values, result$p
# result[[5]] is the adjusted probability values, or result$p.adj

edprs1=edprs[,c(3:12,15,18:37,39)]
exp=t(exp)
# COR=corr.test(x=exp,y=edprs1,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)

# make adjustment for polygenic risk score
edprs.res=sapply(edprs[,3:12],function(x){
  fit=lm(x~Batch+Age+Group+PC1+PC2+PC3+PC4+PC5,data=edprs)
  return(residuals(fit))
})
COR=corr.test(x=exp,y=edprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)

save(COR,exp,edprs,edprs.res,file='/home/data1/R/ED/prs_score/pearsonCor_result.rda')

###### linear model selection #########
#summary(lm(p0.001~Age+Group+Dx+Batch,data=edprs))

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                        0.0426098  0.2734510   0.156   0.8765  
#Age                               -0.0009503  0.0037374  -0.254   0.7998  
#GroupMDD                          -0.0104235  0.3823724  -0.027   0.9783  
#GroupED                           -0.1620933  0.2241936  -0.723   0.4713  
#DxControl                                 NA         NA      NA       NA  
#DxED                               0.0672170  0.2852353   0.236   0.8142  
#DxMDD                             -0.0881937  0.3265371  -0.270   0.7876  
#DxOCD                             -1.6782526  0.6879044  -2.440   0.0164 *
# DxPTSD                             0.5482956  0.3610016   1.519   0.1318  
#BatchIllumina_HumanOmni2-5-8-v1-1  0.0515878  0.1813235   0.285   0.7766  
#BatchIllumina_HumanOmni2-5-8-v1-3  0.2856246  0.1787510   1.598   0.1130  
#BatchIllumina_HumanOmni5-Quad     -0.2686691  0.2725626  -0.986   0.3265  
#BatchMacrogen_11062018            -0.2538390  0.2788527  -0.910   0.3647  
#Batchtao_h650                      0.0690397  0.3043030   0.227   0.8210  
#BatchTao_omni2.5-8_v1.5_56         0.3437933  0.6777697   0.507   0.6130  
#BatchTao_omni2.5v1.4_136           0.7696296  0.3681509   2.091   0.0390 *
######################  
  
  edprs.res=sapply(edprs[,3:12],function(x){
    fit=lm(x~Batch+Age+Group+PC1+PC2+PC3+PC4+PC6+PC7+PC8+PC13,data=edprs)
    return(residuals(fit))
  })
## PCs were selected one by one through linear model: summary(lm(PCx~Group+Batch+Age,data=edprs))
## if the PCx has effect on ED/MDD, it will be chosen into PRS adjusted linear model
COR=corr.test(x=exp,y=edprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
# adjust method could change to 'fdr'
min(COR$p)
# [1] 2.963512e-06 at COR$p[1139,8]
length(which(COR$p[,8]<1e-2))
# [1] 619
length(which(COR$p[,8]<0.005))
# [1] 253

############################## 
## run enrichment analysis ###
##############################
####
library(org.Hs.eg.db)
library(clusterProfiler)
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]

id005=which(COR$p[,8]<0.005)
sig005=rowData(rse_gene)[id005,]
sig005=unique(as.character(sig005$EntrezID[!is.na(sig005$EntrezID)]))
length(sig005)
# [1] 152
# sig01 n423
moduleGeneList =list(sig005=sig005)

goBP <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff  = 1, qvalueCutoff  = 1,
                       readable= TRUE)
goMF <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "MF", pAdjustMethod = "BH",
                       pvalueCutoff  = 1, qvalueCutoff  = 1,
                       readable= TRUE)
goCC <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "CC", pAdjustMethod = "BH",
                       pvalueCutoff  = 1, qvalueCutoff  = 1,
                       readable= TRUE)
kegg <- compareCluster(moduleGeneList, fun = "enrichKEGG",
                       universe = geneUniverse,  pAdjustMethod = "BH",
                       pvalueCutoff  = 1, qvalueCutoff  = 1)
write.csv(goBP$compareClusterResult,file='result.csv')


############# Group Separately ##################
## eating disorder only  ###
edexp=vGene$E[,colData(rse_gene)$Group=='ED']
edexp=t(edexp)
edonlyprs.res=sapply(edprs[,3:12],function(x){
  fit=lm(x~Batch+Age+PC1+PC2+PC3+PC4+PC6+PC7+PC8+PC13,data=edprs,subset=(edprs$Group=='ED'))
  return(residuals(fit))
})
CORed=corr.test(x=edexp,y=edonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)

save(CORed,edexp,edonlyprs.res,file='n38_edonly_pearsonCor_result_adj.rda')

which(CORed$p==min(CORed$p))%%25148
#[1] 4146
CORed$p.adj[4146,]
#p0.001     p0.01     p0.05      p0.1      p0.2      p0.5        p1    p1e-04 
#1.0000000 1.0000000 1.0000000 1.0000000 0.2739007 1.0000000 1.0000000 1.0000000 
#p1e-06    p5e-08 
#1.0000000 1.0000000
id005=which(CORed$p[,5]<0.001)
sig005=rowData(rse_gene)[id005,]
ed001=unique(as.character(sig005$EntrezID[!is.na(sig005$EntrezID)]))
length(ed001)
# [1] 695
moduleGeneList =list(ed=ed001)
goBPed <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "BP", pAdjustMethod = "BH",
                       pvalueCutoff  = 0.1, qvalueCutoff  = 0.1,
                       readable= TRUE)
write.csv(goBPed@compareClusterResult,file='goBP_edonly_less001_result.csv')
save(CORed,edexp,edonlyprs.res,goBPed,file='n38_edonly_pearsonCor_result_adj.rda')

## termsim plot 20220920##
library(ggpubr)
library(enrichplot)
library(ggnewscale)
gobp=enrichGO(gene          = ed001,
              universe      = geneUniverse,
              OrgDb         = org.Hs.eg.db,
              ont           = "BP",
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.1,
              qvalueCutoff  = 0.5,
              readable      = TRUE)
dotplot(gobp,showCategory = 20,title="PRS Biological Process GO")
gobp2<- pairwise_termsim(gobp)
treeplot(gobp2)
jpeg(file='prs_score/prs_GO_EDonly_BP.jpg',height=8,width=12,res=300,units='in')
treeplot(gobp2)
dev.off()


# plot
goMFed <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "MF", pAdjustMethod = "BH",
                       pvalueCutoff  = .1, qvalueCutoff  = .5,
                       readable= TRUE)
goCCed <- compareCluster(moduleGeneList, fun = "enrichGO",
                       universe = geneUniverse, OrgDb = org.Hs.eg.db,
                       ont = "CC", pAdjustMethod = "BH",
                       pvalueCutoff  = .1, qvalueCutoff  = .5,
                       readable= TRUE)
kegged <- compareCluster(moduleGeneList, fun = "enrichKEGG",
                       universe = geneUniverse,  pAdjustMethod = "BH",
                       pvalueCutoff  = .1, qvalueCutoff  = .5)

edo <- enrichDGN(ed001)
pdf('barplot_edonly.pdf',h=12,w=6)
barplot(edo, showCategory=20)
dev.off()

## new version of DOSE plotting, GO over-representation analysis
gene=ed001
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]
egobp <- enrichGO(gene          = gene,
                universe      = geneUniverse,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
p1=barplot(egobp,showCategory = 20,title='ED only Biological Process GO')

egomf <- enrichGO(gene          = gene,
                  universe      = geneUniverse,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
p2=barplot(egomf,showCategory = 20,title='ED only Molecular Function GO')

egocc <- enrichGO(gene          = gene,
                  universe      = geneUniverse,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
p3=barplot(egocc,showCategory = 20,title='ED only Cellular Component GO')

# dotplot
pd1=dotplot(egobp,showCategory = 20,title='ED only Biological Process GO')
pd2=dotplot(egomf,showCategory = 20,title='ED only Molecular Function GO')
pd3=dotplot(egocc,showCategory = 20,title='ED only Cellular Component GO')
pdf('ED only GO dotplot.pdf',height=12,width=24)
cowplot::plot_grid(pd1, pd2, pd3, ncol=3, labels=LETTERS[1:3])
dev.off()

# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# Visualize enriched GO terms as a directed acyclic graph
pdf('ED only GOplot.pdf',height=21,width=10)
p4=goplot(egobp,title="Goplot of Biological Process")
p5=goplot(egomf,title="Goplot of Molecular Function")
p6=goplot(egocc,title="Goplot of Cellular Component")
cowplot::plot_grid(p4, p5, p6, ncol=1,nrow=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
dev.off() # title no show

# Network plot of enriched terms
edox <- setReadable(egobp, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=gene)
p2 <- cnetplot(edox, foldChange = gene,categorySize='pvalue')
p3 <- cnetplot(edox, foldChange = gene,circular=TRUE,colorEdge=TRUE)

pb <- cnetplot(edox, foldChange = gene,circular=TRUE,colorEdge=TRUE)
edoxm <- setReadable(egomf, 'org.Hs.eg.db', 'ENTREZID')
pm <- cnetplot(edoxm,foldChange = gene, circular=TRUE, colorEdge=TRUE)
edoxc <- setReadable(egocc, 'org.Hs.eg.db', 'ENTREZID')
pc <- cnetplot(edoxc,foldChange = gene, circular=TRUE, colorEdge=TRUE)

pdf('gene-concept network a_bp b_mf c_cc.pdf',height=12,width=27)
cowplot::plot_grid(pb, pm, pc, ncol=3, labels=LETTERS[1:3]) #, rel_widths=c(.8, .8, 1.2))
dev.off()

# Network plot with label subset of the nodes
p1 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p2 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p3 <- cnetplot(edox, node_label="all") 
p4 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
pdf('gene-concept network with labels.pdf',height=20,width=20)
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
dev.off()


# Heatmap-like functional classification
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
pdf('heatmap-like classification.pdf',height=8,width=12)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])
dev.off()
pb <- heatplot(edox, foldChange=geneList, showCategory=5)
pm <- heatplot(edoxm, foldChange=geneList, showCategory=5)
pc <- heatplot(edoxc, foldChange=geneList, showCategory=5)

pdf('heatmap-like classification a_bp b_mf c_cc.pdf',height=12,width=16)
cowplot::plot_grid(pb, pm, pc, ncol=1, labels=LETTERS[1:3])
dev.off()

# Tree plot
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
pdf('tree plot.pdf',height=10,width=20)
aplot::plot_list(p1, p2, tag_levels='A')
dev.off()

tb<- treeplot(pairwise_termsim(edox))
tm<- treeplot(pairwise_termsim(edoxm))
tc<- treeplot(pairwise_termsim(edoxc))
pdf('tree plot a_bp b_mf c_cc.pdf',height=10,width=30)
aplot::plot_list(tb, tm, tc, tag_levels='A')
dev.off()


# Biological theme comparison
xx <- pairwise_termsim(egobp)                     
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2) 
p3 <- emapplot(xx, pie="count")
p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
pdf('biological theme comparison.pdf',height=20,width=20)
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])
dev.off()


# upset plot for overexpressed genes
p1=upsetplot(egobp)
p2=upsetplot(egomf)
p3=upsetplot(egocc)
pdf('upset plot for overlapped gene.pdf',height=15,width=12)
cowplot::plot_grid(p1, p2, p3, ncol=1, labels=LETTERS[1:3])
dev.off()



# KEGG
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

# KEGG module over-representation analysis
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)   

# Visualize enriched KEGG pathways
browseKEGG(kk, 'hsa04110')





## control only  ###
ctexp=vGene$E[,colData(rse_gene)$Group=='Control']
ctexp=t(ctexp)
ctonlyprs.res=sapply(edprs[,3:12],function(x){
  fit=lm(x~Batch+Age+PC1+PC2+PC3+PC4+PC6+PC7+PC8+PC13,data=edprs,subset=(edprs$Group=='Control'))
  return(residuals(fit))
})
CORct=corr.test(x=ctexp,y=ctonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
CORct$p[18831,]
#p0.001 (min)       p0.01        p0.05         p0.1         p0.2         p0.5
#3.083648e-06 6.169527e-02 4.298933e-01 2.472294e-01 1.350932e-01 2.067740e-01
#p1       p1e-04       p1e-06       p5e-08
#1.794241e-01 1.008645e-01 1.925969e-01 4.353412e-01
id005=which(CORct$p[,5]<0.001)
sig005=rowData(rse_gene)[id005,]
sig005=unique(as.character(sig005$EntrezID[!is.na(sig005$EntrezID)]))
length(sig005)
# [1] 6
moduleGeneList =list(sig005=sig005)
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]
goBPct <- compareCluster(moduleGeneList, fun = "enrichGO",
                         universe = geneUniverse, OrgDb = org.Hs.eg.db,
                         ont = "BP", pAdjustMethod = "BH",
                         pvalueCutoff  = 1, qvalueCutoff  = 1,
                         readable= TRUE)
write.csv(goBPct@compareClusterResult,file='goBP_ctonly_less001_result.csv')
save(CORct,ctexp,ctonlyprs.res,goBPct,file='n41_ctonly_pearsonCor_result_adj.rda')

### MDD only ###
mdexp=vGene$E[,colData(rse_gene)$Group=='MDD']
mdexp=t(mdexp)
which(rownames(mdexp)=='R17347')
#[1] 5
mdexp=mdexp[-5,]
mdonlyprs.res=sapply(edprs[,3:12],function(x){
  fit=lm(x~Batch+Age+PC1+PC2+PC3+PC4+PC6+PC7+PC8+PC13,data=edprs,subset=(edprs$Group=='MDD'))
  return(residuals(fit))
})
CORmd=corr.test(x=mdexp,y=mdonlyprs.res,use='pairwise',method='pearson',adjust='holm',alpha=0.5,ci=TRUE)
CORmd$p[19501,]
#p0.001        p0.01        p0.05         p0.1         p0.2         p0.5
#1.299864e-05 7.137147e-03 3.258415e-02 6.680986e-02 7.991816e-02 1.387328e-01
#p1       p1e-04       p1e-06       p5e-08
#1.223651e-01 1.725789e-04 2.612781e-02 5.663367e-01
id005=which(CORmd$p[,5]<0.001)
sig005=rowData(rse_gene)[id005,]
sig005=unique(as.character(sig005$EntrezID[!is.na(sig005$EntrezID)]))
length(sig005) #4

moduleGeneList =list(sig005=sig005)
geneUniverse = as.character(sigGeneOV$EntrezID)
geneUniverse = geneUniverse[!is.na(geneUniverse)]
goBPmd <- compareCluster(moduleGeneList, fun = "enrichGO",
                         universe = geneUniverse, OrgDb = org.Hs.eg.db,
                         ont = "BP", pAdjustMethod = "BH",
                         pvalueCutoff  = 1, qvalueCutoff  = 1,
                         readable= TRUE)
write.csv(goBPmd@compareClusterResult,file='goBP_mddonly_less001_result.csv')
save(CORmd,mdexp,mdonlyprs.res,goBPmd,file='n42_mddonly_pearsonCor_result_adj.rda')
