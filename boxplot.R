setwd("/users/cli/ED/gene/")

library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(devtools)
library(limma)
library(edgeR)
library(ggplot2)

## boxplot vGene
load("/users/cli/ED/gene/ED_contrast_gene.Rdata")

focus=sigGeneCNT[order(sigGeneCNT$'q_ED-MDD'),]
mod2 = model.matrix(~Group + AgeDeath + overallMapRate , data=colData(rse_gene) )
mod2=cbind(mod2,qSVs[,1:3])
y=vGene$E['ENSG00000008282.8',]
y_clean_p2 <- cleaningY(t(y), mod2, P = 2)

df=data.frame(t(y_clean_p2),Dx=mod[,2],Group=rse_gene$Group)
dif=mean(y_clean_p2-y)

ggplot(data=df,aes(x=Group,fill=Group,y=df[,1]-dif))+
geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 0.7, aes(colour = Group), position = position_jitterdodge()) + 
 		labs(x = "Group", y = "Adjusted RPKM",title="LMOD1\np(ED-MDD) = 5.70e-05") +
 		scale_colour_brewer(palette = "Set1") +
 		scale_fill_brewer(palette = "Set1") +
 		theme(axis.title.x=element_blank(),
 		legend.background = element_rect(colour = "black"),
 		legend.title = element_blank(),
 		legend.position='none')




ids=which(sigGeneCNT$'q_ED-MDD'<0.05)
box_EM=list()

for(k in 1:length(ids)){
	i=ids[k]
	y=vGene$E[rownames(sigGeneCNT)[i],]
	y_clean_p2 <- cleaningY(t(y), mod2, P = 3)
	df=data.frame(t(y_clean_p2),Group=rse_gene$Group)
	dif=mean(y_clean_p2-y)
	p=ggplot(data=df,aes(x=Group,fill=Group,y=df[,1]-dif))+
	geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 0.7, aes(colour = Group), position = position_jitterdodge()) + 
	 		labs(x = "Group", y = "Adjusted RPKM",title=paste(sigGeneCNT$Symbol[i],signif(sigGeneCNT$'p_ED-MDD'[i], digits=3),sep='\n')) +
	 		scale_colour_brewer(palette = "Set1") +
	 		scale_fill_brewer(palette = "Set1") +
	 		theme(axis.title.x=element_blank(),
	 		legend.background = element_rect(colour = "black"),
	 		legend.title = element_blank(),
	 		legend.position='none')
	jpeg(file=paste0(sigGeneCNT$Symbol[i],'.jpeg'),height=5,width=5,res=300,unit='in')
	print(p)
	dev.off()
}
		
		
## ED vs. CT genes ###
genes=read.delim(pipe('pbpaste'),header=FALSE)
ids=which(sigGeneCNT$Symbol %in% paste(genes$V1))


for(k in 1:length(ids)){
	i=ids[k]
	y=vGene$E[rownames(sigGeneCNT)[i],]
	y_clean_p2 <- cleaningY(t(y), mod2, P = 3)
	df=data.frame(t(y_clean_p2),Group=rse_gene$Group)
	dif=mean(y_clean_p2-y)
		
	p=ggplot(data=df,aes(x=Group,fill=Group,y=df[,1]-dif))+
	geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 0.7, aes(colour = Group), position = position_jitterdodge()) + 
	 		labs(x = "Group", y = "Adjusted RPKM",title=paste(sigGeneCNT$Symbol[i],signif(sigGeneCNT$'p_ED-Control'[i], digits=3),sep='\n')) +
	 		scale_colour_brewer(palette = "Set1") +
	 		scale_fill_brewer(palette = "Set1") +
	 		theme(axis.title.x=element_blank(),
	 		legend.background = element_rect(colour = "black"),
	 		legend.title = element_blank(),
	 		legend.position='none')
	jpeg(file=paste0(sigGeneCNT$Symbol[i],'.jpeg'),height=5,width=5,res=300,unit='in')
	print(p)
	dev.off()

}


		
########################
##### example from nicotine project ######
y_clean_p2 <- cleaningY(y, mod, P = 2)
df=data.frame(t(y_clean_p2),Smoking=mod[,2],dosage=rse_gene_dos$dosage)
df$Smoking=gsub("1","Exposed",paste(df$Smoking))
df$Smoking=gsub("0","Unexposed",paste(df$Smoking))
df$Smoking=factor(df$Smoking,c('Unexposed','Exposed'))
df$dosage=factor(df$dosage,c('Non-Smoker','Light Smoker','Heavy Smoker'))
mean(y_clean_p2[2,]-y[2,])
mean(y_clean_p2[1,]-y[1,])

 MACRO=ggplot(data=df,aes(x=Smoking,fill=Smoking, y= df[,2]+5.58) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 1, aes(colour = Smoking), position = position_jitterdodge()) + 
 		labs(x = "Adult group", y = "Adjusted RPKM",title="MACRO") +
 		scale_colour_brewer(palette = "Set1") +
 		scale_fill_brewer(palette = "Set1") +
 		theme(axis.title.x=element_blank(),
 		legend.background = element_rect(colour = "black"),
 		legend.title = element_blank(),
 		legend.position='none') 
 MACRO_dos=ggplot(data=df,aes(x=dosage,fill=dosage, y= df[,2]+5.58) ) +geom_boxplot(outlier.colour = NA, alpha = 0) + geom_point(alpha = 1, aes(colour = dosage), position = position_jitterdodge()) + 
	    labs(x = "Adult group",y = "Adjusted RPKM",title="MACRO") +
	    scale_colour_brewer(palette = "Set1") +
	    scale_fill_brewer(palette = "Set1") +
	    theme(axis.title.x=element_blank(),
	    legend.background = element_rect(colour = "black"),
	    legend.title = element_blank(),
	    legend.position='none')
		
interaction_boxplot <- arrangeGrob(MACRO,MACRO_dos, 
	ncol=2,nrow=1,widths=c(2,3))
ggsave(interaction_boxplot, 
	file="~/Documents/Nicotine/Adult_sva_boxplot_MACRO.pdf", height=5,width=8) 
 
