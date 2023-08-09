load("./DE/n127_ED_bmi_2qsvsPC_3snpPCs.rda")

ED_MDD <- sigGeneCNT$Symbol[which(sigGeneCNT$'q_ED-MDD'<0.05)]
ED_CT <- sigGeneCNT$Symbol[which(sigGeneCNT$'q_ED-Control'<0.05)]
MDD_CT <- sigGeneCNT$Symbol[which(sigGeneCNT$'q_MDD-Control'<0.05)]


### use library ggvenn
library(ggvenn)
list_venn <- list(`ED-MDD`=ED_MDD,`ED-CT`=ED_CT,`MDD-CT`=MDD_CT)
ggvenn(list_venn,c("ED-MDD","ED-CT","MDD-CT"),fill_color = c("red",'yellow','blue'),show_percentage = FALSE)
