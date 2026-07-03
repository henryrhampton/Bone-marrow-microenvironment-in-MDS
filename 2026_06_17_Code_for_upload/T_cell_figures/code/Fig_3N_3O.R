library(ggpubr)
library(ggplot2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/Users/henryhampton/OneDrive - UNSW/2026_06_17_Code_for_upload/T_cell_figures/code")


df<-read.csv("../data/CD8_GZMK_uniqueclone_anno.csv")
df$percentage<-100-df$scaled #scaled is uniqueclone, 100-scaled is non-unique
df$TP<-factor(df$TP,levels<-c("pre","post")) #plot order pre post

df$ID<-as.character(df$PID)

# Remove patient P11 as it had a small number of cells
df = df[df['PID'] != 'P11', ]
df = df[(df['sample'] != 'P18_C7D22'), ]
df <- df[!(df$sample == 'P18_C1D1' & df$group == 'NR'), ]


to_del = c("P12_C1D1", "P12_C7D1", "P12_C7D22") # there were few cells across timepoints for this patient

df <- df[!df$sample %in% to_del,  ]

options(repr.plot.width = 4, repr.plot.height = 3.5)
ggpaired(df[df$group=="R",], x="TP", y="percentage",
         fill="TP", palette=c("#3a7ab8","#115284"), id="ID",
         line.color="lightgray", line.size=0.15,
         point.size=1.5, shape=1,
         width=0.6, ylab="percent non-unique clonotype", xlab="timepoint") +
  geom_point(aes(x=TP, y=percentage),
             data=df[df$group=="R",],
             shape=21, size=1.5, fill="white", color="black")+
  stat_compare_means(label="p.format", method="t.test", paired=TRUE, label.x=1.8) +
  theme(axis.title=element_text(size=12, face="plain", color="black"),
        axis.text=element_text(size=11, face="plain", color="black"),
        legend.position="right") +
  ylim(0, 100)
ggsave("../results/CD8_cyto_R_percent_timepoint_clonaltype_0208.pdf",width=4,height = 3.5)


ggpaired(df[df$group=="NR",], x="TP", y="percentage",
         fill="TP", palette=c("#ffb347", "#fe9003"), id="ID",
         line.color="lightgray", line.size=0.15,
         point.size=0,
         width=0.6, ylab="percent non-unique clonotype", xlab="timepoint") +
  geom_point(aes(x=TP, y=percentage),
             data=df[df$group=="NR",],
             shape=21, size=1.5, fill="white", color="black")+
  stat_compare_means(label="p.format", method="t.test", paired=TRUE, label.x=1.8) +
  theme(axis.title=element_text(size=12, face="plain", color="black"),
        axis.text=element_text(size=11, face="plain", color="black"),
        legend.position="right") +
  ylim(0, 100)
ggsave("../results/CD8_cyto_R_percent_timepoint_clonaltype_0208.pdf",width=4,height = 3.5)
