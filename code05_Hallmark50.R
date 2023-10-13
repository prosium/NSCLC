library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(gridExtra)

df = read.table("input05_EstrogenIntensity.txt", header=T, sep="\t")
head(df)

df$Color <- factor(df$Color, levels=c("RNA:Protein", "RNA:NA", "NA:Protein", "NA:NA"))




results<-ggplot(df, aes(x=signedRNA, y=signedProtein, label=Label, fill=Color)) +
  geom_point(size=4,shape=21, stroke=1, colour="black")+
  scale_x_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4))+
  scale_y_continuous(limits = c(-4, 4), breaks = c(-4, -2, 0, 2, 4))+
  #geom_text_repel(aes(label=Label), size=3.2, max.overlaps = Inf)+
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line("grey70"),
    panel.grid.minor = element_blank(),
    #panel.border = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust=0.5, face = "bold", size=20),
    axis.text.x = element_text(size=15, face = "bold"),
    axis.text.y = element_text(size=15, face = "bold"),
    axis.title.y = element_text(size=17),
    axis.title.x = element_text(size=17))+
  scale_fill_manual(values=c("#F8766D","#7CAE00","#00BFC4","white"))

results

