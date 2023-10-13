library(ggplot2)
library(ggpubr)
library(ggrepel)

df = read.table("temp2.file", header=T, sep="\t")
head(df)


plot<-ggplot(df, aes(x=SignedP_mRNA, y=SignedP_type, fill=Type, label=Target)) +
  geom_point(size=2.5,shape=21, stroke=0.1, colour="white")+
  scale_x_continuous(limits = c(-6.5, 6.5), breaks = c(-10, -5, 0, 5, 10))+
  scale_y_continuous(limits = c(-6.5, 6.5), breaks = c(-10, -5, 0, 5, 10))+
  geom_text_repel(aes(label=Target), size=3.2, max.overlaps = Inf)+
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line("grey70"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust=0.5, face = "bold", size=20),
    axis.text.x = element_text(size=15, face = "bold"),
    axis.text.y = element_text(size=15, face = "bold"),
    axis.title.y = element_text(size=17),
    axis.title.x = element_text(size=17))+
  scale_fill_manual(values=c("#FFC041","#956ABF"))

plot

