library(survival) 
library(ggplot2) 
library(survminer)
library(RColorBrewer)

df = read.csv("KRAS_OS.txt", header = T, sep="\t")

head(df)
target = df$Group



surv_object_os <- Surv(time=df$OS_follow, event = df$OS)
fit_os <- survfit(surv_object_os ~ target, data=df)
a2<-ggsurvplot(fit_os, pval = T, risk.table = T,
               palette= c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"),
               legend="none",title="ARID1A Mutation",
               xlab="Overall survival (months)",
               ylab="Probability")
a2$plot = a2$plot+
  theme(axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.line=element_line(size=1))
a2
