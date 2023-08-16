library(immunedeconv)

df <- read.table("exp.txt", sep="\t", row.names = 1, header = T)
head(df[1:5,1:5])

results<-immunedeconv::deconvolute(df, "estimate")
