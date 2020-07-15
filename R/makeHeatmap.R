
library("ggplot2")
library("reshape2")

d<-read.table("~/expression_atlas-homo_sapiens.tsv",quote="",sep="\t",header=F,stringsAsFactors = FALSE)
colnames(d1)<-d1[1,]
d<-d1[-c(1),]
d1<-as.data.frame(lapply(d[,-c(1)],as.numeric))
d1[,"Experiment"]<-as.factor(d[,"Experiment"])
d2<-melt(d1,variable.name="Expression")
h<-ggplot(data=d2,mapping=aes(x=Experiment,y=Expression,fill=value))+geom_tile()+xlab(label="Experiment")
