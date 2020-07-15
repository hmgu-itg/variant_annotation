library("ggplot2")
library("reshape2")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==1) {
  args[2] = "out.svg"
}

d<-read.table(args[1],quote="",sep="\t",header=T,stringsAsFactors = FALSE)
names(d)[c(1)]<-c("Experiment")
d[,"Experiment"]<-as.factor(d[,"Experiment"])
d2<-melt(d,variable.name="Expression")
h<-ggplot(data=d2,mapping=aes(x=Experiment,y=Expression,fill=value))+geom_tile()+xlab(label="Experiment")
ggsave(file=args[2],plot=h,width=7,height=12)
