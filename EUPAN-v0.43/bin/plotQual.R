#!/usr/bin/Rscript
argv <- commandArgs(TRUE)
library("ggplot2")
library("reshape2")
d=read.table(argv[1],header=T)
d1=melt(d)
plot1=paste(argv[1],"_heatmap.pdf",sep="")
d1$variable=factor(d1$variable,levels=rev(levels(d1$variable)))
k=ggplot(data=d1)+geom_tile(aes(x=Sample,y=variable,fill=value))+
  scale_x_discrete(name="")+scale_y_discrete(name="")+  
  scale_fill_gradient2(name="",low="red",mid="gold",high="darkgreen",
                       breaks=c(-1,0,1),labels=c("FAIL","WARNING","PASS"))+theme(axis.text.x=element_text(angle=90))
ggsave(plot1,k)

