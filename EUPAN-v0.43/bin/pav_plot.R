#!/usr/bin/Rscript
argv <- commandArgs(TRUE)



#PAV simulation result plotting
#By Yue Zhao and Zhiqiang Hu, 3/30/2016


#packages 
library(ggplot2)


#Function definition

#summarySE
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



#readSim: read simulation data
readSim <- function(file = NULL){
  DATA=read.table(file,sep="\t",header=T)
  d1=summarySE(data = DATA,measurevar ="Pan",groupvars = "Time")
  d1$Ave = d1$Pan
  colnames(d1)[3] = "class"
  d1$class = "pan"
  
  d2=summarySE(data = DATA,measurevar ="Core",groupvars = "Time")
  d2$Ave = d2$Core
  colnames(d2)[3] = "class"
  d2$class = "core"
  
  d = rbind(d1,d2)

  return(d)
}



#start processing data and plotting

all=readSim(file=argv[1])
plot1=paste(argv[1],"_PAV_plot.pdf",sep="")

p<-ggplot(all, aes(x=Time, y=Ave, color = class))
p=p+
geom_errorbar(aes(ymin=Ave-se, ymax=Ave+se), width=.05, size=1, color = "black")+
geom_line(size=1)+geom_point()+
scale_y_continuous(name="Average number")+
scale_x_continuous(name="line number")

ggsave(plot1,p)

