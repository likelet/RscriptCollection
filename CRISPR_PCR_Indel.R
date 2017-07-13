library(ggplot2)
library(ggjoy)
plotdatadf=read.csv("data/indel_crispr.csv",header=T)
#ggjoy
ggplot(plotdatadf, aes(x=pos, y=type, group=type,fill=type, height=..density..)) +
  geom_joy(scale=4) +
  scale_y_discrete(expand=c(0.01, 0)) +
  scale_x_continuous(expand=c(0, 0)) +
  theme_joy()+
  scale_fill_brewer(palette = 4) +
  theme(legend.position="none")+
  ylab("Indel Count")+xlab("Position of Target sequence")+
  ggtitle(paste("Indel distribution of target sequence in ",name,sep=""))
