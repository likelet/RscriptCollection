library("ggplot2")
library(ggsci)
library("ggsignif")


#Read data table ----
rawdata <- read.table("scriptsWithData/TCGAgenePropotion/me1etv4.txt",header=T,row.names=1,sep = " ")

dim(rawdata)
rawdata=data.frame(t(rawdata))


#get Sample size that bigger than median  for each gene
mark1=paste(names(rawdata)[1],"_low",sep = "")
mark2=paste(names(rawdata)[1],"_high",sep = "")
mark3=paste(names(rawdata)[2],"_low",sep = "")
mark4=paste(names(rawdata)[2],"_high",sep = "")

col1=rep(mark1,nrow(rawdata))
col2=rep(mark3,nrow(rawdata))
col1[rawdata[,1]>= median(rawdata[,1])]=mark2
col2[rawdata[,2]>= median(rawdata[,2])]=mark4
df=data.frame(col1,col2)

# statistic test for propotion  using chisq test 
mt=table(df)
sig="NS"
Pvalue=chisq.test(mt)$p.value 
if(Pvalue<0.01){
  sig="**"
}else if(Pvalue<0.05){
  sig="*"
}

#plot function
g=ggplot(df,aes(x = col1,fill = col2)) + 
  geom_bar(position = "fill")+
  xlab('') + 
  ylab("Raletive Fraction")+
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme_classic()+
  theme(legend.position="top")+
  scale_fill_lancet()+
  annotate("rect", xmin = 1, xmax = 2, ymin = 1.1, ymax =1.1, alpha=1,colour = "black")+
  geom_text(x = 1.5, y = 1.11, label = sig)
print(g)
  
save_plot('images/barRatioPlot.png',ggplot2::last_plot(), base_width = 6, base_height = 6)