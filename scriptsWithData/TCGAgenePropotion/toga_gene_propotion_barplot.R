library(ggplot2)
library(ggsci)
library(ggsignif)


#Read data table ----
Gene1Gene2 <- data.frame(t(read.table("scriptsWithData/TCGAgenePropotion/me1etv4.txt",header=T,row.names=1,sep = " ")))
Gene1_mark = rep("Gene1_Low", nrow(Gene1Gene2))
Gene2_mark = rep("Gene2_Low", nrow(Gene1Gene2))
Gene1_mark[Gene1Gene2[, 1]>= median(Gene1Gene2[, 1])] = "Gene1_High"
Gene2_mark[Gene1Gene2[, 2]>= median(Gene1Gene2[, 2])] = "Gene2_High"
Gene1Gene2 = cbind(Gene1Gene2, Gene1_mark, Gene2_mark)

# statistic test for propotion  using chisq test 
mt = table(Gene1Gene2[, 3:4])
sig="NS"
Pvalue=chisq.test(mt)$p.value 
if(Pvalue<0.01){
  sig="**"
}else if(Pvalue<0.05){
  sig="*"
}

#plot function
g = 
  ggplot(Gene1Gene2[, 3:4], aes(x = Gene1_mark,fill = Gene2_mark)) + 
  geom_bar(position = "fill") +
  xlab('') + 
  ylab("Relative Fraction") +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = scales::percent)+
  theme_classic()+
  theme(legend.position="top") +
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text.x = element_text(size = 18, face="bold")) +
  theme(legend.text = element_text(size=18)) + 
  guides(fill = guide_legend(" ")) +
  scale_fill_lancet()+
  annotate("rect", xmin = 1, xmax = 2, ymin = 1.1, ymax = 1.1, alpha = 1, colour = "black")+
  geom_text(x = 1.5, y = 1.11, label = sig, size = 12)
ggsave("images/barRatioPlot.png", width = 8.6, height = 6.6)