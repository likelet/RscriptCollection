library(ggplot2)
args=commandArgs(trailingOnly=TRUE)
#read data from arguements
filestr=args[1]
print(filestr)
genename=strsplit(filestr,split = "\\.")[[1]]
mat <- read.table(filestr,header=T)
submat=mat[mat$Tissue_type=="Cancer",]

medianlist=tapply(submat$Exp_level,INDEX = submat$Cancer_type,median)

order=names(sort(medianlist))
allcancer =levels(mat$Cancer_type)
notorder=allcancer[!allcancer %in% order]
alllevel=c(notorder,order)

mat$Cancer_type<-factor(mat$Cancer_type,levels=order)
g=ggplot(mat, aes(Cancer_type, Exp_level)) + 
  geom_boxplot(aes(colour = Tissue_type)) + 
  theme_classic()+
  ggtitle(paste("Compare Analysis of ", genename,sep=""))+
  theme(axis.text.x = element_text(angle = -45,vjust = 0.6,hjust = 0.3))+
  # scale_y_continuous(limits = quantile(mat$Exp_level, c(0.1, 0.9)))+
  scale_colour_manual(values=c("#663333", "#003366"))
ggsave(g,filename =paste(filestr,".plot.pdf",sep=""),width = 10,height = 6 )
