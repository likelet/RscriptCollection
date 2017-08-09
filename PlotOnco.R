library("ggplot2")
library("cowplot")
library(ggsci)


memoSort <- function(M) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}
#preparefunctions

  data=read.csv("scriptsWithData/tumor_vs_met/plot_mat.csv",header=T,row.names=1)
  #plot barplot
  bardf=data$Barplot
  bardf=bardf/ncol(data)
  bardf=data.frame(Gene=row.names(data),value=bardf)
  percentage=rev(paste(round(100*bardf$value, 1), "%", sep=""))
 
  
  data=data[,-ncol(data)]
  data=memoSort(!is.na(data))
  facnames=colnames(data)
  data=data.frame(data)
  data$Gene=row.names(data)
  plotdf=melt(data,id="Gene")
  names(plotdf)[2]="Sample"
  plotdf=plotdf[!is.na(plotdf$value),]
  plotdf$Gene=factor(plotdf$Gene,levels=data$Gene)
  plotdf$Sample=factor(plotdf$Sample,levels=rev(facnames))
  
  backgroud=expand.grid(Sample=facnames,Gene=row.names(data))
  backgroud$Sample=factor(facnames,levels =facnames )
  backgroud$Gene=factor(backgroud$Gene,levels =rev(data$Gene))
  
  
  #barplot 
  bardf$Gene=factor(bardf$Gene,levels=rev(row.names(data)))
  barg=ggplot(bardf,aes(x=Gene,y=value))+geom_bar(fill=color,stat="identity")+ 
    scale_x_discrete(labels = c(percentage), position = "top")+
    theme(axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          axis.title=element_blank()
           )+
    coord_flip()
  
  
  plotdf=plotdf[plotdf$value,]
  #heatmap 
  heatmap=ggplot()+
    geom_tile(data = backgroud, aes(x = Sample, y = Gene),fill = "grey", width = 0.9, height = 0.9,  size = 1) +
    geom_tile(data = plotdf, aes(x = Sample, y = Gene),fill = color, width = 0.9, height = 0.9,  size = 1) + 
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(angle = 90,hjust = 0, vjust = 1)) +
    theme(axis.text.y = element_text(face="bold.italic",size = rel(0.7),hjust = 1)) +
    theme(axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank())
  
  plot_grid(heatmap,barg,nrow  = 1, align = "h",rel_widths = c(3,1))
  
  
  # save_plot(paste('',final.name,".pdf"), ggplot2::last_plot())


colorlist=scale_fill_npg()
colorvec=colorlist$palette(6)
# run functions 
g1=getHeatmap(filename="bl_scna.txt",color=colorvec[1],rmlable=F)
g2=getHeatmap(filename="bl_snv.txt",color=colorvec[2],rmlable=T)
g3=getHeatmap(filename="bl_sv.txt",color=colorvec[3],rmlable=T)
g4=getHeatmap(filename="pd_scna.txt",color=colorvec[4],rmlable=T)
g5=getHeatmap(filename="pd_snv.txt",color=colorvec[5],rmlable=T)
g6=getHeatmap(filename="mac.txt",color="black",rmlable=T)

#combine plot 
plot_grid(plot_grid(plotlist = g1,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g2,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g3,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g4,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g5,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g6,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          align = "v",
          ncol = 1,
          rel_heights = c(7.5/35,9/35,2.3/35,5.2/35,12/35,2.3/35)
)


save_plot("allfile.pdf", ggplot2::last_plot(),base_width = 6, base_height = 8)
