library("ggplot2")
library("cowplot")
library(ggsci)
#preparefunctions
getHeatmap <- function(filename="scriptsWithData/heatmapForcfDNA/bl_scna.txt",color="red",rmlable=T){
  
 
  file=strsplit(filename,split = "/")[[1]][3]
  final.name=strsplit(file,split = "\\.")[[1]][1]
  data=read.csv(filename,header=T,row.names=1)
  #plot barplot
  bardf=data$Barplot
  bardf=bardf/ncol(data)
  bardf=data.frame(Gene=row.names(data),value=bardf)
  percentage=rev(paste(round(100*bardf$value, 1), "%", sep=""))
  bardf$Gene=factor(bardf$Gene,levels=rev(bardf$Gene))
  barg=ggplot(bardf,aes(x=Gene,y=value))+geom_bar(fill=color,stat="identity")+ 
    scale_x_discrete(labels = c(percentage))+
    theme(axis.ticks=element_blank(),
          axis.line=element_blank(),
          axis.title=element_blank(),
          axis.text.x = element_blank() )+
    coord_flip()
  
  data=data[,-18]
  facnames=colnames(data)
  data$Gene=row.names(data)
  plotdf=melt(data,id="Gene")
  names(plotdf)[2]="Sample"
  plotdf=plotdf[!is.na(plotdf$value),]
  plotdf$Gene=factor(plotdf$Gene,levels=data$Gene)
  plotdf$Sample=factor(plotdf$Sample,levels=rev(facnames))
  
  backgroud=expand.grid(Sample=facnames,Gene=row.names(data))
  backgroud$Sample=factor(facnames,levels =facnames )
  backgroud$Gene=factor(backgroud$Gene,levels =rev(data$Gene))
  heatmapg=ggplot()+
    geom_tile(data = backgroud, aes(x = Sample, y = Gene),fill = "grey", width = 0.9, height = 0.9,  size = 1) +
    geom_tile(data = plotdf, aes(x = Sample, y = Gene),fill = color, width = 0.9, height = 0.9,  size = 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 1,hjust = 1)) +
    theme(axis.text.y = element_text(face="bold.italic",size = rel(0.7))) +
    theme(axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y = element_text(hjust = 0))
    
  if(rmlable){
    heatmapg=heatmapg+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.line.x=element_blank())
  }  

  # save_plot(paste('scriptsWithData/heatmapForcfDNA/',final.name,".pdf"), ggplot2::last_plot())
  return(list(heatmapg,barg))
  
}

colorlist=scale_fill_npg()
colorvec=colorlist$palette(6)
# run functions 
g1=getHeatmap(filename="scriptsWithData/heatmapForcfDNA/bl_scna.txt",color=colorvec[1])
g2=getHeatmap(filename="scriptsWithData/heatmapForcfDNA/bl_snv.txt",color=colorvec[2])
g3=getHeatmap(filename="scriptsWithData/heatmapForcfDNA/bl_sv.txt",color=colorvec[3])
g4=getHeatmap(filename="scriptsWithData/heatmapForcfDNA/pd_scna.txt",color=colorvec[4])
g5=getHeatmap(filename="scriptsWithData/heatmapForcfDNA/pd_snv.txt",color=colorvec[5])
g6=getHeatmap(filename="scriptsWithData/heatmapForcfDNA/mac.txt",color="black",rmlable = T)

#combine plot 
plot_grid(plot_grid(plotlist = g1,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g2,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g3,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g4,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g5,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          plot_grid(plotlist = g6,nrow  = 1, align = "h",rel_widths = c(3/4,1/4)),
          align = "v",
          ncol = 1,
          rel_heights = c(5/35,9/35,2/35,5/35,14/35,2/35)
          )


save_plot("scriptsWithData/heatmapForcfDNA/allfile.pdf", ggplot2::last_plot(),base_width = 6, base_height = 19)
