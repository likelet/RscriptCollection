rm(list = ls())
gc()
#loading packages 
require(ggsci)
require(data.table)
require(cowplot)

#sort function from oncoprint complexHeatmap 
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


#read data 
met <- fread('scriptsWithData/tumor_vs_met/All_metastasis_filter_mut.maf')
tum <- fread('scriptsWithData/tumor_vs_met/All_primary_filter_mut.maf')

#retain mutations potentially functionally 
retaind_features=c("frameshift deletion","frameshift insertion","stopgain","nonsynonymous SNV","splicing","stoploss")
met=met[met$Variant_Classification %in% retaind_features,]
tum=tum[tum$Variant_Classification %in% retaind_features,]



#format data 

tum.VS.met <- data.table(tum$Tumor_Sample_Barcode, tum$Tumor_Sample_Barcode %in% met$Tumor_Sample_Barcode & tum$ChromChange %in% met$ChromChange, tum$Hugo_Symbol)
met.VS.tum <- data.table(met$Tumor_Sample_Barcode, met$Tumor_Sample_Barcode %in% tum$Tumor_Sample_Barcode & met$ChromChange %in% tum$ChromChange, met$Hugo_Symbol)
fuck <- unique(as.vector(rbind(met[,16], tum[,16])))
Initial <- apply(fuck, 1, function(x) table(tum.VS.met[V1 == x, -3, with = FALSE])[1])
Common <- apply(fuck, 1, function(x) unique(table(met.VS.tum[V1 == x, -3, with = FALSE])[2], table(tum.VS.met[V1 == x])[2]))
Recurrence <- apply(fuck, 1, function(x) table(met.VS.tum[V1 == x, -3, with = FALSE])[1])
mutations <- data.table(Initial = Initial, Common = Common, Recurrence = Recurrence, Samples = as.vector(t(fuck)))
rm(Initial, Common, Recurrence, fuck)
#data format for bar plot 
mutations <- melt(mutations, id = 'Samples', variable.name = 'Types')
mutations$Samples=factor(mutations$Samples,levels =unique(mutations$Samples) )

#data format for oncoplot 
met.VS.tum[, V2 := as.character(V2)][V2 == 'FALSE', V2 := "Recurrence"][V2 == 'TRUE', V2 := "Common"]
tum.VS.met[, V2 := as.character(V2)][V2 == 'FALSE', V2 := "Initial"][V2 == 'TRUE', V2 := "Common"]
mutations_onco <- rbind(met.VS.tum, tum.VS.met)

names(mutations_onco) <- c("Samples", "Types", "Genes")

# reorder genes and samples based on mutation numbers 
oncomatrix= data.frame(dcast(mutations_onco, Genes ~ Samples))
row.names(oncomatrix)=oncomatrix$Genes
oncomatrix=oncomatrix[,-1]
oncomatrix=memoSort(oncomatrix!=0)


mutations$Samples=factor(mutations$Samples,levels =unique(colnames(oncomatrix)) )

mbar <- ggplot() + 
  geom_col(data = mutations, aes(x = Samples, y = value, fill = Types)) + 
  ylab('Somantic mutations') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position="top")+
  scale_fill_lancet()


#set background 
backgroud=expand.grid(Samples=unique(mutations_onco$Samples),Genes=unique(mutations_onco$Genes))
backgroud$Samples=factor(backgroud$Samples,levels =colnames(oncomatrix) )
backgroud$Genes=factor(backgroud$Genes,levels =rev(row.names(oncomatrix) ))
mheatmap <- ggplot() + 
  geom_tile(data = backgroud, aes(x = Samples, y = Genes),fill = "Gray", width = 0.9, height = 0.9,  size = 1) +
  geom_tile(data = mutations_onco, aes(x = Samples, y = Genes, fill = Types), width = 0.9, height = 0.9,  size = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)) +
  theme(axis.text.y = element_text(face="bold.italic")) +
  theme(legend.position="none") +
  scale_fill_lancet(limits = c('Initial', 'Common', 'Recurrence'))

plot_grid(mbar,mheatmap,ncol = 1, align = 'v',rel_heights = c(1/4,3/4))
# ggdraw() + draw_plot(mbar, 0.1, 0.7, 1, 0.3) + draw_plot(mheatmap, 0, 0, 1, 0.7) + draw_plot_label(c("a", "b"), c(0, 0), c(1, .5), size = 15)
save_plot('met_VS_tum.png', ggplot2::last_plot(), base_width = 11, base_height = 9)