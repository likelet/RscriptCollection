rm(list = ls())
gc()
require(data.table)
require(cowplot)
df <- fread('stat_final.txt')
fuck <- data.table(table(df[,1]), cumsum(table(df[,1])))
fuck$V2=as.numeric(fuck$V2)
setnames(fuck, c('Occurrency', 'gene number', 'cumsum number'))
fuck <- data.frame(melt(fuck, id = 1, variable.name = "Type"))
fuck$Occurrency=factor(fuck$Occurrency,levels=names(table(df$V1)))
ggplot() +
  geom_bar(data = fuck, aes(x = Occurrency, y = value, fill = Type), position = 'identity', stat = 'identity', alpha = 0.3) + 
  scale_y_log10(expand = c(0.01, 0)) +
  scale_fill_manual(values = c('black', '#666666')) + 
  labs(y = 'Gene number')