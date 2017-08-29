#created by Yu Sun and modified by both Sun and me
# For plotting crispr cleavage plot at genome level 

setwd("~/Rproject/crisprplot")
library(data.table)
# library(gtools)
library(ggsci)
library(cowplot)
score <- fread("s4g-0407_result")
position <- read.table("genome.fa.fai", skip = 1)
fuck <-
  data.frame(do.call('rbind', strsplit(
    as.character(score$Position), ':', fixed = TRUE
  )))
score <- data.table(fuck, score[, 9])
rm(fuck)
position <- data.table(position[, c(1:3)])
position$V3 <- position$V3 - 16915
score[, 2] <-
  apply(score, 1, function(x)
    x[2] <- as.numeric(x[2]) + position[position$V1 == x[1], V3])
names(score) <- c('chr', 'loc', 'score')
score$chr <-
  factor(score$chr, levels = mixedsort(levels(score$chr)))
mut <- fread('s4g_hg19_multianno.txt')
mut <- mut[Func.ensGene == 'exonic', c(1, 2, 7), with = FALSE]
mut[, 2] <-
  apply(mut, 1, function(x)
    x[2] <-
      as.numeric(x[2]) + position[position$V1 == paste('chr', x[1], sep = ""), V3])
mut <-
  merge(
    score,
    mut,
    by.x = 'loc',
    by.y = 'Start',
    all.x = TRUE,
    all.y = TRUE
  )
mut$score=log2(mut$score)
colorlist=scale_fill_uchicago()
colorvec=alpha(rep(colorlist$palette(4),6), 0.4)
ggplot() + geom_col(data = mut,
                    aes(x = loc, y = score, group = chr),
                    width = 5e6) + geom_text(data = mut, aes(
                      x = loc,
                      y = (score + 1),
                      label = Gene.ensGene
                    ),size=rel(2.8)) + scale_y_continuous(limit=c(0,15),expand = c(0,0))+ylab("Log(Clevage Score)")+xlab("Chromesome")+
  scale_x_continuous(expand = c(0,0),breaks = position$V3 + position$V2/2, labels = position$V1) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 0.6))+geom_rect(data = position, aes(xmin = position$V3, xmax = V2+V3), 
                                                                      ymin = -Inf,
                                                                      ymax = Inf,
                                                                      fill = colorvec)


save_plot("cleavage_plot.pdf",plot = ggplot2::last_plot(),base_aspect_ratio=14/3)


