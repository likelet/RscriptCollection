# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")

library("org.Hs.eg.db")
library(data.table)

# setwd("scriptsWithData/convertEnsemblToEntrez/")

fuck <- fread("mRNA.symbol.csv")
xx <- as.list(org.Hs.egENSEMBL2EG)
setnames(fuck, 1:3, c("symbol", "ensembl", "entrez"))

fuck[, entrez := unlist(xx)[ensembl] ]

fwrite(fuck, 'mRNA.entrez.csv')