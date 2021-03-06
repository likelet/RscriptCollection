library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
Input1="Input1_50.bed"
Input2="Input2_50.bed"
JK1="JK-1_L1_A006.bed"
JK2="JK-2_L1_A012.bed"
JK3="JK-3_L1_A002.bed"
JK4="JK-4_L3_A004.bed"
JK5igg="JK-5-IgG_L2_A006.bed"
JK5="JK-5_L3_A005.bed"
JK6igg="JK-6-IgG_L2_A012.bed"
JK6="JK-6_L2_A007.bed"
sample1="Sample1_50.bed"
sample2="Sample2_50.bed"
BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=1e-3
extend=300
shift=0
ws=100

#read data bed or bam file 

Input1_set=MEDIPS.createSet(file = Input1, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
Input2_set=MEDIPS.createSet(file = Input2, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
JK1_set=MEDIPS.createSet(file = JK1, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
JK2_set=MEDIPS.createSet(file = JK2, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
JK3_set=MEDIPS.createSet(file = JK3, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
JK4_set=MEDIPS.createSet(file = JK4, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
JK5igg_set=MEDIPS.createSet(file = JK5igg, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
JK5_set=MEDIPS.createSet(file = JK5, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
JK6igg_set=MEDIPS.createSet(file = JK6igg, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
JK6_set=MEDIPS.createSet(file = JK6, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws)
sample1_set=MEDIPS.createSet(file = sample1, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 
sample2_set=MEDIPS.createSet(file = sample2, BSgenome = BSgenome,extend = extend, shift = shift, uniq = uniq, window_size = ws) 

#correlation analysis
cor.matrix = MEDIPS.correlation(MSets = c(Input1_set,Input2_set,sample1_set,sample2_set,JK1_set,JK2_set,JK3_set,JK4_set,JK5igg_set,JK5_set,JK6igg_set,JK6_set), plot = F, method = "pearson")

#GC enrichment anlysis
JK1_er=MEDIPS.CpGenrich(file =JK1, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK2_er=MEDIPS.CpGenrich(file =JK2, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK3_er=MEDIPS.CpGenrich(file =JK3, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK4_er=MEDIPS.CpGenrich(file =JK4, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK5igg_er=MEDIPS.CpGenrich(file =JK5igg, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK5_er=MEDIPS.CpGenrich(file =JK5, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK6igg_er=MEDIPS.CpGenrich(file =JK6igg, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)
JK6_er=MEDIPS.CpGenrich(file =JK6, BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq)

samplenames=c("JK1","JK2","JK3","JK4","JK5igg","JK5","JK6igg","JK6")
enrichment.score.relH_vector=c(JK1_er$enrichment.score.relH,JK2_er$enrichment.score.relH,JK3_er$enrichment.score.relH,JK4_er$enrichment.score.relH,JK5igg_er$enrichment.score.relH,JK5_er$enrichment.score.relH,JK6igg_er$enrichment.score.relH,JK6_er$enrichment.score.relH)
enrichment.score.GoGe_vector=c(JK1_er$enrichment.score.GoGe,JK2_er$enrichment.score.GoGe,JK3_er$enrichment.score.GoGe,JK4_er$enrichment.score.GoGe,JK5igg_er$enrichment.score.GoGe,JK5_er$enrichment.score.GoGe,JK6igg_er$enrichment.score.GoGe,JK6_er$enrichment.score.GoGe)
CpgEnrichTable=data.frame(sample=samplenames,enrichment.score.relH=enrichment.score.relH_vector,enrichment.score.GoGe=enrichment.score.GoGe_vector)

#saturation analysis 

JK6_sr = MEDIPS.saturation(file = JK6, BSgenome = BSgenome,
                       uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

JK1_sr = MEDIPS.saturation(file = JK1, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

JK2_sr = MEDIPS.saturation(file = JK2, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

JK3_sr = MEDIPS.saturation(file = JK3, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

JK4_sr = MEDIPS.saturation(file = JK4, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)
JK5_sr = MEDIPS.saturation(file = JK5, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

inpu1_sr = MEDIPS.saturation(file = Input1, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)
inpu2_sr = MEDIPS.saturation(file = Input2, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)
sample1_sr = MEDIPS.saturation(file = sample1, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)
sample2_sr = MEDIPS.saturation(file = sample2, BSgenome = BSgenome,
                           uniq = uniq, extend = extend, shift = shift, window_size = ws, nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)

jpeg(file = "input1.jpeg")
MEDIPS.plotSaturation(inpu1_sr)
dev.off()
jpeg(file = "input1.jpeg")
MEDIPS.plotSaturation(inpu2_sr)
dev.off()

jpeg(file = "sample1.jpeg")
MEDIPS.plotSaturation(sample1_sr)
dev.off()

jpeg(file = "sample2.jpeg")
MEDIPS.plotSaturation(sample2_sr)
dev.off()


#DMR analysis using replicates
treatset=c(JK2_set,JK4_set,JK6_set)
koset=c(JK1_set,JK3_set,JK5_set)
CS=MEDIPS.couplingVector(pattern = "CG", refObj = treatset[[1]])
mr.edgeR = MEDIPS.meth(MSet1 = koset, MSet2 = treatset,
                       CSet = CS, ISet1 = JK5igg_set, ISet2 = JK6igg_set,
                       p.adj = "BH", diff.method = "edgeR", MeDIP = T, CNV = F, minRowSum = 6)

mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = 0.1, adj = T,
                               ratio = NULL, bg.counts = NULL, CNV = F)
