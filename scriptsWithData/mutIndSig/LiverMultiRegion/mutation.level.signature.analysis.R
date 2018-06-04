suppressMessages(library(Palimpsest))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
ref_genome <- BSgenome.Hsapiens.UCSC.hg19

#set input parameter
cancertype="liver"
sampleID="313544"

#load ref data 
load("Palimpsest-master/data/ensgene_hg19.RData")
load("Palimpsest-master/data/cytoband_hg19.RData")
# all cancer.signature 
cancer.signature = list(`adrenocortical carcinoma` = c(1,2,4,5,6,13,18), `ALL` = c(1,2,5,13), `AML` = c(1,5), `bladder` = c(1,2,5,10,13), `breast` = c(1,2,3,5,6,8,10,13,17,18,20,26), `cervix` = c(1,2,5,6,10,13,26), `chondrosarcoma` = c(1,5), `CLL` = c(1,2,5,9,13), `colorectum` = c(1,5,6,10), `glioblastoma` = c(1,5,11), `glioma low grade` = c(1,5,6,14), `head and neck` = c(1,2,4,5,7,13), `kidney chromophobe` = c(1,5,6), `kidney clear cell` = c(1,5,6,27), `kidney papillary` = c(1,2,5,13), `liver` = c(1,4,5,6,12,16,17,22,23,24), `lung adeno` = c(1,2,4,5,6,13,17), `lung small cell` = c(1,4,5,15), `lung squamous` = c(1,2,4,5,13), `lymphoma B-cell` = c(1,2,5,9,13,17), `lymphoma hodgkin` = c(1,5,25), `medulloblastma` = c(1,5,8), `melanoma` = c(1,5,7,11,17), `myeloma` = c(1,2,5,13),`nospharyngeal carcinoma` = c(1,2,5,6,13), `neuroblastoma` = c(1,5,18), `oesophagus` = c(1,2,4,5,6,13,17), `oral gingivo-bucca squamous` = c(1,2,5,7,13,29), `osteosarcoma` = c(1,2,5,6,13,30),`ovary` = c(1,3,5), `pancreas` = c(1,2,3,5,6,13), `paraganglioma` = c(1,5),`pilocytic astrocytoma` = c(1,5,19), `prostate` = c(1,5,6), `stomach` = c(1,2,5,13,15,17,18,20,21,26,28), `thyroid` = c(1,2,5,13), `urothelial carcinoma` = c(1,2,5,13,22), `uterine carcinoma` = c(1,2,5,6,10,13,14,26), `uterine carcinosarcoma` = c(1,2,5,6,10,13), `uvearl melanoma` = c(1,5,6), `small cell oesophagus carcinoma` = c(1,5,13))
save(cancer.signature,file = "cancer.signature.cosmic.RData")
names(cancer.signature)

my.sig<-paste("Signature.",cancer.signature[[cancertype]],sep="")


#result DIR
resdir <- "Results";if(!file.exists(resdir)) dir.create(resdir)

#read data 
mut_data<-read.delim(paste("data/",sampleID,"/",sampleID,".mut.txt",sep=""))

#set head name 
Sample.col = "Sample"
CHROM.col = "CHROM"
POS.col = "POS"
REF.col = "REF"
ALT.col = "ALT"

# Annotating the mutation data with necessary information
chroms <- unique(mut_data[, CHROM.col])
if (1 %in% chroms == TRUE) 
  mut_data[, CHROM.col] <- paste("chr", mut_data[, CHROM.col], sep = "")
"%ni%" <- Negate("%in%")
remove_chrs <- c("chrM", "chrMT")
mut_data <- mut_data[which(mut_data[, CHROM.col] %ni% remove_chrs), ]
mut_data$strand.mut <- "+"
mut_data$strand.gene <- NA
vcf <- palimpsest_dfPosXSegm(mut_data, dfPos.chrom.col = CHROM.col, 
                             dfPos.pos.col = POS.col, ensgene, dfSegm.chrom.col = "Chromosome.Name", 
                             dfSegm.start.col = "Gene.Start..bp.", dfSegm.end.col = "Gene.End..bp.", 
                             colsToAdd = "Associated.Gene.Name", namesColsToAdd = "Associated.Gene.Name")

vcf$strand.gene <- ensgene[match(vcf$Associated.Gene.Name, ensgene$Associated.Gene.Name), "Strand"]
vcf$strand.gene[vcf$strand.gene=="1"] = "+"
vcf$strand.gene[vcf$strand.gene=="-1"] = "-"
print("Adding mutation categories:")
vcf.snv <- palimpsest_addMutationContextToVcf(vcf[which(vcf$Type == "SNV"), ], ref_genome, chrom.col = CHROM.col, start.col = POS.col, 
                                              end.col = POS.col, ref.col = REF.col, alt.col = ALT.col, 
                                              strand.mut.col = "strand.mut", strand.gene.col = "strand.gene", 
                                              sample.col = Sample.col)
#vcf <- merge(vcf, vcf.snv, by = setdiff(intersect(names(vcf), names(vcf.snv)), "strand.mut"), all = TRUE, sort = FALSE)
vcf = vcf.snv
vcf$strand.mut <- vcf$strand.mut.x
#ind <- which(vcf$Type == "SNP")
#vcf[ind, "strand.mut"] <- vcf[ind, "strand.mut.y"]
vcf$strand.mut <- vcf$strand.mut.y
vcf <- vcf[, setdiff(names(vcf), c("strand.mut.x", "strand.mut.y"))]

# Processing input for mutational signature extraction
propMutsByCat <- palimpsestInput(vcf = vcf,
                                 type = "SNV",
                                 sample.col = "Sample",
                                 mutcat.col = "mutcat3",
                                 proportion = TRUE)
# Calculating contributions (exposures) of COSMIC signatures operative in NPC 
# Selecting relevant signatures

cos_signature_for_certain_cancer <- COSMIC_Signatures[my.sig,]
# Calculating contributions (exposures) of known NPC signatures
mycol <- c(1:nrow(cos_signature_for_certain_cancer))
signatures_exp <- deconvolution_fit(vcf=vcf,
                                    type = "SNV",
                                    input_data = propMutsByCat,
                                    threshold = 5,
                                    input_signatures = cos_signature_for_certain_cancer,
                                    sig_cols = mycol,
                                    plot = T,
                                    resdir = resdir)
# Plotting the exposures of signatures across the series:
deconvolution_exposure(mutSign_nums = signatures_exp$sig_nums,
                       mutSign_props = signatures_exp$sig_props,
                       sig_cols = mycol)

# Assigning the most likely signature at the origin of each mutation
vcf <- palimpsestOrigin(vcf,
                        type = "SNV",
                        sample.col = "Sample",
                        mutcat.col = "mutcat3",
                        signature_contribution = signatures_exp$sig_nums,
                        input_signatures = cos_signature_for_certain_cancer)

write.table(vcf, paste(resdir,"/",sampleID, "_All_mut_add_signature.maf", sep = ""), sep = "\t", quote = F, row.names = F)
