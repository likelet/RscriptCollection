#!/usr/bin/env Rscript
message("Checking packages")
packageStartupMessage("Initializing ...", appendLF = FALSE)
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(ShortRead))
packageStartupMessage(" done")


parser <- OptionParser()
parser <- add_option(parser, c("--min-armlength"), action="store", default=4, 
                     help="An integer giving the minimum length of the arms of the palindromes to search for [default %default]")
parser <- add_option(parser, c("--max-looplength"), action="store", default=1, 
                     help="An integer giving the maximum length of \"the loop\" (i.e the sequence separating the 2 arms) of the palindromes to search for [default %default]. 
                     Note that by default (max.looplength=1), findPalindromes will search for strict palindromes only")
parser <- add_option(parser, c("--min-looplength"), action="store", default=0, 
                     help="An integer giving the minimum length of \"the loop\" of the palindromes to search for [default %default]")
parser <- add_option(parser, c("--max-mismatch"), action="store", default=0, 
                     help="The maximum number of mismatching letters allowed between the 2 arms of the palindromes to search for [default %default]")

rfq <- readFastq('.', pattern = '*.fastq|*.fq')

for (i in 1:length(rfq)) {
  print(findPalindromes(sread(rfq)[[i]], 
                        min.armlength = parse_args(parser)$`min-armlength`, 
                        max.looplength = parse_args(parser)$`max-looplength`,
                        min.looplength = parse_args(parser)$`min-looplength`,
                        max.mismatch = parse_args(parser)$`max-mismatch`))
}
