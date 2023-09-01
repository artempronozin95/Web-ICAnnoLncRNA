library(seqinr)
library(LncFinder)

args <- commandArgs(TRUE)
strc <- args[2]
print(strc)
if (strc == "SS") {

frequencies <- readRDS(paste(args[1], 'frequencies.rds', sep=""))
model <- readRDS(paste(args[1], 'model.rds', sep=""))

seq <- read.fasta(args[3])
test <- run_RNAfold(seq, RNAfold.path = RNAfold.path, parallel.cores = 40)
result_2 <- lnc_finder(test, SS.features = TRUE, format = "SS",
frequencies.file = frequencies, svm.model = model,
parallel.cores = 2)
write.table(result_2, file=paste(args[4],'/lncFinder.csv',sep=""), col.names=F, sep = ',')
} else {

frequencies <- readRDS(paste(args[1], 'frequencies.rds', sep=""))
model <- readRDS(paste(args[1], 'model.rds', sep=""))

seq <- read.fasta(args[3])
result_2 <- lnc_finder(seq, SS.features = FALSE, format = "DNA",
frequencies.file = frequencies, svm.model = model,
parallel.cores = 2)
write.table(result_2, file=paste(args[4],'/lncFinder.csv',sep=""), col.names=F, sep = ',')
}
