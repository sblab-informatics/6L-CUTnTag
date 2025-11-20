library("methylKit")
library("genomation")
library("GenomicRanges")

file.list <- list("E14CnT4X6LH3K4Me1B1.mC.cov","E14CnT4X6LH3K4Me1B2.mC.cov","E14CnT4X6LH3K4Me1B3.mC.cov")


myobj <- methRead(file.list,sample.id=list("H3K4Me1_rep1","H3K4Me1_rep2","H3K4Me1_rep3"),pipeline = "bismarkCoverage",assembly="mm10",treatment=c(1,1,1),mincov = 10)

meth=unite(myobj, destrand=FALSE)

png("H3K4Me1_Correlation_mC.pdf")
getCorrelation(meth_tiles,plot=TRUE)
dev.off()



file.list <- list("E14CnT4X6LH3K4Me3B1.mC.cov","E14CnT4X6LH3K4Me3B2.mC.cov","E14CnT4X6LH3K4Me3B3.mC.cov")


myobj <- methRead(file.list,sample.id=list("H3K4Me3_rep1","H3K4Me3_rep2","H3K4Me3_rep3"),pipeline = "bismarkCoverage",assembly="mm10",treatment=c(1,1,1),mincov = 10)

meth=unite(myobj, destrand=FALSE)

png("H3K4Me3_Correlation_mC.pdf")
getCorrelation(meth_tiles,plot=TRUE)
dev.off()


file.list <- list("E14CnT4X6LH3K27Me3B2.mC.cov","E14CnT4X6LH3K27Me3B3.mC.cov")


myobj <- methRead(file.list,sample.id=list("H3K27Me3_rep1","H3K27Me3_rep2"),pipeline = "bismarkCoverage",assembly="mm10",treatment=c(1,1),mincov = 10)

meth=unite(myobj, destrand=FALSE)

png("H3K27Me3_Correlation_mC.pdf")
getCorrelation(meth_tiles,plot=TRUE)
dev.off()


################################


file.list <- list("E14CnT4X6LH3K4Me1B1.hmC.cov","E14CnT4X6LH3K4Me1B2.hmC.cov","E14CnT4X6LH3K4Me1B3.hmC.cov")


myobj <- methRead(file.list,sample.id=list("H3K4Me1_rep1","H3K4Me1_rep2","H3K4Me1_rep3"),pipeline = "bismarkCoverage",assembly="mm10",treatment=c(1,1,1),mincov = 10)

meth=unite(myobj, destrand=FALSE)

png("H3K4Me1_Correlation_hmC.pdf")
getCorrelation(meth_tiles,plot=TRUE)
dev.off()



file.list <- list("E14CnT4X6LH3K4Me3B1.hmC.cov","E14CnT4X6LH3K4Me3B2.hmC.cov","E14CnT4X6LH3K4Me3B3.hmC.cov")


myobj <- methRead(file.list,sample.id=list("H3K4Me3_rep1","H3K4Me3_rep2","H3K4Me3_rep3"),pipeline = "bismarkCoverage",assembly="mm10",treatment=c(1,1,1),mincov = 10)

meth=unite(myobj, destrand=FALSE)

png("H3K4Me3_Correlation_hmC.pdf")
getCorrelation(meth_tiles,plot=TRUE)
dev.off()


file.list <- list("E14CnT4X6LH3K27Me3B2.hmC.cov","E14CnT4X6LH3K27Me3B3.hmC.cov")


myobj <- methRead(file.list,sample.id=list("H3K27Me3_rep1","H3K27Me3_rep2"),pipeline = "bismarkCoverage",assembly="mm10",treatment=c(1,1),mincov = 10)

meth=unite(myobj, destrand=FALSE)

png("H3K27Me3_Correlation_hmC.pdf")
getCorrelation(meth_tiles,plot=TRUE)
dev.off()
