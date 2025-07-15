library(ggplot2)
dat<- read.table("mc_active.txt", fill=TRUE, header=T)
df = data.frame(dat)

require(reshape2)
sample_main <- melt(df)
levels(sample_main$variable) <- c("H3K4me1","H3K27ac")

pdf("mc_active.pdf")
p <- ggplot(sample_main, aes(variable, value)) + geom_violin((aes(fill=variable))) + scale_fill_manual(values = c("H3K4me1" = 'darkgreen',"H3K27ac" = 'darkblue'))
p + theme_classic()
dev.off()


library(ggplot2)
dat<- read.table("hmc_active.txt", fill=TRUE, header=T)
df = data.frame(dat)

require(reshape2)
sample_main <- melt(df)
levels(sample_main$variable) <- c("H3K4me1","H3K27ac")

pdf("hmc_active.pdf")
p <- ggplot(sample_main, aes(variable, value)) + geom_violin((aes(fill=variable))) + scale_fill_manual(values = c("H3K4me1" = 'palegreen3',"H3K27ac" = 'skyblue3'))
p + theme_classic()
dev.off()


============

library(ggplot2)
dat<- read.table("mc_poised.txt", fill=TRUE, header=T)
df = data.frame(dat)

require(reshape2)
sample_main <- melt(df)
levels(sample_main$variable) <- c("H3K4me1","H3K27me3")

pdf("mc_poised.pdf")
p <- ggplot(sample_main, aes(variable, value)) + geom_violin((aes(fill=variable))) + scale_fill_manual(values = c("H3K4me1" = 'darkgreen',"H3K27me3" = 'red'))
p + theme_classic()
dev.off()




dat<- read.table("hmc_poised.txt", fill=TRUE, header=T)
df = data.frame(dat)

require(reshape2)
sample_main <- melt(df)
levels(sample_main$variable) <- c("H3K4me1","H3K27me3")

pdf("hmc_poised.pdf")
p <- ggplot(sample_main, aes(variable, value)) + geom_violin((aes(fill=variable))) + scale_fill_manual(values = c("H3K4me1" = 'palegreen3',"H3K27me3" = 'indianred1'))
p + theme_classic()
dev.off()


============

library(ggplot2)
dat<- read.table("mc_primed.txt", fill=TRUE, header=T)
df = data.frame(dat)

require(reshape2)
sample_main <- melt(df)
levels(sample_main$variable) <- c("H3K4me1")

pdf("mc_primed.pdf")
p <- ggplot(sample_main, aes(variable, value)) + geom_violin(width=0.3,(aes(fill=variable))) + scale_fill_manual(values = c("H3K4me1" = 'darkgreen'))
p + theme_classic()
dev.off()



dat<- read.table("hmc_primed.txt", fill=TRUE, header=T)
df = data.frame(dat)

require(reshape2)
sample_main <- melt(df)
levels(sample_main$variable) <- c("H3K4me1")

pdf("hmc_primed.pdf")
p <- ggplot(sample_main, aes(variable, value)) + geom_violin(width=0.3, (aes(fill=variable))) + scale_fill_manual(values = c("H3K4me1" = 'palegreen3'))
p + theme_classic()
dev.off()


