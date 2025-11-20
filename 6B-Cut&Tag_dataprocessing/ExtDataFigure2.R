
dat<- read.table("H3K27ac_CnT_shared_cpg_mc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27ac_6L_shared_cpg_mc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "6L mc (CpGs shared)",xlab = "H3K27ac enrichment", ylab = "6L (mc Methylation %)", col = "blue3", pch = ".")
mtext("n=42313", side=3)
dev.off()



dat<- read.table("H3K27ac_CnT_shared_cpg_hmc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27ac_6L_shared_cpg_hmc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "6L hmc (CpGs shared)",xlab = "H3K27ac enrichment", ylab = "6L (hmc Methylation %)", col = "dodgerblue3", pch = ".")
mtext("n=42313", side=3)
dev.off()




dat<- read.table("H3K27ac_WG_shared_cpg_mc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27ac_WG_shared_cpg_mc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "WG mc (CpGs shared with 6L)",xlab = "H3K27ac enrichment", ylab = "WG (mc Methylation %)", col = "blue3", pch = ".")
mtext("n=42313", side=3)
dev.off()



dat<- read.table("H3K27ac_WG_shared_cpg_hmc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27ac_WG_shared_cpg_hmc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "WG hmc (CpGs shared with 6L)",xlab = "H3K27ac enrichment", ylab = "WG (hmc Methylation %)", col = "dodgerblue3", pch = ".")
mtext("n=42313", side=3)
dev.off()

#####################################################  H3K4Me1  ########################################




dat<- read.table("H3K4Me1_CnT_shared_cpg_mc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K4Me1_6L_shared_cpg_mc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "6L mc (CpGs shared)",xlab = "H3K4Me1 enrichment", ylab = "6L (mc Methylation %)", col = "green3", pch = ".")
mtext("n=75292", side=3)
dev.off()



dat<- read.table("H3K4Me1_CnT_shared_cpg_hmc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K4Me1_6L_shared_cpg_hmc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "6L hmc (CpGs shared)",xlab = "H3K4Me1 enrichment", ylab = "6L (hmc Methylation %)", col = "palegreen1", pch = ".")
mtext("n=75292", side=3)
dev.off()




dat<- read.table("H3K4Me1_WG_shared_cpg_mc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K4Me1_WG_shared_cpg_mc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "WG mc (CpGs shared with 6L)",xlab = "H3K4Me1 enrichment", ylab = "WG (mc Methylation %)", col = "green3", pch = ".")
mtext("n=75292", side=3)
dev.off()



dat<- read.table("H3K4Me1_WG_shared_cpg_hmc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K4Me1_WG_shared_cpg_hmc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "WG hmc (CpGs shared with 6L)",xlab = "H3K4Me1 enrichment", ylab = "WG (hmc Methylation %)", col = "palegreen1", pch = ".")
mtext("n=75292", side=3)
dev.off()

#################################################  H3K27Me3  ############################################




dat<- read.table("H3K27Me3_CnT_shared_cpg_mc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27Me3_6L_shared_cpg_mc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "6L mc (CpGs shared)",xlab = "H3K27Me3 enrichment", ylab = "6L (mc Methylation %)", col = "indianred1", pch = ".")
mtext("n=27381", side=3)
dev.off()



dat<- read.table("H3K27Me3_CnT_shared_cpg_hmc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27Me3_6L_shared_cpg_hmc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "6L hmc (CpGs shared)",xlab = "H3K27Me3 enrichment", ylab = "6L (hmc Methylation %)", col = "lightsalmon", pch = ".")
mtext("n=27381", side=3)
dev.off()




dat<- read.table("H3K27Me3_WG_shared_cpg_mc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27Me3_WG_shared_cpg_mc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "WG mc (CpGs shared with 6L)",xlab = "H3K27Me3 enrichment", ylab = "WG (mc Methylation %)", col = "indianred1", pch = ".")
mtext("n=27381", side=3)
dev.off()



dat<- read.table("H3K27Me3_WG_shared_cpg_hmc.bed", fill=TRUE, header=T)
df = data.frame(dat)
head(df)
y <- log2(df["Val1"])
x <- (df["Val2"])*100
df2 = data.frame(y,x)
png("H3K27Me3_WG_shared_cpg_hmc.png", res = 300, width = 5, height = 4, units = 'in')
plot(df2, main = "WG hmc (CpGs shared with 6L)",xlab = "H3K27Me3 enrichment", ylab = "WG (hmc Methylation %)", col = "lightsalmon", pch = ".")
mtext("n=27381", side=3)
dev.off()
