# dss_like_from_wide_with_scatter_white_bed_raster_jitter.R
# Usage:
#   Rscript Supp.Figure_8.R input.tsv output.tsv scatter.pdf peaks.bed
#
# Description:
#   - Reads a wide-format table with 5 replicate blocks (each block: Chr Start End Meth_Prop #mC #C)
#     * First 3 blocks = treatment 1  → mu2
#     * Last 2 blocks  = treatment 0  → mu1
#   - Computes DSS-like mu1/mu2, p/q, meth.diff (= 100*(mu2 - mu1))
#   - Marks loci as "significant" (red) IFF they overlap any interval in the supplied BED file
#     (no p-value or effect-size threshold used for significance)
#   - Saves a jittered, rasterized scatter plot (mu1 vs mu2) with white background to PDF
#   - Writes full results to output.tsv and significant rows to tmp.txt

suppressWarnings(suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrastr)   # for geom_point_rast()
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) 
{
  stop("Usage: Rscript Supp.Figure_SW.R input.tsv output.tsv scatter.pdf peaks.bed")
}

infile   <- args[1]
outfile  <- args[2]
plotfile <- args[3]
bedfile  <- args[4]

dt <- fread(infile)

# ---------- Helper functions (DSS-like pieces) ----------
rowVars <- function (x, center = NULL, ...) {
  n <- !is.na(x); n <- rowSums(n); n[n <= 1] <- NA
  if (is.null(center)) center <- rowMeans(x, ...)
  x <- x - center; x <- x * x; x <- rowSums(x, ...); x/(n - 1)
}
compute_mean_vec <- function(X, N) (rowSums(X) + 0.5) / (rowSums(N) + 1)
dbb <- function (size, x, mu, phi, log=TRUE)  {
  tmp <- 1/phi - 1
  alpha <- mu * tmp; beta <- tmp - alpha
  v <- lchoose(size, x) - lbeta(beta, alpha) + lbeta(size - x + beta, x + alpha)
  if (!log) exp(v) else v
}
dispersion_shrinkage <- function(X, N, prior, estprob) {
  plik <- function(size, X, mu, m0, tau, phi_log) {
    -(sum(dbb(size, X, mu, exp(phi_log))) + dnorm(phi_log, mean=m0, sd=tau, log=TRUE))
  }
  shrk_phi <- rep(exp(prior[1]), nrow(N))
  ix <- rowSums(N > 0) > 0
  if (any(ix)) {
    X2 <- X[ix,,drop=FALSE]; N2 <- N[ix,,drop=FALSE]; mu2 <- estprob[ix]
    out <- numeric(nrow(X2))
    for (i in seq_len(nrow(X2))) {
      opt <- optimize(plik, interval=c(-5, log(0.99)),
                      size=N2[i,], X=X2[i,], mu=mu2[i],
                      m0=prior[1], tau=prior[2], tol=1e-4)
      out[i] <- exp(opt$minimum)
    }
    shrk_phi[ix] <- out
  }
  shrk_phi
}
est_prior_logN <- function(X, N) {
  ix <- rowMeans(N > 10) == 1 & rowSums(N == 0) == 0
  X <- X[ix,,drop=FALSE]; N <- N[ix,,drop=FALSE]
  if (!nrow(X)) return(c(log(0.05), 1))
  p  <- X / N
  mm <- rowMeans(p); mm[mm==0] <- 1e-5; mm[mm==1] <- 1-1e-5
  vv <- rowVars(p)
  phi <- vv/mm/(1-mm)
  phi <- phi[vv > 0 & is.finite(phi) & phi > 0]
  if (!length(phi)) return(c(log(0.05), 1))
  lphi <- log(phi)
  c(median(lphi, na.rm=TRUE), IQR(lphi, na.rm=TRUE)/1.39)
}
wald_stats <- function(mu1, mu2, n1, n2, phi1, phi2) {
  dif <- mu1 - mu2
  n1m <- rowSums(n1); n2m <- rowSums(n2)
  var1 <- rowSums(n1 * mu1 * (1 - mu1) * (1 + (n1 - 1) * phi1)) / (n1m^2)
  var2 <- rowSums(n2 * mu2 * (1 - mu2) * (1 + (n2 - 1) * phi2)) / (n2m^2)
  vv <- var1 + var2; vv[vv < 1e-5] <- 1e-5
  se <- sqrt(vv); stat <- dif / se
  pval <- 2 * (1 - pnorm(abs(stat)))
  data.table(diff = dif, diff.se = se, stat = stat, pval = pval)
}

# ---------- Parse wide-format data (5 blocks x 6 cols) ----------
block_size <- 6L
if (ncol(dt) %% block_size != 0L) stop("Column count not a multiple of 6.")
n_blocks <- ncol(dt) %/% block_size
if (n_blocks != 5L) stop(sprintf("Expected 5 replicate blocks, found %d.", n_blocks))

blk_idx <- function(b) {
  s <- (b - 1L) * block_size
  list(chr = 1L + s, start = 2L + s, end = 3L + s, mC = 5L + s, C = 6L + s)
}
i1 <- blk_idx(1L)
Chr   <- dt[[i1$chr]]
Start <- dt[[i1$start]]
End   <- dt[[i1$end]]

get_XN <- function(b) {
  i <- blk_idx(b)
  X <- as.numeric(dt[[i$mC]])
  N <- as.numeric(dt[[i$mC]] + dt[[i$C]])
  list(X = X, N = N)
}

# Treatment 1 (blocks 1–3) -> mu2
b1 <- get_XN(1L); b2 <- get_XN(2L); b3 <- get_XN(3L)
x_trt1 <- cbind(b1$X, b2$X, b3$X)
n_trt1 <- cbind(b1$N, b2$N, b3$N)

# Treatment 0 (blocks 4–5) -> mu1
b4 <- get_XN(4L); b5 <- get_XN(5L)
x_trt0 <- cbind(b4$X, b5$X)
n_trt0 <- cbind(b4$N, b5$N)

# ---------- DSS-like computation ----------
mu1 <- compute_mean_vec(x_trt0, n_trt0)   # treatment 0
mu2 <- compute_mean_vec(x_trt1, n_trt1)   # treatment 1
prior0 <- est_prior_logN(x_trt0, n_trt0)
prior1 <- est_prior_logN(x_trt1, n_trt1)
phi0 <- dispersion_shrinkage(x_trt0, n_trt0, prior0, mu1)
phi1 <- dispersion_shrinkage(x_trt1, n_trt1, prior1, mu2)
phi0m <- matrix(phi0, nrow=length(phi0), ncol=ncol(n_trt0))
phi1m <- matrix(phi1, nrow=length(phi1), ncol=ncol(n_trt1))
wd <- wald_stats(mu1, mu2, n_trt0, n_trt1, phi0m, phi1m)
qval <- p.adjust(wd$pval, method = "fdr")

meth_diff <- -100 * (mu1 - mu2)  # == 100*(mu2 - mu1)
out <- data.table(
  Chr   = Chr,
  Start = Start,
  End   = End,
  mu1   = mu1,
  mu2   = mu2,
  pvalue = wd$pval,
  qvalue = qval,
  meth.diff = meth_diff
)

# ---------- Read BED file and flag overlaps as significant ----------
bed <- fread(bedfile, header = FALSE, fill = TRUE, data.table = TRUE)
if (ncol(bed) < 3) stop("BED file must have at least 3 columns: chrom, start, end.")
bed <- bed[, .(V1, V2, V3)]
setnames(bed, c("Chr","start","end"))
bed <- bed[!grepl("^(track|browser|#)", Chr, ignore.case = TRUE)]
bed[, `:=`(
  start = suppressWarnings(as.integer(start)),
  end   = suppressWarnings(as.integer(end))
)]
bed <- bed[!is.na(start) & !is.na(end) & end >= start]

# If BED is 0-based half-open and sites are 1-based closed, you may need:
# bed[, start := start + 1L]

# Overlaps (data.table::foverlaps requires keys on Chr,start,end)
out[, idx := .I]
gr_out <- out[, .(idx, Chr, start = Start, end = End)]
setkey(gr_out, Chr, start, end)
setkey(bed,    Chr, start, end)

ov <- foverlaps(bed, gr_out, nomatch = 0L)
sig_idx <- unique(ov$idx)

out[, sig := FALSE]
if (length(sig_idx)) out[sig_idx, sig := TRUE]

# ---------- Write outputs ----------
fwrite(out[, -"idx"], outfile, sep = "\t", quote = FALSE)
fwrite(out[sig == TRUE, -"idx"], "tmp.txt", sep = "\t", quote = FALSE)
cat(sprintf("Wrote %d total rows to %s\n", nrow(out), outfile))
cat(sprintf("Wrote %d BED-overlap (red) rows to tmp.txt\n", nrow(out[sig == TRUE])))

# ---------- Plot (rasterized, smaller points with jitter) ----------
p <- ggplot(out, aes(x = mu1, y = mu2)) +
  geom_point_rast(
    data = out[sig == FALSE],
    color = "grey70",
    alpha = 0.3,
    size = 0.3,
    position = position_jitter(width = 0.005, height = 0.005)
  ) +
  geom_point_rast(
    data = out[sig == TRUE],
    color = "#DF7D2E",#"#F4BC43", # "#296219", #"#65C84B",
    alpha = 0.8,
    size = 0.3,
    position = position_jitter(width = 0.005, height = 0.005)
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = expression("Scatter of H3K4me3 5mC  Vs DUET 5mC"),
    x = "DUET 5mC",
    y = "H3K4Me3 5mC"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.3)
  )

# ---------- Save PDF ----------
# If user passed a non-PDF name, convert to .pdf automatically
if (!grepl("\\.pdf$", plotfile, ignore.case = TRUE)) {
  plotfile <- sub("\\.[^.]+$", ".pdf", plotfile)
}

ggsave(
  filename = plotfile,
  plot = p,
  width = 6,
  height = 6,
  device = cairo_pdf   # high-quality vector PDF
)

cat(sprintf("Saved jittered rasterized scatter plot as PDF to %s\n", plotfile))

