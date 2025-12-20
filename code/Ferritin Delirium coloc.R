library(data.table)
library(coloc)
library(knitr)

# Load GWAS Datasets
ferr <- fread("Ferritin_AF0p005.mr_ready.tsv.gz")
del  <- fread("Delirium_AF0p005.mr_ready.tsv.gz")

# 1) Delirium case fraction
N_cases <- 8461
N_ctrls <- 449979
case_frac <- N_cases / (N_cases + N_ctrls)

# 2) Clean basic types
for (dtname in c("ferr","del")) {
  dt <- get(dtname)
  if (!"chr" %in% names(dt)) stop(sprintf("%s is missing 'chr'", dtname))
  dt[, chr := as.integer(sub("^chr","", as.character(chr)))]
  for (nm in intersect(c("beta","se","eaf","pval","samplesize","pos"), names(dt))) {
    dt[, (nm) := as.numeric(get(nm))]
  }
  assign(dtname, dt)
}

# 3) Priors
p1 <- 1e-4; p2 <- 1e-4; p12 <- 1e-5

# Region A: APOE/TOMM40/CEACAM19 (chr19:44,415,533–45,415,533) 
chr <- 19L; start <- 44415533L; end <- 45415533L
cat(sprintf("\n== APOE/TOMM40/CEACAM19: chr%d:%d-%d ==\n", chr, start, end))

# Ferritin slice
f_win <- ferr[chr == chr & pos >= start & pos <= end]
f_win[, maf := pmin(eaf, 1 - eaf)]
f_win <- f_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]
f_pos <- if ("pval" %in% names(f_win)) f_win[order(pval)][, .SD[1], by = .(chr, pos)] else f_win[order(se)][, .SD[1], by = .(chr, pos)]
f_pos[, key := paste(chr, pos, sep=":")]

# Delirium slice
d_win <- del[chr == chr & pos >= start & pos <= end]
d_win[, maf := pmin(eaf, 1 - eaf)]
d_win <- d_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]
d_pos <- if ("pval" %in% names(d_win)) d_win[order(pval)][, .SD[1], by = .(chr, pos)] else d_win[order(se)][, .SD[1], by = .(chr, pos)]
d_pos[, key := paste(chr, pos, sep=":")]

aln <- merge(f_pos, d_pos, by = c("chr","pos","key"), suffixes = c(".ferr",".del"))
cat(sprintf("Ferritin positions: %d | Delirium positions: %d | Overlap: %d\n", nrow(f_pos), nrow(d_pos), nrow(aln)))

if (nrow(aln) > 0) {
  d1 <- list(snp=aln$key, beta=aln$beta.ferr, varbeta=aln$se.ferr^2, MAF=aln$maf.ferr, type="quant", N=aln$samplesize.ferr, sdY=1)
  d2 <- list(snp=aln$key, beta=aln$beta.del,  varbeta=aln$se.del^2,  MAF=aln$maf.del,  type="cc",    s=case_frac,          N=N_cases+N_ctrls)
  cc <- coloc.abf(d1, d2, p1=p1, p2=p2, p12=p12)
  print(cc$summary)
  cat(sprintf("PP.H4 (shared causal): %.3e (%.3e%%)\n", cc$summary["PP.H4.abf"], 100*cc$summary["PP.H4.abf"]))
  res <- as.data.table(cc$results)
  if ("SNP.PP.H4" %in% names(res) && nrow(res)>0) {
    setorder(res, -SNP.PP.H4)
    cat("\nTop 10 SNPs by SNP.PP.H4:\n"); print(res[1:min(10,.N), .(snp, PP.H4=signif(SNP.PP.H4,3))])
    res[, cum := cumsum(SNP.PP.H4)]; cs_idx <- which(res$cum >= 0.95)[1]; cs_size <- ifelse(length(cs_idx), cs_idx, nrow(res))
    cat(sprintf("\n95%% credible set size: %d | Lead SNP: %s\n", cs_size, res[1, snp]))
  } else {
    cat("\n(No per-SNP PP.H4 returned.)\n")
  }
} else {
  cat("No overlapping SNPs after QC.\n")
}

# Region B: SLC11A2 (chr12:50,760,339–51,260,339)
chr <- 12L; start <- 50760339L; end <- 51260339L
cat(sprintf("\n== SLC11A2: chr%d:%d-%d ==\n", chr, start, end))

f_win <- ferr[chr == chr & pos >= start & pos <= end]
f_win[, maf := pmin(eaf, 1 - eaf)]
f_win <- f_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]
f_pos <- if ("pval" %in% names(f_win)) f_win[order(pval)][, .SD[1], by = .(chr, pos)] else f_win[order(se)][, .SD[1], by = .(chr, pos)]
f_pos[, key := paste(chr, pos, sep=":")]

d_win <- del[chr == chr & pos >= start & pos <= end]
d_win[, maf := pmin(eaf, 1 - eaf)]
d_win <- d_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]
d_pos <- if ("pval" %in% names(d_win)) d_win[order(pval)][, .SD[1], by = .(chr, pos)] else d_win[order(se)][, .SD[1], by = .(chr, pos)]
d_pos[, key := paste(chr, pos, sep=":")]

aln <- merge(f_pos, d_pos, by = c("chr","pos","key"), suffixes = c(".ferr",".del"))
cat(sprintf("Ferritin positions: %d | Delirium positions: %d | Overlap: %d\n", nrow(f_pos), nrow(d_pos), nrow(aln)))

if (nrow(aln) > 0) {
  d1 <- list(snp=aln$key, beta=aln$beta.ferr, varbeta=aln$se.ferr^2, MAF=aln$maf.ferr, type="quant", N=aln$samplesize.ferr, sdY=1)
  d2 <- list(snp=aln$key, beta=aln$beta.del,  varbeta=aln$se.del^2,  MAF=aln$maf.del,  type="cc",    s=case_frac,          N=N_cases+N_ctrls)
  cc <- coloc.abf(d1, d2, p1=p1, p2=p2, p12=p12)
  print(cc$summary)
  cat(sprintf("PP.H4 (shared causal): %.3e (%.3e%%)\n", cc$summary["PP.H4.abf"], 100*cc$summary["PP.H4.abf"]))
  res <- as.data.table(cc$results)
  if ("SNP.PP.H4" %in% names(res) && nrow(res)>0) {
    setorder(res, -SNP.PP.H4)
    cat("\nTop 10 SNPs by SNP.PP.H4:\n"); print(res[1:min(10,.N), .(snp, PP.H4=signif(SNP.PP.H4,3))])
    res[, cum := cumsum(SNP.PP.H4)]; cs_idx <- which(res$cum >= 0.95)[1]; cs_size <- ifelse(length(cs_idx), cs_idx, nrow(res))
    cat(sprintf("\n95%% credible set size: %d | Lead SNP: %s\n", cs_size, res[1, snp]))
  } else {
    cat("\n(No per-SNP PP.H4 returned.)\n")
  }
} else {
  cat("No overlapping SNPs after QC.\n")
}

# Region C: TF (chr3:133,515,185–134,015,185) 
chr <- 3L; start <- 133515185L; end <- 134015185L
cat(sprintf("\n== TF: chr%d:%d-%d ==\n", chr, start, end))

f_win <- ferr[chr == chr & pos >= start & pos <= end]
f_win[, maf := pmin(eaf, 1 - eaf)]
f_win <- f_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]
f_pos <- if ("pval" %in% names(f_win)) f_win[order(pval)][, .SD[1], by = .(chr, pos)] else f_win[order(se)][, .SD[1], by = .(chr, pos)]
f_pos[, key := paste(chr, pos, sep=":")]

d_win <- del[chr == chr & pos >= start & pos <= end]
d_win[, maf := pmin(eaf, 1 - eaf)]
d_win <- d_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]
d_pos <- if ("pval" %in% names(d_win)) d_win[order(pval)][, .SD[1], by = .(chr, pos)] else d_win[order(se)][, .SD[1], by = .(chr, pos)]
d_pos[, key := paste(chr, pos, sep=":")]

aln <- merge(f_pos, d_pos, by = c("chr","pos","key"), suffixes = c(".ferr",".del"))
cat(sprintf("Ferritin positions: %d | Delirium positions: %d | Overlap: %d\n", nrow(f_pos), nrow(d_pos), nrow(aln)))

if (nrow(aln) > 0) {
  d1 <- list(snp=aln$key, beta=aln$beta.ferr, varbeta=aln$se.ferr^2, MAF=aln$maf.ferr, type="quant", N=aln$samplesize.ferr, sdY=1)
  d2 <- list(snp=aln$key, beta=aln$beta.del,  varbeta=aln$se.del^2,  MAF=aln$maf.del,  type="cc",    s=case_frac,          N=N_cases+N_ctrls)
  cc <- coloc.abf(d1, d2, p1=p1, p2=p2, p12=p12)
  print(cc$summary)
  cat(sprintf("PP.H4 (shared causal): %.3e (%.3e%%)\n", cc$summary["PP.H4.abf"], 100*cc$summary["PP.H4.abf"]))
  res <- as.data.table(cc$results)
  if ("SNP.PP.H4" %in% names(res) && nrow(res)>0) {
    setorder(res, -SNP.PP.H4)
    cat("\nTop 10 SNPs by SNP.PP.H4:\n"); print(res[1:min(10,.N), .(snp, PP.H4=signif(SNP.PP.H4,3))])
    res[, cum := cumsum(SNP.PP.H4)]; cs_idx <- which(res$cum >= 0.95)[1]; cs_size <- ifelse(length(cs_idx), cs_idx, nrow(res))
    cat(sprintf("\n95%% credible set size: %d | Lead SNP: %s\n", cs_size, res[1, snp]))
  } else {
    cat("\n(No per-SNP PP.H4 returned.)\n")
  }
} else {
  cat("No overlapping SNPs after QC.\n")
}
