library(data.table)
library(coloc)

del      <- fread("Delirium_AF0p005.mr_ready.tsv.gz")
eqtlgen  <- fread("2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

# APOE window (GRCh38): chr19:44,421,094–45,421,094
chr_apo   <- 19L
start_apo <- 44421094L
end_apo   <- 45421094L

gwas <- as.data.table(del)
stopifnot("SNP" %in% names(gwas))

# chr parsing
if (!"chr" %in% names(gwas)) stop("No 'chr' column in Delirium file.")
gwas[, chr := sub("^chr", "", as.character(chr))]
gwas[, chr := as.integer(chr)]

# coerce numerics & compute MAF from EAF
for (nm in intersect(c("beta","se","eaf","pval","samplesize","pos"), names(gwas))) {
  gwas[, (nm) := as.numeric(get(nm))]
}
gwas[, maf := pmin(eaf, 1 - eaf)]

# window + basic QC
gwas_win <- gwas[chr == chr_apo & pos >= start_apo & pos <= end_apo]
gwas_win <- gwas_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1]
# one row per rsID
gwas_rs <- unique(gwas_win, by = "SNP")

## eQTLGen: APOE rows; column detection & standardisation 

eg <- copy(eqtlgen)
pick1 <- function(cols, choices) { x <- intersect(choices, cols); if (length(x)) x[1] else NA_character_ }
cn <- names(eg)

snp_col  <- pick1(cn, c("SNP","snp","rsid","RSID","MarkerName"))
gene_col <- pick1(cn, c("geneSymbol","GeneSymbol","HUGO","Symbol","Gene","gene"))
beta_col <- pick1(cn, c("beta","Beta","BETA","b","B","Effect","effect","NES"))
se_col   <- pick1(cn, c("se","SE","se_beta","sebeta","StdErr","stderr","SE.Effect"))
p_col    <- pick1(cn, c("p","P","pval","pvalue","Pvalue","p_val","P-value","PValue"))
N_col    <- pick1(cn, c("N","samplesize","SampleSize","NrSamples","n_samples","n"))

stopifnot(!is.na(snp_col), !is.na(gene_col))
setnames(eg, snp_col,  "SNP")
setnames(eg, gene_col, "geneSymbol")

have_beta_se <- !is.na(beta_col) && !is.na(se_col)
if (have_beta_se) {
  setnames(eg, beta_col, "beta")
  setnames(eg, se_col,   "se")
} else {
  stopifnot(!is.na(p_col))         # need p if no beta/se
  setnames(eg, p_col, "pval")
  if (!is.na(N_col)) setnames(eg, N_col, "N")
}

eg_apoe <- eg[geneSymbol == "APOE"]
stopifnot(nrow(eg_apoe) > 0)

## Align strictly by rsID
aln <- merge(
  gwas_rs[, .(SNP, beta_gwas = beta, se_gwas = se, maf_gwas = maf)],
  if (have_beta_se) {
    eg_apoe[, .(SNP, beta_eqtl = beta, se_eqtl = se)]
  } else {
    eg_apoe[, .(SNP, pval_eqtl = pval, N_eqtl = if ("N" %in% names(eg_apoe)) as.numeric(N) else NA_real_)]
  },
  by = "SNP"
)

# Drop incomplete GWAS rows
aln <- aln[!is.na(beta_gwas) & !is.na(se_gwas)]

# If p-value path, require N; fill remaining NA N with median of observed N
if (!have_beta_se) {
  aln <- aln[!is.na(pval_eqtl)]
  if (!("N_eqtl" %in% names(aln))) stop("eQTLGen sample size (N) column not found in your file.")
  if (any(is.na(aln$N_eqtl))) {
    N_med <- suppressWarnings(as.numeric(median(aln$N_eqtl, na.rm = TRUE)))
    if (is.finite(N_med)) aln[is.na(N_eqtl), N_eqtl := N_med]
  }
  aln <- aln[!is.na(N_eqtl) & is.finite(N_eqtl) & N_eqtl > 0]
}

stopifnot(nrow(aln) > 0)
cat(sprintf("Aligned rsID overlap: %d SNPs\n", nrow(aln)))

## Build coloc datasets
N_cases <- 8461
N_ctrls <- 449979

d1 <- list(
  snp     = aln$SNP,
  beta    = aln$beta_gwas,
  varbeta = aln$se_gwas^2,
  MAF     = gwas_rs[J(aln$SNP), on="SNP"]$maf,  # ensure same length; equivalent to aln$maf_gwas
  type    = "cc",
  s       = N_cases / (N_cases + N_ctrls),
  N       = N_cases + N_ctrls
)

if (have_beta_se) {
  d2 <- list(
    snp     = aln$SNP,
    beta    = aln$beta_eqtl,
    varbeta = aln$se_eqtl^2,
    type    = "quant",
    sdY     = 1
  )
} else {
  d2 <- list(
    snp     = aln$SNP,
    pvalues = aln$pval_eqtl,
    N       = as.numeric(aln$N_eqtl),
    type    = "quant",
    sdY     = 1
  )
  # supply MAF to help sdY estimation; reuse GWAS MAFs if eQTL MAF absent
  d2$MAF <- d1$MAF
}

## Run coloc
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

cat("\n=== Delirium × APOE — eQTLGen Whole Blood (rsID merge) — coloc.abf ===\n")
print(cc$summary)
cat(sprintf("PP.H4 (shared causal): %.3f%%\n", 100 * cc$summary["PP.H4.abf"]))

res <- as.data.table(cc$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res)) {
  setorder(res, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4:\n")
  print(res[1:min(.N,10), .(SNP = snp, PP.H4 = SNP.PP.H4)])
} else {
  cat("\n(No per-SNP PP.H4 returned — likely p-value-only mode or tiny overlap.)\n")
}