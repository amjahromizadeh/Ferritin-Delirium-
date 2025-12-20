library(data.table)
library(coloc)
library(arrow)

del <- fread("Delirium_AF0p005.mr_ready.tsv.gz") 
sqtl_ctx  <- read_parquet("Brain_Cortex.v10.sQTLs.signif_pairs.parquet")
sqtl_ctx  <- as.data.table(sqtl_ctx)

## Case-control fraction for Delirium (NFE): 8,461 cases / 449,979 controls
case_frac <- 8461 / (8461 + 449979)

## APOE region (GRCh38), 500 kb window
chr_apo   <- 19L
start_apo <- 44421094L
end_apo   <- 45421094L

## 1) Prep Delirium GWAS in window
keep_cols <- intersect(
  c("SNP","chr","pos","beta","se","eaf","pval","samplesize"),
  names(del)
)
del <- del[, ..keep_cols]
# types
for (nm in intersect(c("chr","pos","beta","se","eaf","pval","samplesize"), names(del))) {
  del[, (nm) := as.numeric(get(nm))]
}
del_win <- del[chr == chr_apo & pos >= start_apo & pos <= end_apo]
del_win[, maf := pmin(eaf, 1 - eaf)]
del_win <- del_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1]
del_win[, key := paste(chr, pos)]
# collapse to one row per key (keep smallest p if present)
if (!"pval" %in% names(del_win)) del_win[, pval := (beta/se)^2]
setorder(del_win, key, pval)
del_best <- del_win[!duplicated(key)]

## 2) Prep sQTL (Brain Cortex v10) in window
parts <- tstrsplit(sqtl_ctx$variant_id, "_", fixed = TRUE)
sqtl_ctx[, `:=`(
  chr   = as.integer(sub("^chr","", parts[[1]])),
  pos   = as.integer(parts[[2]]),
  ref   = parts[[3]],
  alt   = parts[[4]],
  build = parts[[5]]
)]

# keep APOE window and SNPs only
ctx <- sqtl_ctx[chr == chr_apo & pos >= start_apo & pos <= end_apo]
ctx <- ctx[nchar(ref) == 1 & nchar(alt) == 1]

# ensure slope / slope_se exist
if (!"slope" %in% names(ctx))     setnames(ctx, "beta", "slope",      skip_absent = TRUE)
if (!"slope_se" %in% names(ctx))  setnames(ctx, "se",   "slope_se",   skip_absent = TRUE)
# MAF
if (!"af" %in% names(ctx) && "maf" %in% names(ctx)) ctx[, af := maf]
ctx[, maf := pmin(af, 1 - af)]
ctx <- ctx[!is.na(slope) & !is.na(slope_se) & !is.na(maf) & maf > 0 & maf < 1]
ctx[, key := paste(chr, pos)]

# collapse to one row per key (pick most significant; use any available p-value column)
if (!("pval_nominal" %in% names(ctx))) ctx[, pval_nominal := pmin(pval_perm, pval_beta, na.rm = TRUE)]
setorder(ctx, key, pval_nominal)
ctx_best <- ctx[!duplicated(key)]

## 3) Align by key (inner join → identical rows/order for both traits)
aln <- merge(
  del_best[, .(key, SNP, beta_gwas = beta, se_gwas = se, maf_gwas = maf, N_gwas = samplesize)],
  ctx_best[, .(key, slope, slope_se, maf_sqtl = maf)],
  by = "key"
)

# numeric clean
num_cols <- c("beta_gwas","se_gwas","maf_gwas","N_gwas","slope","slope_se","maf_sqtl")
for (nm in num_cols) aln[, (nm) := as.numeric(get(nm))]
A <- aln[complete.cases(beta_gwas, se_gwas, maf_gwas, N_gwas, slope, slope_se, maf_sqtl)]
cat(sprintf("Aligned, clean rows: %d\n", nrow(A)))
stopifnot(nrow(A) > 0)

## 4) Build coloc inputs from the SAME table
d1 <- list(
  snp     = A$SNP,
  beta    = A$beta_gwas,
  varbeta = A$se_gwas^2,
  MAF     = A$maf_gwas,
  N       = A$N_gwas,
  type    = "cc",
  s       = case_frac
)
d2 <- list(
  snp     = A$SNP,
  beta    = A$slope,
  varbeta = A$slope_se^2,
  MAF     = A$maf_sqtl,
  type    = "quant",
  sdY     = 1
)
stopifnot(length(d1$beta) == length(d2$beta))

## 5) Run coloc
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_ctx <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

cat("\n=== Delirium × GTEx v10 sQTL (Brain Cortex, APOE region) — coloc.abf ===\n")
print(cc_ctx$summary)
cat(sprintf("PP.H4 (shared causal): %.2f%%\n", 100 * cc_ctx$summary["PP.H4.abf"]))

# Top SNPs by per-SNP PP.H4 (if available)
res_ctx <- as.data.table(cc_ctx$results)
if ("SNP.PP.H4" %in% names(res_ctx) && nrow(res_ctx) > 0) {
  setorder(res_ctx, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4:\n")
  print(res_ctx[1:min(.N, 10), .(snp, PP.H4 = SNP.PP.H4)])
}