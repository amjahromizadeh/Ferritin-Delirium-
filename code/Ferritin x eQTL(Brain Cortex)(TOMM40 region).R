library(data.table)
library(coloc) 

Ferr <- fread("Ferritin_AF0p005.mr_ready.tsv.gz")
cor_sig <- fread("Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz")
cor_genes <- fread("Brain_Cortex.v8.egenes.txt.gz")

## 1) Define TOMM40 window (GRCh38)
chr_tomm  <- 19L
start_bp  <- 44401715L   # 44,401,715
end_bp    <- 45401715L   # 45,401,715

## 2) Clean Ferritin GWAS in window
Ferr_dt <- as.data.table(Ferr)
# make sure numeric
for (nm in intersect(c("chr","pos","beta","se","eaf","samplesize"), names(Ferr_dt))) {
  Ferr_dt[, (nm) := as.numeric(get(nm))]
}
Ferr_win <- Ferr_dt[chr == chr_tomm & pos >= start_bp & pos <= end_bp]
Ferr_win[, maf := pmin(eaf, 1 - eaf)]
Ferr_win <- Ferr_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf >= 0.005]
Ferr_win <- unique(Ferr_win, by = "SNP")

## 3) Pull TOMM40 eQTL rows from GTEx v8 Brain Cortex
gid <- cor_genes[gene_name == "TOMM40", unique(gene_id)]
if (length(gid) == 0L) stop("TOMM40 not found in Brain Cortex eGenes.")
eqtl <- cor_sig[gene_id %in% gid]

# parse variant_id like 'chr19_44903689_A_G_b38'
parts <- tstrsplit(eqtl$variant_id, "_", fixed = TRUE)
eqtl[, `:=`(
  chr = as.integer(sub("^chr","", parts[[1]])),
  pos = as.integer(parts[[2]]),
  ref = parts[[3]],
  alt = parts[[4]]
)]
# keep SNPs only; limit to TOMM40 window
eqtl <- eqtl[nchar(ref) == 1 & nchar(alt) == 1]
eqtl <- eqtl[chr == chr_tomm & pos >= start_bp & pos <= end_bp]

# ensure slopes/se are present (GTEx v8 has slope/slope_se)
if (!"slope" %in% names(eqtl) && "beta" %in% names(eqtl)) setnames(eqtl, "beta", "slope")
if (!"slope_se" %in% names(eqtl) && "se" %in% names(eqtl)) setnames(eqtl, "se", "slope_se")

# get allele frequency → maf (prefer 'maf', else compute from ma_count/ma_samples, else 'af')
if (!"maf" %in% names(eqtl)) {
  if (all(c("ma_count","ma_samples") %in% names(eqtl))) {
    eqtl[, af := as.numeric(ma_count) / as.numeric(ma_samples)]
  } else if ("af" %in% names(eqtl)) {
    eqtl[, af := as.numeric(af)]
  }
  if ("af" %in% names(eqtl)) eqtl[, maf := pmin(af, 1 - af)]
}
eqtl <- eqtl[!is.na(slope) & !is.na(slope_se) & !is.na(maf) & maf >= 0.005]
eqtl <- unique(eqtl, by = "variant_id")

## 4) Align by position 
Ferr_win[, key := paste(chr, pos)]
eqtl[,     key := paste(chr, pos)]
keys <- intersect(Ferr_win$key, eqtl$key)
Ferr_aln <- Ferr_win[key %in% keys]
eQTL_aln <- eqtl[key %in% keys]
setorder(Ferr_aln, key); setorder(eQTL_aln, key)

if (nrow(Ferr_aln) == 0L || nrow(eQTL_aln) == 0L) {
  stop("No overlapping SNPs in the TOMM40 window after QC.")
}

## 5) Build coloc datasets
d1 <- list(
  snp     = Ferr_aln$SNP,
  beta    = Ferr_aln$beta,
  varbeta = Ferr_aln$se^2,
  MAF     = Ferr_aln$maf,
  N       = Ferr_aln$samplesize,
  type    = "quant"
)

# QTL effect sizes are on (approximately) standardized scale → set sdY=1
d2 <- list(
  snp     = Ferr_aln$SNP,             # match order/IDs with d1
  beta    = eQTL_aln$slope,
  varbeta = eQTL_aln$slope_se^2,
  MAF     = eQTL_aln$maf,
  type    = "quant",
  sdY     = 1
)

## 6) Run COLOC
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_tomm40_ctx <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

print(cc_tomm40_ctx$summary)

## 7) top per-SNP H4
res <- as.data.table(cc_tomm40_ctx$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res) > 0) {
  setorder(res, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4 (Brain Cortex):\n")
  print(res[1:min(.N, 10), .(snp, PP.H4 = SNP.PP.H4)])
} else {
  cat("\nNo per-SNP PP.H4 returned (very small overlap or model degenerate).\n")
}