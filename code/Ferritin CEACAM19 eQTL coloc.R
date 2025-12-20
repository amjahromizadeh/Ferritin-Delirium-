library(data.table)
library(coloc) 

Ferr <- fread("Ferritin_AF0p005.mr_ready.tsv.gz")
cor_sig <- fread("Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz")
cor_genes <- fread("Brain_Cortex.v8.egenes.txt.gz")

## ---- CEACAM19 window (GRCh38) ----
chr_ceac   <- 19L
start_ceac <- 44119993L
end_ceac   <- 45119993L

## ---- 1) Ferritin GWAS in window ----
Ferr_dt <- as.data.table(Ferr)
for (nm in intersect(c("beta","se","eaf","samplesize","chr","pos"), names(Ferr_dt))) {
  Ferr_dt[, (nm) := as.numeric(get(nm))]
}
gwas_win <- Ferr_dt[chr == chr_ceac & pos >= start_ceac & pos <= end_ceac]
gwas_win[, maf := pmin(eaf, 1 - eaf)]
gwas_win <- gwas_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf >= 0.005]
gwas_win <- unique(gwas_win, by = "SNP")

## ---- 2) CEACAM19 eQTL (GTEx v8 Brain Cortex) in window ----
# Map gene symbol -> Ensembl IDs using eGenes
gid <- cor_genes[gene_name == "CEACAM19", unique(gene_id)]
if (length(gid) == 0L) stop("CEACAM19 not found in Brain Cortex v8 eGenes.")

eqtl <- as.data.table(cor_sig)[gene_id %in% gid]
if (nrow(eqtl) == 0L) stop("CEACAM19 has 0 significant eQTL pairs in Brain Cortex v8 (signif pairs file).")

# Parse variant_id: chr_pos_ref_alt_b38
parts <- tstrsplit(eqtl$variant_id, "_", fixed = TRUE)
if (length(parts) < 5) stop("Unexpected variant_id format; expected chr_pos_ref_alt_b38.")
eqtl[, `:=`(
  chr   = as.integer(sub("^chr","", parts[[1]])),
  pos   = as.integer(parts[[2]]),
  ref   = parts[[3]],
  alt   = parts[[4]],
  build = parts[[5]]
)]

# Window & SNP-only
eqtl <- eqtl[chr == chr_ceac & pos >= start_ceac & pos <= end_ceac]
eqtl <- eqtl[nchar(ref) == 1 & nchar(alt) == 1]

# Make sure effect columns are present (v8 uses slope/slope_se)
if (!"slope" %in% names(eqtl) && "beta" %in% names(eqtl))    setnames(eqtl, "beta", "slope",     skip_absent = TRUE)
if (!"slope_se" %in% names(eqtl) && "se" %in% names(eqtl))   setnames(eqtl, "se",   "slope_se",  skip_absent = TRUE)

# Minor allele frequency (prefer 'af'; fallback to ma_count/ma_samples if needed)
if (!"af" %in% names(eqtl) && all(c("ma_count","ma_samples") %in% names(eqtl))) {
  eqtl[, af := as.numeric(ma_count) / (2 * as.numeric(ma_samples))]
}
if (!"af" %in% names(eqtl)) stop("No allele frequency in eQTL table (need 'af' or ma_count/ma_samples).")
eqtl[, maf := pmin(af, 1 - af)]
eqtl <- eqtl[!is.na(slope) & !is.na(slope_se) & !is.na(maf) & maf >= 0.005]
eqtl <- unique(eqtl, by = "variant_id")

## ---- 3) Position-align (chr:pos key) ----
gwas_win[, key := paste(chr, pos)]
eqtl[,     key := paste(chr, pos)]

keys <- intersect(gwas_win$key, eqtl$key)
gw  <- gwas_win[key %in% keys]
eq  <- eqtl[key %in% keys]

data.table::setorder(gw, key)
data.table::setorder(eq, key)

if (nrow(gw) == 0L || nrow(eq) == 0L) stop("No overlapping SNPs after window/MAF/SNP-only filters.")

## ---- 4) Build coloc datasets ----
d1 <- list(
  snp     = gw$SNP,
  beta    = gw$beta,
  varbeta = gw$se^2,
  MAF     = gw$maf,
  N       = gw$samplesize,
  type    = "quant"
)

d2 <- list(
  snp     = gw$SNP,              # keep same order/IDs
  beta    = eq$slope,
  varbeta = eq$slope_se^2,
  MAF     = eq$maf,
  type    = "quant",
  sdY     = 1                    # reasonable when expression units are scaled
)

## ---- 5) COLOC ----
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_ceacam19_ctx <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

print(cc_ceacam19_ctx$summary)

# Per-SNP H4 weights (if available)
res <- as.data.table(cc_ceacam19_ctx$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res) > 0) {
  data.table::setorder(res, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4 (Brain Cortex):\n")
  print(res[1:min(.N, 10), .(snp, PP.H4 = SNP.PP.H4)])
} else {
  cat("\nNo per-SNP PP.H4 returned (small overlap or degenerate model).\n")
}