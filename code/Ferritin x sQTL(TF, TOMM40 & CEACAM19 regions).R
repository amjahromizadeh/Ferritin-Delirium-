library(data.table)
library(coloc) 
library(arrow)


Ferr <- fread("Ferritin_AF0p005.mr_ready.tsv.gz")

# GTEx v10 sQTL — Liver (for TF, ±250 kb)
liver_sig <- read_parquet("Liver.v10.sQTLs.signif_pairs.parquet" )
liver_genes <- fread("Liver.v10.sGenes.txt.gz") 

# GTEx v10 sQTL — Brain (Cortex) (for TOMM40 and CEACAM19, ±500 kb) 
braincor_sig <- read_parquet("Brain_Cortex.v10.sQTLs.signif_pairs.parquet")
braincor_genes <- fread("Brain_Cortex.v10.sGenes.txt.gz")

# GTEx v10 sQTL — Brain (Frontal Cortex BA9) (also for TOMM40 and CEACAM19, ±500 kb) 
brainfro_sig <- read_parquet("Brain_Frontal_Cortex_BA9.v10.sQTLs.signif_pairs.parquet")
brainfro_genes <- fread("Brain_Frontal_Cortex_BA9.v10.sGenes.txt.gz")
#-------------------------------------------------------------------------------
## ---- 1) Define TF window (GRCh38) ----
chr_tf   <- 3L
start_tf <- 133515185L
end_tf   <- 134015185L

## ---- 2) Prep Ferritin GWAS in window ----
Ferr_dt <- as.data.table(Ferr)

# Make sure numeric types and compute MAF
numcols <- intersect(c("beta","se","eaf","samplesize","chr","pos"), names(Ferr_dt))
for (nm in numcols) Ferr_dt[, (nm) := as.numeric(get(nm))]

Ferr_win <- Ferr_dt[chr == chr_tf & pos >= start_tf & pos <= end_tf]
Ferr_win[, maf := pmin(eaf, 1 - eaf)]
Ferr_win <- Ferr_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf >= 0.005]
Ferr_win <- unique(Ferr_win, by = "SNP")

# Build coloc dataset for Ferritin (quantitative)
d_ferr <- list(
  snp     = Ferr_win$SNP,
  beta    = Ferr_win$beta,
  varbeta = Ferr_win$se^2,
  MAF     = Ferr_win$maf,
  N       = Ferr_win$samplesize,   # per-SNP OK; scalar also fine
  type    = "quant"
)

## ---- 3) Prep TF sQTL (Liver) in window ----
tf_ids <- liver_genes[gene_name == "TF", unique(gene_id)]
if (length(tf_ids) == 0) stop("TF not found in liver sGenes.txt.gz")

# Keep rows for TF only
liver_tf <- as.data.table(liver_sig)[gene_id %in% tf_ids]

# Parse variant_id -> chr/pos/ref/alt/build
parts <- tstrsplit(liver_tf$variant_id, "_", fixed = TRUE)
liver_tf[, `:=`(
  chr   = as.integer(sub("^chr", "", parts[[1]])),
  pos   = as.integer(parts[[2]]),
  ref   = parts[[3]],
  alt   = parts[[4]],
  build = parts[[5]]
)]

# Restrict to TF window and biallelic SNPs (length 1 ref/alt)
liver_tf <- liver_tf[chr == chr_tf & pos >= start_tf & pos <= end_tf]
liver_tf <- liver_tf[nchar(ref) == 1 & nchar(alt) == 1]

# Harmonise column names (sQTL has slope/slope_se; sometimes beta/se)
if (!"slope" %in% names(liver_tf))  setnames(liver_tf, "beta", "slope",      skip_absent = TRUE)
if (!"slope_se" %in% names(liver_tf)) setnames(liver_tf, "se",   "slope_se",   skip_absent = TRUE)

# Minor allele frequency for sQTL side 
if (!"af" %in% names(liver_tf) && "maf" %in% names(liver_tf)) liver_tf[, af := maf]
liver_tf[, maf := pmin(af, 1 - af)]
liver_tf <- liver_tf[!is.na(slope) & !is.na(slope_se) & !is.na(maf) & maf >= 0.005]

# Deduplicate by variant_id
liver_tf <- unique(liver_tf, by = "variant_id")

# Build coloc dataset for TF sQTL (quantitative). 
# GTEx sQTL slopes are on ~standardised scale; use sdY = 1.
d_sqtl <- list(
  snp     = liver_tf$variant_id,     # we'll align by position, so snp IDs can differ
  beta    = liver_tf$slope,
  varbeta = liver_tf$slope_se^2,
  MAF     = liver_tf$maf,
  type    = "quant",
  sdY     = 1
)

## ---- 4) Align datasets by position (safe cross-panel key) ----
# Create a position key for both (chr:pos)
Ferr_win[, key := paste(chr, pos)]
liver_tf[, key := paste(chr, pos)]

# Intersect keys
keys <- intersect(Ferr_win$key, liver_tf$key)
Ferr_aln  <- Ferr_win[key %in% keys]
sQTL_aln  <- liver_tf[key %in% keys]

# Order both by key to keep rows aligned
setorder(Ferr_aln, key); setorder(sQTL_aln, key)

# Rebuild coloc lists with aligned SNP sets (use Ferritin rsID as the 'snp' ID for readability)
d1 <- list(
  snp     = Ferr_aln$SNP,
  beta    = Ferr_aln$beta,
  varbeta = Ferr_aln$se^2,
  MAF     = Ferr_aln$maf,
  N       = Ferr_aln$samplesize,
  type    = "quant")

d2 <- list(
  snp     = Ferr_aln$SNP,           # use same order/IDs as d1 to be safe
  beta    = sQTL_aln$slope,
  varbeta = sQTL_aln$slope_se^2,
  MAF     = sQTL_aln$maf,
  type    = "quant",
  sdY     = 1)

## ---- 5) Run coloc.abf ----
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_tf_liver <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

print(cc_tf_liver$summary)

# Directionality check (Ferritin × TF Liver sQTL) 
# 1) Lead SNP by posterior contribution to H4
res <- as.data.table(cc_tf_liver$results)
if (!"SNP.PP.H4" %in% names(res)) stop("Expected column SNP.PP.H4 not found in coloc results")
setorder(res, -SNP.PP.H4)
lead <- res[1, snp]

# 2) Merge alleles/frequencies for overlapping SNPs
dir_df <- merge(
  Ferr_aln[, .(key, SNP, effect_allele, other_allele, eaf, beta_gwas = beta)],
  sQTL_aln[, .(key, ref, alt, af, slope)],
  by = "key"
)

# 3) Align sQTL slope to GWAS effect allele (Ferritin):
#    - If GWAS effect_allele == sQTL ALT -> same counted allele -> keep sign
#    - If GWAS effect_allele == sQTL REF -> opposite counted allele -> flip sQTL sign
dir_df[, slope_aligned := fifelse(effect_allele == alt,  slope,
                                  fifelse(effect_allele == ref, -slope, NA_real_))]
dir_df[, allele_match := fifelse(effect_allele == alt, "GWAS EA = ALT (no flip)",
                                 fifelse(effect_allele == ref, "GWAS EA = REF (flip sQTL)",
                                         "allele mismatch"))]

# Optional cross-check using frequencies (which allele the EAF refers to)
dir_df[, d_alt := abs(eaf - af)]
dir_df[, d_ref := abs(eaf - (1 - af))]
dir_df[, freq_support := fifelse(d_alt < d_ref, "EAF≈ALT", "EAF≈REF")]

# 4) Lead SNP summary (what allele does what?)
lead_row <- dir_df[SNP == lead]
cat("\n=== Lead SNP directionality ===\n")
print(lead_row[, .(SNP, effect_allele, other_allele,
                   ref, alt, eaf, af,
                   beta_gwas, slope, slope_aligned,
                   allele_match, freq_support)])

# 5) Global sign-concordance across SNPs
dir_df2 <- dir_df[!is.na(slope_aligned)]
dir_df2[, sign_agree := sign(beta_gwas) == sign(slope_aligned)]
agree_rate <- mean(dir_df2$sign_agree)

# Weighted by SNP colocalisation mass (PP.H4)
dir_df2 <- merge(dir_df2, res[, .(SNP = snp, w = SNP.PP.H4)], by = "SNP", all.x = TRUE)
dir_df2[is.na(w), w := 0]
w_agree  <- dir_df2[, sum(w * as.numeric(sign_agree))]
w_total  <- dir_df2[, sum(w)]
w_agree_pct <- ifelse(w_total > 0, 100 * w_agree / w_total, NA_real_)

cat(sprintf("\nSign concordance (unaligned count): %d/%d (%.1f%%)\n",
            sum(dir_df2$sign_agree), nrow(dir_df2), 100*agree_rate))
cat(sprintf("Sign concordance weighted by PP.H4 mass: %s%%\n",
            ifelse(is.na(w_agree_pct), "NA", sprintf("%.1f", w_agree_pct)))) 
#-------------------------------------------------------------------------------
# ---- 0) Define window ----
chr_ceac   <- 19L
start_ceac <- 44119993L
end_ceac   <- 45119993L

# ---- 1) Ferritin GWAS in window (reuse Ferr_dt) ----
Ferr_ceac <- copy(Ferr_dt)[chr == chr_ceac & pos >= start_ceac & pos <= end_ceac]
Ferr_ceac[, maf := pmin(eaf, 1 - eaf)]
Ferr_ceac <- Ferr_ceac[!is.na(beta) & !is.na(se) & !is.na(maf) & maf >= 0.005]
Ferr_ceac <- unique(Ferr_ceac, by = "SNP")

# ---- 2) Helper: robust parse of GTEx variant_id ("chr19_..._b38") ----
parse_variant_id <- function(dt, id_col = "variant_id") {
  stopifnot(id_col %in% names(dt))
  ss <- strsplit(dt[[id_col]], "_", fixed = TRUE)
  grab <- function(v, k) if (length(v) >= k) v[[k]] else NA_character_
  dt[, `:=`(
    chr_str = vapply(ss, grab, "", 1),
    pos_str = vapply(ss, grab, "", 2),
    ref     = vapply(ss, grab, "", 3),
    alt     = vapply(ss, grab, "", 4),
    build   = vapply(ss, grab, "", 5)
  )]
  dt[, chr := as.integer(sub("^chr","", chr_str))]
  dt[, pos := as.integer(pos_str)]
  dt[]
}

# ---- 3) Choose tissue: Cortex first; if empty, try BA9 ----
braincor_sig <- as.data.table(braincor_sig)
brainfro_sig <- as.data.table(brainfro_sig)

# CEACAM19 gene IDs
ceac_ctx_ids <- braincor_genes[gene_name == "CEACAM19", unique(gene_id)]
ceac_ba9_ids <- brainfro_genes[gene_name == "CEACAM19", unique(gene_id)]

# Cortex subset
sqtl_dt   <- braincor_sig[gene_id %in% ceac_ctx_ids]
genes_dt  <- braincor_genes
tissue_lbl <- "Brain Cortex"

if (nrow(sqtl_dt) == 0L) {
  # Fallback to BA9
  sqtl_dt    <- brainfro_sig[gene_id %in% ceac_ba9_ids]
  genes_dt   <- brainfro_genes
  tissue_lbl <- "Brain Frontal Cortex BA9"
}
if (nrow(sqtl_dt) == 0L) {
  stop("CEACAM19 has 0 significant sQTL pairs in both Cortex and BA9 'signif_pairs'. Use sQTL all-pairs/nominal or brain eQTL instead.")
}

# ---- 4) Parse, window, QC on sQTL side ----
sqtl_dt <- parse_variant_id(copy(sqtl_dt))                         # adds chr/pos/ref/alt
sqtl_dt <- sqtl_dt[chr == chr_ceac & pos >= start_ceac & pos <= end_ceac]
sqtl_dt <- sqtl_dt[nchar(ref) == 1 & nchar(alt) == 1]              # SNP-only
# Standardise effect columns (GTEx v10 sQTL gives slope/slope_se; but guard anyway)
if (!"slope"    %in% names(sqtl_dt) && "beta" %in% names(sqtl_dt)) sqtl_dt[, slope := beta]
if (!"slope_se" %in% names(sqtl_dt) && "se"   %in% names(sqtl_dt)) sqtl_dt[, slope_se := se]
# MAF
if (!"af" %in% names(sqtl_dt) && "maf" %in% names(sqtl_dt)) sqtl_dt[, af := maf]
sqtl_dt[, maf := pmin(af, 1 - af)]
sqtl_dt <- sqtl_dt[!is.na(slope) & !is.na(slope_se) & !is.na(maf) & maf >= 0.005]
sqtl_dt <- unique(sqtl_dt, by = "variant_id")

cat(sprintf("%s: CEACAM19 sQTL rows after QC in window: %d\n", tissue_lbl, nrow(sqtl_dt)))

# ---- 5) Align by position (key = chr:pos) ----
Ferr_ceac[, key := paste(chr, pos)]
sqtl_dt[,  key := paste(chr, pos)]
keys <- intersect(Ferr_ceac$key, sqtl_dt$key)

Ferr_aln <- Ferr_ceac[key %in% keys]; setorder(Ferr_aln, key)
sQTL_aln <- sqtl_dt[key %in% keys];   setorder(sQTL_aln, key)

cat(sprintf("Overlap SNPs after QC: %d\n", nrow(Ferr_aln)))
if (nrow(Ferr_aln) == 0L) stop("No overlapping SNPs after alignment. Try all-pairs sQTL or a broader window.")

# ---- 6) Build coloc datasets (quantitative) ----
d1 <- list(
  snp     = Ferr_aln$SNP,
  beta    = Ferr_aln$beta,
  varbeta = Ferr_aln$se^2,
  MAF     = Ferr_aln$maf,
  N       = Ferr_aln$samplesize,
  type    = "quant"
)
d2 <- list(
  snp     = Ferr_aln$SNP,                 # reuse GWAS rsIDs/order
  beta    = sQTL_aln$slope,
  varbeta = sQTL_aln$slope_se^2,
  MAF     = sQTL_aln$maf,
  type    = "quant",
  sdY     = 1
)

# ---- 7) coloc.abf ----
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_ceac <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)
cat(sprintf("\n== Ferritin × CEACAM19 sQTL (%s) ==\n", tissue_lbl))
print(cc_ceac$summary)

# ---- 8) Directionality sanity check (like liver) ----
res <- as.data.table(cc_ceac$results)
if (nrow(res)) {
  setorder(res, -SNP.PP.H4)
  lead <- res[1, snp]
  
  dir_df <- merge(
    Ferr_aln[, .(key, SNP, effect_allele, other_allele, eaf, beta_gwas = beta)],
    sQTL_aln[, .(key, ref, alt, af, slope)],
    by = "key"
  )
  dir_df[, slope_aligned := fifelse(effect_allele == alt, slope,
                                    fifelse(effect_allele == ref, -slope, NA_real_))]
  dir_df[, allele_match := fifelse(effect_allele == alt, "GWAS EA = ALT (no flip)",
                                   fifelse(effect_allele == ref, "GWAS EA = REF (flip sQTL)", "allele mismatch"))]
  dir_df[, d_alt := abs(eaf - af)]
  dir_df[, d_ref := abs(eaf - (1 - af))]
  dir_df[, freq_support := fifelse(d_alt < d_ref, "EAF≈ALT", "EAF≈REF")]
  
  cat("\n=== Lead SNP directionality ===\n")
  print(dir_df[SNP == lead,
               .(SNP, effect_allele, other_allele, ref, alt, eaf, af,
                 beta_gwas, slope, slope_aligned, allele_match, freq_support)][1])
  
  dir_df2 <- dir_df[!is.na(slope_aligned)]
  dir_df2[, sign_agree := sign(beta_gwas) == sign(slope_aligned)]
  agree_rate <- mean(dir_df2$sign_agree)
  
  dir_df2 <- merge(dir_df2, res[, .(SNP = snp, w = SNP.PP.H4)], by = "SNP", all.x = TRUE)
  dir_df2[is.na(w), w := 0]
  w_agree  <- dir_df2[, sum(w * as.numeric(sign_agree))]
  w_total  <- dir_df2[, sum(w)]
  w_agree_pct <- ifelse(w_total > 0, 100 * w_agree / w_total, NA_real_)
  
  cat(sprintf("\nSign concordance (unaligned count): %d/%d (%.1f%%)\n",
              sum(dir_df2$sign_agree), nrow(dir_df2), 100*agree_rate))
  cat(sprintf("Sign concordance weighted by PP.H4 mass: %s%%\n",
              ifelse(is.na(w_agree_pct), "NA", sprintf("%.1f", w_agree_pct))))
}


