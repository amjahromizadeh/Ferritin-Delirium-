library(data.table)
library(coloc)
library(knitr)

## files
del       <- fread("Delirium_AF0p005.mr_ready.tsv.gz")
ctx_sig   <- fread("Brain_Cortex.v8.signif_variant_gene_pairs.txt.gz")
ctx_genes <- fread("Brain_Cortex.v8.egenes.txt.gz")

## Region (GRCh38): CEACAM19 ±500 kb  -> chr19:44,162,645–45,162,645
chr_ce   <- 19L
start_ce <- 44162645L
end_ce   <- 45162645L

## Delirium GWAS: types, slice window, 1 row / position
gwas <- as.data.table(del)

# chr parsing 
if (!"chr" %in% names(gwas)) stop("No 'chr' column in Delirium file.")
gwas[, chr := sub("^chr", "", as.character(chr))]
gwas[, chr := as.integer(chr)]

# coerce numerics
for (nm in intersect(c("beta","se","eaf","pval","samplesize","pos"), names(gwas))) {
  gwas[, (nm) := as.numeric(get(nm))]
}

# window + MAF filter (>=0.005 as per your spec)
gwas_win <- gwas[chr == chr_ce & pos >= start_ce & pos <= end_ce]
gwas_win[, maf := pmin(eaf, 1 - eaf)]
gwas_win <- gwas_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005]

# one row per position: pick strongest (smallest p if available, else smallest se)
if ("pval" %in% names(gwas_win)) {
  gwas_pos <- gwas_win[order(pval)][, .SD[1], by = .(chr, pos)]
} else {
  gwas_pos <- gwas_win[order(se)][,   .SD[1], by = .(chr, pos)]
}

cat(sprintf("GWAS in window: %d rows -> %d unique positions\n", nrow(gwas_win), nrow(gwas_pos)))

## GTEx v8 Brain Cortex eQTL: keep CEACAM19 rows, parse positions, QC
# find gene_id(s) for CEACAM19
ce_ids <- ctx_genes[gene_name == "CEACAM19", unique(gene_id)]
if (length(ce_ids) == 0) stop("CEACAM19 not found in Brain_Cortex.v8.egenes.txt.gz")

eqtl_all <- as.data.table(ctx_sig)[gene_id %in% ce_ids]

# parse 'variant_id' like 'chr19_44910503_A_G_b38'
if (!"variant_id" %in% names(eqtl_all)) stop("No 'variant_id' column in v8 signif pairs file.")
parts <- tstrsplit(eqtl_all$variant_id, "_", fixed = TRUE)
if (length(parts) < 5) stop("Unexpected variant_id format in v8 signif pairs.")

eqtl_all[, `:=`(
  chr = as.integer(sub("^chr", "", parts[[1]])),
  pos = as.integer(parts[[2]]),
  ref = parts[[3]],
  alt = parts[[4]])
]

# slopes/SE columns (v8 uses slope/slope_se)
if (!"slope" %in% names(eqtl_all) && "beta" %in% names(eqtl_all))    setnames(eqtl_all, "beta", "slope")
if (!"slope_se" %in% names(eqtl_all) && "se"   %in% names(eqtl_all)) setnames(eqtl_all, "se",   "slope_se")

# MAF on eQTL side: prefer 'maf'; else derive from 'af' if present
if ("maf" %in% names(eqtl_all)) {
  eqtl_all[, maf := pmin(as.numeric(maf), 1 - as.numeric(maf))]
} else if ("af" %in% names(eqtl_all)) {
  eqtl_all[, maf := pmin(as.numeric(af), 1 - as.numeric(af))]
} else {
  stop("No 'maf' or 'af' column in Brain Cortex v8 signif pairs for CEACAM19.")
}

# window + QC: biallelic SNPs, complete stats, MAF ≥ 0.005
eqtl_win <- eqtl_all[
  chr == chr_ce & pos >= start_ce & pos <= end_ce &
    !is.na(slope) & !is.na(slope_se) & !is.na(maf) & maf > 0 & maf < 1 & maf >= 0.005 &
    nchar(ref) == 1 & nchar(alt) == 1
]

# one row per position: prefer smallest pval_nominal if available
if ("pval_nominal" %in% names(eqtl_win)) {
  eqtl_pos <- eqtl_win[order(pval_nominal)][, .SD[1], by = .(chr, pos)]
} else {
  eqtl_pos <- eqtl_win[order(slope_se)][,       .SD[1], by = .(chr, pos)]
}

cat(sprintf("CEACAM19 eQTL in window: %d rows -> %d unique positions\n", nrow(eqtl_win), nrow(eqtl_pos)))

## Align by genomic position
gwas_pos[, key := paste(chr, pos, sep=":")]
eqtl_pos[, key := paste(chr, pos, sep=":")]

aln <- merge(gwas_pos, eqtl_pos, by = c("chr","pos","key"))
cat(sprintf("Overlapping positions (chr:pos) = %d\n", nrow(aln)))
if (nrow(aln) == 0) {
  cat("Example GWAS keys:\n"); print(head(gwas_pos$key, 5))
  cat("Example eQTL keys:\n"); print(head(eqtl_pos$key, 5))
  stop("No overlapping positions after cleaning. Check build or allele filters.")
}

## Build coloc datasets (same SNP vector/order)
# Delirium case/control
N_cases <- 8461
N_ctrls <- 449979
d1 <- list(
  snp     = aln$key,
  beta    = aln$beta,         # GWAS beta (on effect_allele)
  varbeta = aln$se^2,
  MAF     = aln$maf.x,        # from GWAS
  type    = "cc",
  s       = N_cases / (N_cases + N_ctrls),
  N       = N_cases + N_ctrls
)

# CEACAM19 expression (quantitative)
d2 <- list(
  snp     = aln$key,
  beta    = aln$slope,        # eQTL slope
  varbeta = aln$slope_se^2,
  MAF     = aln$maf.y,        # from eQTL
  type    = "quant",
  sdY     = 1
)

## Run coloc
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_ctx <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

cat("\n=== Delirium × CEACAM19 — GTEx v8 eQTL (Brain Cortex) — coloc.abf ===\n")
print(cc_ctx$summary)
cat(sprintf("PP.H4 (shared causal): %.3f%%\n", 100 * cc_ctx$summary["PP.H4.abf"]))

# Optional: top SNPs by PP.H4 if available
res <- as.data.table(cc_ctx$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res) > 0) {
  data.table::setorder(res, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4:\n")
  print(res[1:min(.N, 10), .(snp = snp, PP.H4 = SNP.PP.H4)])
}

## ==== Directionality sanity check: Delirium GWAS vs CEACAM19 eQTL (Cortex v8) ====

stopifnot(exists("aln"), exists("cc_ctx"))
if (!exists("res")) res <- data.table::as.data.table(cc_ctx$results)

# 1) Build an allele-aware table for overlapping SNPs
dir_df <- data.table::data.table(
  SNP           = aln$key,                    # "chr:pos"
  effect_allele = aln$effect_allele,          # GWAS effect allele
  other_allele  = aln$other_allele,           # GWAS non-effect
  ref           = aln$ref,                    # eQTL REF
  alt           = aln$alt,                    # eQTL ALT
  eaf           = aln$eaf,                    # GWAS effect allele freq (if available)
  beta_gwas     = aln$beta,                   # GWAS effect (per GWAS EA)
  slope         = aln$slope,                  # eQTL effect (per ALT by GTEx convention)
  af_alt        = if ("af" %in% names(aln)) aln$af else NA_real_  # ALT freq if present
)

# 2) Align eQTL slope to the GWAS effect allele
dir_df[, slope_aligned := fifelse(effect_allele == alt,  slope,
                                  fifelse(effect_allele == ref, -slope, NA_real_))]
dir_df[, allele_match := fifelse(effect_allele == alt, "GWAS EA = ALT (no flip)",
                                 fifelse(effect_allele == ref, "GWAS EA = REF (flip eQTL)", 
                                         "allele mismatch"))]

# frequency cross-check when ALT freq is available
if (!all(is.na(dir_df$af_alt))) {
  dir_df[, d_alt := abs(eaf - af_alt)]
  dir_df[, d_ref := abs(eaf - (1 - af_alt))]
  dir_df[, freq_support := fifelse(d_alt < d_ref, "EAF≈ALT", "EAF≈REF")]
} else {
  dir_df[, freq_support := NA_character_]
}

# 3) Lead SNP by per-SNP H4
data.table::setorder(res, -SNP.PP.H4)
lead <- res[1, snp]
lead_row <- merge(
  dir_df[SNP == lead],
  res[, .(SNP = snp, PP.H4 = SNP.PP.H4)],
  by = "SNP", all.x = TRUE
)

# 4) Global sign concordance
keep <- dir_df[!is.na(slope_aligned)]
keep[, sign_agree := sign(beta_gwas) == sign(slope_aligned)]
agree_rate <- mean(keep$sign_agree)

# PP.H4-weighted concordance
keep <- merge(keep, res[, .(SNP = snp, w = SNP.PP.H4)], by = "SNP", all.x = TRUE)
keep[is.na(w), w := 0]
w_agree <- keep[, sum(w * as.numeric(sign_agree))]
w_tot   <- keep[, sum(w)]
w_pct   <- ifelse(w_tot > 0, 100 * w_agree / w_tot, NA_real_)

# 5) Print a compact summary + a clean lead-SNP table for results log
cat(sprintf("\nSign concordance (count): %d/%d (%.1f%%)\n",
            sum(keep$sign_agree), nrow(keep), 100*agree_rate))
cat(sprintf("Sign concordance (PP.H4-weighted): %s%%\n",
            ifelse(is.na(w_pct), "NA", sprintf('%.1f', w_pct))))

# Pretty table for results
lead_tbl <- lead_row[, .(SNP, effect_allele, other_allele, ref, alt,
                         eaf, af_alt, beta_gwas, slope, slope_aligned,
                         allele_match, freq_support, PP.H4)]
knitr::kable(lead_tbl, caption = "Lead SNP directionality — Delirium × CEACAM19 (Cortex v8)")  
#-------------------------------------------------------------------------------
