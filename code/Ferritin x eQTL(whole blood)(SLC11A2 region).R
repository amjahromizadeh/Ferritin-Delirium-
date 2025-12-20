library(data.table) 
library(coloc) 

Ferr <- fread("Ferritin_AF0p005.mr_ready.tsv.gz")

# Opening primary panel 
primary <- fread("2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

# Opening the proxy panel 
proxy1 <- fread("GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
proxy2 <- fread("GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt.gz")
#-------------------------------------------------------------------------------
# ---- helpers to be robust to column naming in eQTLGen ----
first_present <- function(nm, choices) {
  hit <- choices[choices %in% names(nm)]
  if (length(hit) == 0) stop("None of these columns found: ", paste(choices, collapse=", "))
  hit[1]
}

# 1) Define GRCh38 SLC11A2 window (±250 kb)
chr38   <- 12L
start38 <- 50725057L   # chr12:50,725,057
end38   <- 51225057L   # chr12:51,225,057

# 2) Prep Ferritin GWAS in window (GRCh38)
Ferr_dt <- as.data.table(Ferr)
# make sure numeric
for (nm in intersect(c("beta","se","eaf","samplesize","chr","pos"), names(Ferr_dt))) {
  Ferr_dt[, (nm) := as.numeric(get(nm))]
}
Ferr_win <- Ferr_dt[chr == chr38 & pos >= start38 & pos <= end38]
Ferr_win[, maf := pmin(eaf, 1 - eaf)]
Ferr_win <- Ferr_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf >= 0.005]
Ferr_win <- unique(Ferr_win, by = "SNP")

if (nrow(Ferr_win) < 50) {
  message("Note: very few GWAS SNPs in the window; coloc may be underpowered.")
}

# 3) Subset eQTLGen to SLC11A2 and keep only columns we need
# Try common gene-symbol columns in eQTLGen:
gene_col <- first_present(primary, c("GeneSymbol","HGNCName","Gene","Symbol","gene","gene_symbol"))
snp_col  <- first_present(primary, c("SNP","rsid","RSID","variant","Variant"))
p_col    <- first_present(primary, c("PValue","Pvalue","pvalue","p","P"))
# Sample-size column (row-level); if absent we’ll fallback to a single N
n_col    <- if (any(c("NrSamples","N","n") %in% names(primary))) {
  first_present(primary, c("NrSamples","N","n"))
} else NA_character_

eqtl_gene <- primary[get(gene_col) == "SLC11A2", 
                     .SD, .SDcols = unique(na.omit(c(snp_col, p_col, n_col)))]

setnames(eqtl_gene, old = snp_col, new = "SNP")
setnames(eqtl_gene, old = p_col,   new = "pval_eQTL")
if (!is.na(n_col)) setnames(eqtl_gene, old = n_col, new = "N_eQTL")

if (nrow(eqtl_gene) == 0L) {
  stop("eQTLGen primary file has 0 rows for SLC11A2. ",
       "This is expected in some cases (no probe or no passing hits in this release).")
}

# 4) Align by rsID within the biological region:
# Intersect eQTLGen SNPs with the GWAS SNPs that lie inside the GRCh38 window.
eqtl_gene <- unique(eqtl_gene, by = "SNP")
Ferr_aln  <- Ferr_win[SNP %in% eqtl_gene$SNP]
eqtl_aln  <- eqtl_gene[SNP %in% Ferr_aln$SNP]

# Re-match to ensure identical order
setkey(Ferr_aln, SNP); setkey(eqtl_aln, SNP)
stopifnot(identical(Ferr_aln$SNP, eqtl_aln$SNP))

if (nrow(Ferr_aln) < 100) {
  message("Note: only ", nrow(Ferr_aln), " overlapping SNPs by rsID; coloc may be noisy.")
}

# 5) Build coloc datasets
N_gwas <- if ("samplesize" %in% names(Ferr_aln) && all(!is.na(Ferr_aln$samplesize))) Ferr_aln$samplesize else max(Ferr_dt$samplesize, na.rm=TRUE)

d1 <- list(
  snp     = Ferr_aln$SNP,
  beta    = Ferr_aln$beta,
  varbeta = Ferr_aln$se^2,
  MAF     = Ferr_aln$maf,
  N       = N_gwas,
  type    = "quant"
)

# eQTLGen
N_eqtl_scalar <- if ("N_eQTL" %in% names(eqtl_aln) && all(!is.na(eqtl_aln$N_eQTL))) {
  NULL
} else {
  # conservative scalar fallback — if you know the exact meta N, replace here
  30000L
}

d2 <- list(
  snp     = eqtl_aln$SNP,
  pvalues = eqtl_aln$pval_eQTL,
  MAF     = Ferr_aln$maf,                 # reuse GWAS MAF as proxy (EUR ancestry match)
  N       = if (is.null(N_eqtl_scalar))   eqtl_aln$N_eQTL else N_eqtl_scalar,
  type    = "quant"
)

# 6) Run coloc
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc_slc11a2_eqtlgen <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

# 7) Summarise
print(cc_slc11a2_eqtlgen$summary)

# lead SNP and quick sanity
res <- as.data.table(cc_slc11a2_eqtlgen$results)

# Inspect available columns (optional, for your log)
print(names(res))

# Sort by per-SNP H4 weight
stopifnot("SNP.PP.H4" %in% names(res))
setorder(res, -SNP.PP.H4)

cat("\nLead SNP by PP.H4:\n")
lead_row <- res[1, .(snp, PP.H4 = SNP.PP.H4)]
print(lead_row)

cat("\nTop 10 SNPs by PP.H4:\n")
print(res[1:10, .(snp, PP.H4 = SNP.PP.H4)])
#-------------------------------------------------------------------------------
# Directionality sanity check 
gene_symbol <- "SLC11A2"

# 1) SNPs actually used by coloc (from its results table)
res <- as.data.table(cc_slc11a2_eqtlgen$results)
stopifnot("snp" %in% names(res))
snps_used <- unique(res$snp)

# 2) Slice Ferritin at those SNPs and rename for clarity
gw <- Ferr[SNP %in% snps_used,
           .(SNP, effect_allele, other_allele, eaf, beta_gwas = beta)]

# 3) Pull SLC11A2 rows from eQTLGen and standardize columns
#    (This block auto-detects key columns across common eQTLGen formats.)
eqtl <- copy(primary)

# Identify the gene column
gene_col <- intersect(c("GeneSymbol","Gene","hgnc_symbol","gene"), names(eqtl))[1]
if (is.na(gene_col)) stop("Could not find a gene column in eQTLGen table.")

eqtl <- eqtl[get(gene_col) == gene_symbol]

# Identify SNP and allele columns
snp_col  <- intersect(c("SNP","rsid","variant","MarkerName"), names(eqtl))[1]
EA_col   <- intersect(c("AssessedAllele","Assessed_Allele","EA","effect_allele","AlleleAssessed","A1"), names(eqtl))[1]
OA_col   <- intersect(c("OtherAllele","NEA","non_effect_allele","A2","other_allele"), names(eqtl))[1]
beta_col <- intersect(c("beta","Beta","slope","NES","Effect","effect_size"), names(eqtl))[1]
z_col    <- intersect(c("Zscore","Z","z"), names(eqtl))[1]
af_col   <- intersect(c("AF","EAF","MAF","af","maf"), names(eqtl))[1]

if (is.na(snp_col))  stop("No SNP/rsid column found in eQTLGen table.")
if (is.na(EA_col))   stop("No effect/assessed allele column (EA) found in eQTLGen table.")
# OA can be missing in some exports; we handle that below.

# Choose an eQTL effect measure for sign (prefer beta/slope, else Z)
eff_vec <- NULL
eff_name <- NA_character_
if (!is.na(beta_col)) { eff_vec <- eqtl[[beta_col]]; eff_name <- beta_col }
if (is.null(eff_vec) && !is.na(z_col)) { eff_vec <- eqtl[[z_col]]; eff_name <- z_col }
if (is.null(eff_vec)) stop("Neither beta/slope nor Zscore found in eQTLGen table.")

eqtl_std <- data.table(
  SNP = eqtl[[snp_col]],
  EA  = toupper(eqtl[[EA_col]]),
  OA  = if (!is.na(OA_col)) toupper(eqtl[[OA_col]]) else NA_character_,
  eff_eqtl_raw = as.numeric(eff_vec),
  af  = if (!is.na(af_col)) as.numeric(eqtl[[af_col]]) else NA_real_
)

# Keep only the coloc SNPs
eqtl_std <- eqtl_std[SNP %in% snps_used]

# 4) Align eQTL effect to the GWAS effect allele (Ferritin)
dir_df <- merge(gw, eqtl_std, by = "SNP", all = FALSE)

# If OA is missing in eQTLGen, try to infer flip using GWAS other_allele
# (We only flip when we see a certain allele match; else set NA)
dir_df[, slope_aligned := fifelse(toupper(effect_allele) == EA,  eff_eqtl_raw,
                                  fifelse(!is.na(OA) & toupper(effect_allele) == OA, -eff_eqtl_raw,
                                          fifelse(!is.na(OA) & toupper(other_allele)  == EA, -eff_eqtl_raw,
                                                  fifelse(!is.na(OA) & toupper(other_allele)  == OA,  eff_eqtl_raw, NA_real_))))]

dir_df[, allele_match := fifelse(toupper(effect_allele) == EA, "GWAS EA = eQTL EA (no flip)",
                                 fifelse(!is.na(OA) & toupper(effect_allele) == OA, "GWAS EA = eQTL OA (flip)",
                                         fifelse(!is.na(OA) & toupper(other_allele)  == EA, "GWAS OA = eQTL EA (flip)",
                                                 fifelse(!is.na(OA) & toupper(other_allele)  == OA, "GWAS OA = eQTL OA (no flip)",
                                                         "allele mismatch/unknown"))))]

# Optional frequency sanity hint (if AF available): which allele is close to AF?
if ("af" %in% names(dir_df)) {
  dir_df[, freq_support := {
    ea_diff <- abs(eaf - af)                 # if eQTLGen AF is for EA
    oa_diff <- abs(eaf - (1 - af))           # if eQTLGen AF is for OA
    fifelse(ea_diff < oa_diff, "EAF≈EA_AF", "EAF≈OA_AF")
  }]
}

# 5) Sign concordance
keep <- dir_df[!is.na(slope_aligned)]
keep[, sign_agree := sign(beta_gwas) == sign(slope_aligned)]
agree_rate <- mean(keep$sign_agree)

# Weighted by SNP.PP.H4 from coloc
wtab <- res[, .(SNP = snp, w = SNP.PP.H4)]
keep <- merge(keep, wtab, by = "SNP", all.x = TRUE)
keep[is.na(w), w := 0]
w_agree <- keep[, sum(w * as.numeric(sign_agree))]
w_tot   <- keep[, sum(w)]
w_pct   <- ifelse(w_tot > 0, 100 * w_agree / w_tot, NA_real_)

cat(sprintf("\nSign concordance (count): %d/%d (%.1f%%)\n",
            sum(keep$sign_agree), nrow(keep), 100*agree_rate))
cat(sprintf("Sign concordance (PP.H4-weighted): %s%%\n",
            ifelse(is.na(w_pct), "NA", sprintf('%.1f', w_pct))))

# 6) Lead SNP inspection (by per-SNP H4 mass)
setorder(res, -SNP.PP.H4)
lead <- res[1, snp]
lead_row <- keep[SNP == lead,
                 .(SNP, effect_allele, other_allele,
                   EA, OA, eaf, af,
                   beta_gwas, eff_eqtl_raw,
                   slope_aligned, allele_match)]
cat("\n== Lead SNP (by SNP.PP.H4) ==\n")
print(lead_row)
