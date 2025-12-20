library(data.table)
library(coloc) 

Ferr <- fread("Ferritin_AF0p005.mr_ready.tsv.gz") 
pqtl  <- fread("prot-a-131.vcf.gz", sep = "\t", header = TRUE, skip = "#CHROM", showProgress = FALSE)
setnames(pqtl, sub("^#", "", names(pqtl))) 

# ---- Define APOE window (GRCh38) ----
chr_apo   <- 19L
start_apo <- 44409976L
end_apo   <- 45409976L

# ---- Subset Ferritin to window + prep ----
Ferr_dt <- as.data.table(Ferr)
# Ensure numeric types
numcols <- intersect(c("beta","se","eaf","samplesize","chr","pos"), names(Ferr_dt))
for (nm in numcols) Ferr_dt[, (nm) := as.numeric(get(nm))]

Ferr_win <- Ferr_dt[chr == chr_apo & pos >= start_apo & pos <= end_apo]
Ferr_win[, maf := pmin(eaf, 1 - eaf)]
Ferr_win <- Ferr_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1]
Ferr_win <- unique(Ferr_win, by = "SNP")  # one row per rsID

# ---- Parse pQTL sample column (FORMAT = ES:SE:LP:AF:ID) ----
# Identify the single sample column (everything that's not a core VCF field)
core_cols   <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
sample_cols <- setdiff(names(pqtl), core_cols)
stopifnot(length(sample_cols) == 1)  # should be exactly one: "PROT-a-131"
sample_col  <- sample_cols[1]

# Split FORMAT keys once 
fmt_keys <- strsplit(pqtl$FORMAT[1], ":", fixed = TRUE)[[1]]    # e.g., c("ES","SE","LP","AF","ID")

# Split the sample field into columns by ':'
parts <- tstrsplit(pqtl[[sample_col]], ":", fixed = TRUE)
stopifnot(length(parts) == length(fmt_keys))
for (i in seq_along(fmt_keys)) {
  # Coerce numeric where appropriate; 'ID' stays character
  if (fmt_keys[i] %in% c("ES","SE","LP","AF")) {
    pqtl[, (fmt_keys[i]) := as.numeric(parts[[i]])]
  } else {
    pqtl[, (fmt_keys[i]) := parts[[i]]]
  }
}

# Keep biallelic SNPs with sane alleles
pqtl <- pqtl[
  nchar(REF) == 1 & nchar(ALT) == 1 &
    REF %in% c("A","C","G","T") & ALT %in% c("A","C","G","T")
]

# Compute MAF
pqtl[, maf := pmin(AF, 1 - AF)]
pqtl <- pqtl[!is.na(ES) & !is.na(SE) & !is.na(maf) & maf > 0 & maf < 1]

# ---- Cross-build alignment by rsID: intersect with Ferritin window rsIDs ----
ids <- intersect(Ferr_win$SNP, pqtl$ID)
if (length(ids) == 0L) stop("No overlapping rsIDs between Ferritin window and prot-a-131. Widen window or check files.")

Ferr_aln <- Ferr_win[SNP %in% ids]
pqtl_aln <- pqtl[ID %in% ids]

# Order both identically by rsID
setorder(Ferr_aln, SNP)
setorder(pqtl_aln, ID)
stopifnot(identical(Ferr_aln$SNP, pqtl_aln$ID))

# ---- Build COLOC datasets ----
# d1: Ferritin (quantitative)
d1 <- list(
  snp     = Ferr_aln$SNP,
  beta    = Ferr_aln$beta,
  varbeta = Ferr_aln$se^2,
  MAF     = Ferr_aln$maf,
  N       = Ferr_aln$samplesize,   # per-SNP is fine; coloc accepts vector N
  type    = "quant"
)

# d2: APOE pQTL (quantitative) — use sdY=1 to avoid needing N
d2 <- list(
  snp     = pqtl_aln$ID,
  beta    = pqtl_aln$ES,           # effect relative to ALT allele
  varbeta = pqtl_aln$SE^2,
  MAF     = pqtl_aln$maf,
  type    = "quant",
  sdY     = 1
)

# ---- Run COLOC ----
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

cat("\n=== COLOC summary (Ferritin × APOE pQTL) ===\n")
print(cc$summary)

# ---- Top SNPs by PP.H4 ----
res <- as.data.table(cc$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res)) {
  setorder(res, -SNP.PP.H4)
  cat("\nTop 10 SNPs by PP.H4:\n")
  print(res[1:min(10, .N), .(snp, PP.H4 = SNP.PP.H4)])
}

# ---- Directionality sanity check (optional, recommended) ----
# pQTL ES is per ALT allele. Align its sign to Ferritin effect_allele.
dir_df <- merge(
  Ferr_aln[, .(SNP, effect_allele, other_allele, eaf, beta_gwas = beta)],
  pqtl_aln[, .(ID, REF, ALT, AF, ES)],
  by.x = "SNP", by.y = "ID", all = FALSE
)

dir_df[, slope_aligned := fifelse(effect_allele == ALT,  ES,
                                  fifelse(effect_allele == REF, -ES, NA_real_))]
dir_df[, allele_match := fifelse(effect_allele == ALT, "GWAS EA = ALT (no flip)",
                                 fifelse(effect_allele == REF, "GWAS EA = REF (flip pQTL)",
                                         "allele mismatch"))]

# Quick summary
keep <- dir_df[!is.na(slope_aligned)]
keep[, sign_agree := sign(beta_gwas) == sign(slope_aligned)]
agree_rate <- mean(keep$sign_agree)
cat(sprintf("\nSign concordance (count): %d/%d (%.1f%%)\n",
            sum(keep$sign_agree), nrow(keep), 100*agree_rate))

# Weight by PP.H4
if ("SNP.PP.H4" %in% names(res)) {
  keep <- merge(keep, res[, .(SNP = snp, w = SNP.PP.H4)], by = "SNP", all.x = TRUE)
  keep[is.na(w), w := 0]
  w_agree <- keep[, sum(w * as.numeric(sign_agree))]
  w_tot   <- keep[, sum(w)]
  w_pct   <- ifelse(w_tot > 0, 100 * w_agree / w_tot, NA_real_)
  cat(sprintf("Sign concordance (PP.H4-weighted): %s%%\n",
              ifelse(is.na(w_pct), "NA", sprintf('%.1f', w_pct))))
}

# Lead SNP summary for logs
if (exists("res") && nrow(res)) {
  lead <- res[1, snp]
  lead_row <- keep[SNP == lead]
  if (nrow(lead_row)) {
    cat("\nLead SNP directionality (effect is per Ferritin GWAS effect allele):\n")
    print(lead_row[, .(SNP, effect_allele, other_allele,
                       REF, ALT, eaf, AF,
                       beta_gwas, ES, slope_aligned, allele_match)])
  }
}
