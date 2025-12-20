library(data.table)
library(coloc) 
library(arrow)

del <- fread("Delirium_AF0p005.mr_ready.tsv.gz") 
pqtl  <- fread("prot-a-131.vcf.gz", sep = "\t", header = TRUE, skip = "#CHROM", showProgress = FALSE)
#-------------------------------------------------------------------------------
# Set case–control meta for delirium
n_cases    <- 8461
n_controls <- 449979
N_cc       <- n_cases + n_controls
s_cc       <- n_cases / N_cc

# APOE window (GRCh38)
chr_apo <- 19L; start_apo <- 44421094L; end_apo <- 45421094L

# Keep only chr19 APOE-window SNPs from Delirium (and clean to SNP-only)
del <- del[chr == chr_apo & pos >= start_apo & pos <= end_apo]
del <- del[nchar(effect_allele)==1 & nchar(other_allele)==1] 

## ==== Delirium (cc) × APOE pQTL (prot-a-131) — coloc ====

# 0) Small helpers
split_one <- function(x, sep=":") strsplit(x, sep, fixed=TRUE)[[1]]
info_tag  <- function(info, tag) {
  # Pull a single numeric tag (e.g. AF=) from INFO; returns NA if absent
  v <- sub(paste0(".*(?:^|;)", tag, "=([^;]+).*"), "\\1", info)
  v[ v == info ] <- NA_character_
  suppressWarnings(as.numeric(v))
}

# 1) Make pQTL tidy: rename "#CHROM" -> CHROM and parse FORMAT fields (ES, SE, AF)
setnames(pqtl, old = "#CHROM", new = "CHROM", skip_absent = TRUE)

# Check FORMAT is present and looks like "ES:SE:LP:AF:ID"
stopifnot("FORMAT" %in% names(pqtl), any(grepl("ES", pqtl$FORMAT, fixed=TRUE)))

fmt_keys <- unique(pqtl$FORMAT)
if (length(fmt_keys) != 1L) {
  # If multiple FORMAT strings appear, take the first but warn
  warning("Multiple FORMAT strings detected; using the first: ", fmt_keys[1])
}
keys <- split_one(fmt_keys[1])

# Split the sample column using FORMAT order
samp_col <- setdiff(names(pqtl), c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
stopifnot(length(samp_col) == 1L)  # should be "PROT-a-131"
parts <- tstrsplit(pqtl[[samp_col]], ":", fixed = TRUE, fill = NA_character_)
if (length(parts) != length(keys)) {
  stop("Unexpected number of fields in sample column vs FORMAT. Check the VCF.")
}
pqtl[, (keys) := parts]

# Coerce numeric where relevant
num_keys <- intersect(keys, c("ES","SE","AF","LP","SS"))
for (k in num_keys) pqtl[, (k) := as.numeric(get(k))]

# If AF missing in FORMAT, try INFO tag AF
if (!"AF" %in% keys || all(is.na(pqtl$AF))) {
  pqtl[, AF := info_tag(INFO, "AF")]
}

# Keep biallelic SNPs only
pqtl <- pqtl[nchar(REF) == 1 & nchar(ALT) == 1]

# 2) Restrict del to SNPs that have an rsID and to biallelic SNPs 
del <- del[grepl("^rs", SNP)]
del <- del[nchar(effect_allele)==1 & nchar(other_allele)==1]
del[, maf := pmin(eaf, 1 - eaf)]

# 3) Join by rsID, align pQTL beta to GWAS effect allele
ov <- merge(
  del[, .(SNP, beta_cc = beta, se_cc = se, eaf_cc = eaf,
          EA = effect_allele, OA = other_allele,
          N_cc_snp = samplesize)],
  pqtl[, .(ID, CHROM, POS, REF, ALT, ES, SE, AF)],
  by.x = "SNP", by.y = "ID"
)

# Drop rows lacking ES/SE/AF
ov <- ov[!is.na(ES) & !is.na(SE) & !is.na(AF)]

# Align pQTL effect (ES is per ALT allele). If GWAS effect allele is REF, flip pQTL sign.
ov[, beta_pqtl := fifelse(EA == ALT, ES,
                          fifelse(EA == REF, -ES, NA_real_))]
ov[, allele_ok := (EA == ALT) | (EA == REF)]
ov <- ov[allele_ok & !is.na(beta_pqtl)]

# Minor allele frequencies
ov[, maf_cc   := pmin(eaf_cc, 1 - eaf_cc)]
ov[, maf_pqtl := pmin(AF,      1 - AF)]

# 4) Build coloc datasets
# Case–control (Delirium)
d1 <- list(
  snp     = ov$SNP,
  beta    = ov$beta_cc,
  varbeta = ov$se_cc^2,
  MAF     = ov$maf_cc,
  N       = if ("N_cc_snp" %in% names(ov) && all(!is.na(ov$N_cc_snp))) ov$N_cc_snp else N_cc,
  type    = "cc",
  s       = s_cc
)

# Quantitative (APOE pQTL)
d2 <- list(
  snp     = ov$SNP,
  beta    = ov$beta_pqtl,
  varbeta = ov$SE^2,
  MAF     = ov$maf_pqtl,
  type    = "quant",
  sdY     = 1
)

# 5) Run coloc
priors <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
cc <- coloc.abf(d1, d2, p1 = priors$p1, p2 = priors$p2, p12 = priors$p12)

cat("\n=== Delirium × APOE pQTL — coloc.abf SUMMARY ===\n")
print(cc$summary)
cat(sprintf("PP.H4 (shared causal): %.1f%%\n", 100 * cc$summary["PP.H4.abf"]))

# 6) Top SNPs by PP.H4 (if available)
res <- as.data.table(cc$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res) > 0) {
  data.table::setorder(res, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4:\n")
  print(res[1:min(.N, 10), .(snp, PP.H4 = SNP.PP.H4)])
}

