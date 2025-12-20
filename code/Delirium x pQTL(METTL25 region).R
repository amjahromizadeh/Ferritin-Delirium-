library(data.table)
library(coloc)

## 0) Inputs
del <- fread("Delirium_AF0p005.mr_ready.tsv.gz")  
pqtl_path <- "prot-a-2958.vcf.gz"                 

chr_met   <- 17L          # <-- REPLACE if METTL25 is on a different chr in GRCh38
start_met <- 42900000L    # <-- REPLACE start (e.g., TSS-250kb)
end_met   <- 43400000L    # <-- REPLACE end  (e.g., TES+250kb)

## 2) Delirium GWAS: slice window, QC, one row per position
gwas <- as.data.table(del)
stopifnot("chr" %in% names(gwas), "pos" %in% names(gwas))
gwas[, chr := as.integer(sub("^chr","", as.character(chr)))]
for (nm in intersect(c("beta","se","eaf","pval","samplesize"), names(gwas))) gwas[, (nm) := as.numeric(get(nm))]

gwas_win <- gwas[chr == chr_met & pos >= start_met & pos <= end_met]
gwas_win[, maf := pmin(eaf, 1 - eaf)]
gwas_win <- gwas_win[!is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1]
gwas_pos <- if ("pval" %in% names(gwas_win)) gwas_win[order(pval)][, .SD[1], by=.(chr,pos)] else gwas_win[order(se)][, .SD[1], by=.(chr,pos)]
cat(sprintf("GWAS in METTL25 window: %d rows -> %d unique positions\n", nrow(gwas_win), nrow(gwas_pos)))

## 3) Read SomaLogic pQTL VCF 
# We read past header meta-lines by skipping to '#CHROM'
pqtl <- fread(pqtl_path, sep="\t", header=TRUE, skip="#CHROM", showProgress=FALSE)
stopifnot(all(c("#CHROM","POS","REF","ALT","FILTER","INFO","FORMAT") %in% names(pqtl)))
sample_col <- setdiff(names(pqtl), c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
stopifnot(length(sample_col) == 1)  # should be exactly one: e.g., PROT-a-2958

# Parse per-variant FORMAT fields (ES:SE:LP:AF:ID) from the sample column
fmt <- tstrsplit(pqtl[[sample_col]], ":", fixed=TRUE)
names(fmt) <- c("ES","SE","LP","AF","RSID")  # typical SomaLogic order
pq <- data.table(
  chr  = as.integer(pqtl[[ "#CHROM" ]]),
  pos  = as.integer(pqtl[[ "POS" ]]),
  ref  = pqtl[["REF"]],
  alt  = pqtl[["ALT"]],
  pass = pqtl[["FILTER"]] == "PASS",
  ES   = as.numeric(fmt$ES),
  SE   = as.numeric(fmt$SE),
  LP   = as.numeric(fmt$LP),
  AF   = as.numeric(fmt$AF),
  RSID = fmt$RSID
)

# QC + derive MAF (from AF) and p-value from -log10(p)
pq[, maf := pmin(AF, 1 - AF)]
pq[, p := 10^(-LP)]
pq <- pq[pass == TRUE & !is.na(ES) & !is.na(SE) & !is.na(maf) & maf > 0 & maf < 1 & nchar(ref) == 1 & nchar(alt) == 1]
pq <- pq[!is.na(RSID) & RSID != "."]

## 4) Align by rsID (avoids liftover issues)
# Keep only GWAS SNPs in the METTL25 window, then join by rsID
stopifnot("SNP" %in% names(gwas_pos))
aln <- merge(
  gwas_pos[, .(SNP, chr, pos, effect_allele, other_allele, eaf, beta, se, maf)],
  pq[, .(SNP = RSID, ES, SE, AF, maf_p = maf)],
  by = "SNP", all = FALSE
)
cat(sprintf("Overlapping rsIDs (GWAS × %s): %d\n", sample_col, nrow(aln)))
if (nrow(aln) < 10) {
  warning("Very few overlaps; coloc may be unstable. Consider trying the other pQTL VCF (e.g., prot-a-2959).")
}

## 5) Build coloc datasets (
# Delirium GWAS case-control
N_cases <- 8461
N_ctrls <- 449979

d1 <- list(
  snp     = aln$SNP,
  beta    = aln$beta,
  varbeta = aln$se^2,
  MAF     = aln$maf,
  type    = "cc",
  s       = N_cases / (N_cases + N_ctrls),
  N       = N_cases + N_ctrls
)
# pQTL 
d2 <- list(
  snp     = aln$SNP,
  beta    = aln$ES,        # ES is effect size on ALT allele
  varbeta = aln$SE^2,
  MAF     = aln$maf_p,
  type    = "quant",
  sdY     = 1
)

## 6) Run coloc
priors <- list(p1=1e-4, p2=1e-4, p12=1e-5)
cc_mettl25_pqtl <- coloc.abf(d1, d2, p1=priors$p1, p2=priors$p2, p12=priors$p12)

cat("\n=== Delirium × METTL25 — SomaLogic pQTL coloc ===\n")
print(cc_mettl25_pqtl$summary)
cat(sprintf("PP.H4 (shared causal): %.3f%%\n", 100 * cc_mettl25_pqtl$summary["PP.H4.abf"]))

# Optional: show top SNPs by per-SNP PP.H4 (if present)
res <- as.data.table(cc_mettl25_pqtl$results)
if ("SNP.PP.H4" %in% names(res) && nrow(res)) {
  data.table::setorder(res, -SNP.PP.H4)
  cat("\nTop SNPs by PP.H4:\n")
  print(res[1:min(.N, 10), .(snp, PP.H4 = SNP.PP.H4)])
}
