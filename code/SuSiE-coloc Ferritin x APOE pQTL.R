library(data.table)
library(susieR)
library(coloc)
library(SNPRelate)
library(gdsfmt)

## ---------- inputs ----------
fer_file   <- "Ferritin_AF0p005.mr_ready.tsv.gz"
pqtl_file  <- "prot-a-131.vcf.gz"

## APOE ±500 kb (GRCh38)
chr_apo   <- 19L
start_apo <- 44421094L
end_apo   <- 45421094L

## sample sizes
N_fer   <- 246000          
N_pqtl  <- 3300

## ---------- ferritin GWAS: slice & prep ----------
fer <- fread(fer_file)
stopifnot(all(c("SNP","chr","pos","beta","se","eaf") %in% names(fer)))
fer[, chr := as.integer(sub("^chr","", as.character(chr)))]
fer[, `:=`(beta = as.numeric(beta), se = as.numeric(se), eaf = as.numeric(eaf), pos = as.integer(pos))]
fer[, maf := pmin(eaf, 1 - eaf)]
gwas <- fer[chr == chr_apo & pos >= start_apo & pos <= end_apo &
              !is.na(beta) & !is.na(se) & !is.na(maf) & maf > 0 & maf < 1]
gwas <- unique(gwas, by = "SNP")
gwas <- gwas[nchar(effect_allele) == 1 & nchar(other_allele) == 1]

## ---------- APOE pQTL VCF: parse ES/SE/AF ----------
pqtl <- fread(pqtl_file, sep = "\t", header = TRUE, skip = "#CHROM", showProgress = FALSE)
setnames(pqtl, c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","PROT-a-131"),
         c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","STAT"))
fmt_names <- strsplit(pqtl$FORMAT[1], ":", fixed = TRUE)[[1]]
parts     <- tstrsplit(pqtl$STAT, ":", fixed = TRUE)
stopifnot(length(parts) == length(fmt_names))
pf <- as.data.table(setNames(parts, fmt_names))
need <- c("ES","SE","AF")
stopifnot(all(need %in% names(pf)))
pq <- cbind(pqtl[, .(CHROM, POS, ID, REF, ALT)], pf[, ..need])
pq[, `:=`(ES = as.numeric(ES), SE = as.numeric(SE), AF = as.numeric(AF))]
pq <- pq[!is.na(ES) & !is.na(SE) & !is.na(AF)]
pq <- pq[nchar(REF) == 1 & nchar(ALT) == 1]

## ---------- harmonise by rsID & align alleles ----------
aln <- merge(gwas[, .(SNP, effect_allele, other_allele, beta_gwas = beta, se_gwas = se, maf_gwas = maf)],
             pq[, .(SNP = ID, REF, ALT, beta_pqtl = ES, se_pqtl = SE, af_alt = AF)],
             by = "SNP")

pal <- aln[(effect_allele %in% c("A","T") & other_allele %in% c("A","T")) |
             (effect_allele %in% c("C","G") & other_allele %in% c("C","G"))]
aln <- rbind(
  aln[!(SNP %in% pal$SNP)],
  pal[abs(af_alt - 0.5) > 0.05]
)

aln <- aln[
  (effect_allele == ALT & other_allele == REF) |
    (effect_allele == REF & other_allele == ALT)
]

aln[, beta_pqtl_aligned :=
      fifelse(effect_allele == ALT,  beta_pqtl,
              fifelse(effect_allele == REF, -beta_pqtl, NA_real_))]
aln <- aln[!is.na(beta_pqtl_aligned)]

## ---------- restrict to SNPs present in 1000G & build LD ----------
rs_use <- unique(aln$SNP)

vcf_file <- "chr19.EUR.vcf.gz"   # same panel you used in Test 1
gds_file <- paste0(vcf_file, ".gds")

if (!file.exists(gds_file)) {
  tmp_gds <- tempfile(fileext = ".gds")
  snpgdsVCF2GDS(vcf.fn = vcf_file, out.fn = tmp_gds, method = "biallelic.only")
  file.rename(tmp_gds, gds_file)
}

genofile <- snpgdsOpen(gds_file)
rs_all   <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
id_all   <- read.gdsn(index.gdsn(genofile, "snp.id"))

idx <- match(rs_use, rs_all)
ok  <- !is.na(idx)
if (sum(ok) < 3L) {
  snpgdsClose(genofile)
  stop(sprintf("Only %d/%d requested rsIDs found in %s; need >=3 for a stable LD matrix.",
               sum(ok), length(ok), vcf_file))
}

# keep only the SNPs present in the panel, in the same order as rs_use[ok]
aln    <- aln[ok]
rs_use <- rs_use[ok]
snp_id <- id_all[idx[ok]]

ld_obj <- snpgdsLDMat(genofile, snp.id = snp_id, method = "corr", slide = -1)
R      <- ld_obj$LD
snpgdsClose(genofile)

genofile2 <- snpgdsOpen(gds_file)
fr <- snpgdsSNPRateFreq(genofile2, snp.id = snp_id)
snpgdsClose(genofile2)

panel_maf  <- pmin(fr$MinorFreq, 1 - fr$MinorFreq)
keep_mafcr <- is.finite(panel_maf) & (panel_maf > 0) & (fr$MissingRate <= 0.05)

cat(sprintf("Panel filter (MAF>0 & callrate>=0.95): kept %d/%d SNPs\n",
            sum(keep_mafcr), length(keep_mafcr)))

R      <- R[keep_mafcr, keep_mafcr, drop = FALSE]
aln    <- aln[keep_mafcr]
rs_use <- rs_use[keep_mafcr]

# Light numeric tidy (keeps R a correlation matrix)
R <- (R + t(R)) / 2
diag(R) <- 1
stopifnot(all(is.finite(R)), nrow(R) >= 3L)

## ---------- z-scores ----------
z1 <- aln$beta_gwas        / aln$se_gwas      # Ferritin
z2 <- aln$beta_pqtl_aligned / aln$se_pqtl     # APOE pQTL
names(z1) <- names(z2) <- rs_use
rownames(R) <- colnames(R) <- rs_use

## ---------- run SuSiE per trait ----------
Lmax <- 5  # keep same as your working run
fit1 <- susie_rss(z = z1, R = R, n = N_fer,  L = Lmax,
                  estimate_residual_variance = FALSE,
                  max_iter = 600, tol = 1e-3,
                  refine = FALSE, verbose = TRUE, track_fit = TRUE)

fit2 <- susie_rss(z = z2, R = R, n = N_pqtl, L = Lmax,
                  estimate_residual_variance = FALSE,
                  max_iter = 600, tol = 1e-3,
                  refine = FALSE, verbose = TRUE, track_fit = TRUE)

## ---------- coloc.susie ----------
cs <- coloc.susie(fit1, fit2)

print(cs$summary)
if (!is.null(cs$results)) {
  print(cs$results[order(-cs$results$SNP.PP.H4), ][1:min(10, .N)])
}


library(knitr)
library(data.table)

## -------- helpers --------
snp_ids <- colnames(R)  # susie_rss used names(z*) <- rs IDs; R has same order

comp_table <- function(fit, snp_ids, R) {
  a <- fit$alpha
  L <- nrow(a)
  lead_idx <- apply(a, 1, which.max)
  cs       <- susie_get_cs(fit, Xcorr = R)$cs
  data.table(
    component = seq_len(L),
    lead_snp  = snp_ids[lead_idx],
    pip_lead  = round(apply(a, 1, max), 6),
    cs_size   = vapply(seq_len(L), function(j) if (j <= length(cs) && length(cs[[j]])>0) length(cs[[j]]) else 0L, integer(1))
  )
}

posterior_overlap <- function(fit1, fit2, snp_ids, z1=NULL, z2=NULL, R=NULL) {
  a1 <- fit1$alpha; a2 <- fit2$alpha
  L1 <- nrow(a1);   L2 <- nrow(a2)
  PP <- a1 %*% t(a2)  # dot products (proxy for sharedness)
  # best SNP per pair by elementwise product
  best_idx <- vapply(seq_len(L1*L2), function(ix){
    j <- ((ix - 1L) %% L1) + 1L
    k <- ((ix - 1L) %/% L1) + 1L
    which.max(a1[j,] * a2[k,])
  }, integer(1))
  dt <- CJ(j = seq_len(L1), k = seq_len(L2))
  dt[, `:=`(PP_shared = as.vector(PP),
            best_snp  = snp_ids[best_idx])]
  if (!is.null(z1) && !is.null(z2)) {
    eg1 <- apply(a1, 1, which.max)
    eg2 <- apply(a2, 1, which.max)
    dt[, direction := ifelse(sign(z1[eg1])[j] * sign(z2[eg2])[k] > 0, "concordant", "discordant")]
  }
  setorder(dt, -PP_shared)
  dt
}

## -------- tables --------
tab_fer   <- comp_table(fit1, snp_ids, R)  # Ferritin SuSiE
tab_pqtl  <- comp_table(fit2, snp_ids, R)  # APOE pQTL SuSiE

cat("\n### Ferritin components (SuSiE)\n")
print(kable(tab_fer, format = "markdown"))

cat("\n### APOE pQTL components (SuSiE)\n")
print(kable(tab_pqtl, format = "markdown"))

cat("\n### coloc.susie summary (sorted by PP.H4)\n")
if (!is.null(cs$summary)) {
  sum_sorted <- as.data.table(cs$summary)[order(-PP.H4.abf)]
  print(kable(sum_sorted, format = "markdown"))
} else {
  cat("(No cs$summary returned)\n")
}

cat("\n### Per-SNP PP.H4\n")
if (!is.null(cs$results) && any(!is.na(cs$results$SNP.PP.H4))) {
  res_sorted <- cs$results[order(-cs$results$SNP.PP.H4), ][1:min(20, nrow(cs$results))]
  print(kable(res_sorted, format = "markdown"))
} else {
  cat("coloc.susie did not return per-SNP PP.H4; using SuSiE×SuSiE posterior-overlap proxy.\n")
  po <- posterior_overlap(fit1, fit2, snp_ids, z1, z2, R)
  print(kable(head(po, 20), format = "markdown",
              caption = "Top SuSiE×SuSiE component pairs by PP_shared (proxy)"))
}