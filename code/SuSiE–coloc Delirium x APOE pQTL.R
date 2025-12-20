library(data.table)
library(susieR)
library(coloc)
library(SNPRelate)
library(gdsfmt)

## ---------- inputs ----------
del_file   <- "Delirium_AF0p005.mr_ready.tsv.gz"
pqtl_file  <- "prot-a-131.vcf.gz"

## APOE ±500 kb (GRCh38)
chr_apo   <- 19L
start_apo <- 44421094L
end_apo   <- 45421094L

## sample sizes
N_cases   <- 8461
N_ctrls   <- 449979
N_del_eff <- 4 / (1/N_cases + 1/N_ctrls)
N_pqtl  <- 3300

## ---------- delirium GWAS: slice & prep ----------
del <- fread(del_file)
stopifnot(all(c("SNP","chr","pos","beta","se","eaf") %in% names(del)))
del[, chr := as.integer(sub("^chr","", as.character(chr)))]
del[, `:=`(beta = as.numeric(beta), se = as.numeric(se), eaf = as.numeric(eaf), pos = as.integer(pos))]
del[, maf := pmin(eaf, 1 - eaf)]
gwas <- del[chr == chr_apo & pos >= start_apo & pos <= end_apo &
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

vcf_file <- "chr19.EUR.vcf.gz"
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

# Give all objects consistent SNP names
rownames(R) <- colnames(R) <- rs_use
names(z1)   <- rs_use
names(z2)   <- rs_use

# Light numeric tidy (keeps R a correlation matrix)
R <- (R + t(R)) / 2
diag(R) <- 1

# Quick sanity check
stopifnot(all(is.finite(R)), nrow(R) >= 3L)

## ---------- z-scores ----------
z1 <- aln$beta_gwas / aln$se_gwas
z2 <- aln$beta_pqtl_aligned / aln$se_pqtl

## ---------- run SuSiE per trait ----------
Lmax <- 5  # was 10
fit1 <- susie_rss(z = z1, R = R, n = N_del_eff, L = Lmax,
                  estimate_residual_variance = FALSE,
                  max_iter = 600, tol = 1e-3,
                  refine = FALSE, verbose = TRUE, track_fit = TRUE)

fit2 <- susie_rss(z = z2, R = R, n = N_pqtl,    L = Lmax,
                  estimate_residual_variance = FALSE,
                  max_iter = 600, tol = 1e-3,
                  refine = FALSE, verbose = TRUE, track_fit = TRUE)
## ---------- coloc.susie ----------
cs <- coloc.susie(fit1, fit2)

print(cs$summary)
if (!is.null(cs$results)) {
  print(cs$results[order(-cs$results$SNP.PP.H4), ][1:min(10, .N)])
}
