library(data.table)

thr_p_smr   <- 1e-8
thr_p_heidi <- 0.01
thr_nsnp    <- 3

f_blood <- fread("Ferritin and Blood.zip")
f_brain <- fread("Ferritin and Brain.zip")
f_liver <- fread("Ferritin and Liver.zip")

mark_pass <- function(dt, panel){
  setDT(dt)
  setnames(dt, old = intersect(names(dt), c("nsnp_HEIDI","Nsnp(HEIDI)","Nsnp")),
           new = "Nsnp", skip_absent = TRUE)
  dt[, pass := (p_SMR <= thr_p_smr) & (p_HEIDI >= thr_p_heidi) & (Nsnp >= thr_nsnp)]
  dt[, tissue := panel]
  dt[]
}

b <- mark_pass(f_blood, "Blood")
r <- mark_pass(f_brain, "Brain")
l <- mark_pass(f_liver, "Liver")

all_smr <- rbindlist(list(b,r,l), use.names = TRUE, fill = TRUE)
fwrite(all_smr, "ferritin_all_panels.annotated.txt", sep="\t")


panel_counts <- all_smr[, .(
  n_sig = sum(pass),
  n_total = .N
), by = tissue][, pct := round(100*n_sig/n_total, 2)]
panel_counts

passed <- all_smr[pass == TRUE, .(gene = index, tissue, p_SMR, p_HEIDI, topSNP)]

mat_pass <- dcast(passed[, .(gene, tissue, pass=TRUE)], gene ~ tissue, value.var = "pass", fun.aggregate = any)
fwrite(mat_pass, "ferritin_gene_by_tissue.tsv", sep="\t")

mat_pass[order(gene)][1:30]

lab <- copy(mat_pass)
for(col in c("Blood","Brain","Liver")) if(!col %in% names(lab)) lab[, (col) := FALSE]

lab[, class := fifelse(Brain & Blood, "Shared (Brain+Blood)",
                       fifelse(Brain & !Blood & !Liver, "Suggestive Brain-specific",
                               fifelse(Blood & !Brain & !Liver, "Suggestive Blood-specific",
                                       "Other/Unclear")))]
fwrite(lab, "ferritin_gene_by_tissue.classified.tsv", sep="\t")

targets <- c("SLC11A2","ORMDL1")
dossiers <- all_smr[index %in% targets & pass==TRUE,
                    .(index, tissue, topSNP, p_SMR, p_HEIDI, Nsnp)][order(index, tissue)]
fwrite(dossiers, "ferritin_locus_dossiers.tsv", sep="\t")

library(pheatmap)

dt <- copy(all_smr)

mat_p <- dcast(
  dt, index ~ tissue, value.var = "p_SMR",
  fun.aggregate = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
mat_h <- dcast(
  dt, index ~ tissue, value.var = "p_HEIDI",
  fun.aggregate = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))

P <- as.matrix(mat_p[, -1]); rownames(P) <- mat_p$index
L <- -log10(P); L[is.infinite(L)] <- NA; L[L > 30] <- 30  # cap for visibility

H <- as.matrix(mat_h[, -1]); rownames(H) <- mat_h$index
lab <- matrix("", nrow=nrow(L), ncol=ncol(L), dimnames=dimnames(L))
fail_cells <- (!is.na(P)) & (P <= thr_p_smr) & (!is.na(H)) & (H < thr_p_heidi)
lab[fail_cells] <- "HEIDI×"


finite_mask <- is.finite(L2)
keep_rows <- rowSums(finite_mask) >= 1
keep_cols <- colSums(finite_mask) >= 1

L2   <- L2[keep_rows, keep_cols, drop = FALSE]
lab2 <- lab2[keep_rows, keep_cols, drop = FALSE]

# now plot
pheatmap(
  L2,
  color = colorRampPalette(c("#F7FBFF","#6BAED6","#08306B"))(100),
  na_col = "grey90",
  display_numbers = lab2, number_color = "black", fontsize_number = 8,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, fontsize_row = 8, fontsize_col = 10,
  main = "-log10(p_SMR) with HEIDI failures marked"
)
