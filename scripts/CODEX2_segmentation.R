library(doParallel)
library(CODEX2)
library(BSgenome.Hsapiens.UCSC.hg38)

# usage
# Rscript --vanilla CODEX2_segmentation.R ALA ALA_new ALA_new_new


args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

DIAGNOSE = paste(args, collapse="_")

PROJECT_DIR=file.path("/storage","brno3-cerit", "home", "xputer00", "000020-Shares", "999990-rcxag", "Projects", "Ostrava_JP", "Ostrava_data")

output_file = file.path(PROJECT_DIR, "analysis", DIAGNOSE, "clonality", "copynumbers", "CODEX2", paste0(DIAGNOSE, "_workspace.RData"))
chromosomes <- c(paste0('chr',1:22))

# Running CODEX2 with negative control samples
load(output_file)

print("Rdata loaded.")

# Below are pre-computed demo dataset, stored as part of the CODEX2 R-package.

cl <- makeForkCluster(4)
registerDoParallel(cl)

results <- foreach(chr=chromosomes, .packages='CODEX2', .inorder=TRUE) %dopar% {
  chr.index <- which(seqnames(ref_qc)==chr)
  normObj <- normalize_codex2_ns(Y_qc = Y_qc[chr.index,],
                                 gc_qc = gc_qc[chr.index],
                                 K = 1:5, norm_index = norm_samples_index, N = N)

  Yhat.ns <- normObj$Yhat; fGC.hat.ns <- normObj$fGC.hat;
  beta.hat.ns <- normObj$beta.hat; g.hat.ns <- normObj$g.hat; h.hat.ns <- normObj$h.hat
  AIC.ns <- normObj$AIC; BIC.ns <- normObj$BIC; RSS.ns <- normObj$RSS

  choiceofK_pdf = file.path(PROJECT_DIR,  "analysis", DIAGNOSE, "clonality", "copynumbers", "CODEX2", paste0("codex2_ns_choiceofK", chr, ".pdf"))
  choiceofK(AIC.ns, BIC.ns, RSS.ns, K = 1:5 , filename = choiceofK_pdf)

  # running segmentation
  finalcall.CBS <- segmentCBS(Y_qc[chr.index,],  # recommended
                              Yhat.ns, optK = which.max(BIC.ns),
                              K = 1:5,
                              sampname_qc = colnames(Y_qc),
                              ref_qc = ranges(ref_qc)[chr.index],
                              chr = chr, lmax = 400, mode = "fraction")

  # Post segmentation filtering
  filter1 <- finalcall.CBS$length_kb <= 200 # less than 200 kb
  filter2 <- finalcall.CBS$length_kb/(finalcall.CBS$ed_exon-finalcall.CBS$st_exon+1) < 50 # exon length
  finalcall.CBS.filter <- finalcall.CBS[filter1 & filter2, ]

  filter3 <- finalcall.CBS.filter$lratio>40 # likelihood ratio ????
  filter4 <- (finalcall.CBS.filter$ed_exon-finalcall.CBS.filter$st_exon) > 1 # number of exons
  final <- finalcall.CBS.filter[filter3|filter4,]
  final
}

# ---------------- END PARALELIZE ------------------------------------------------


stopCluster(cl)

copynumbers_file <- file.path(PROJECT_DIR, "analysis", DIAGNOSE, "clonality", "copynumbers", "CODEX2", paste0(DIAGNOSE, '_copynumbers.csv')) 
print(copynumbers_file)
print("Merging data frames....") 
results_df <- do.call("rbind", results)
print("Saving workspace....")
save.image(file= file.path(PROJECT_DIR, "analysis", DIAGNOSE, "clonality","copynumbers", "CODEX2", paste0(DIAGNOSE, "_copynumbers.Rdata"))) 
print("Writing output...")
write.table(results_df, file = copynumbers_file, sep = '\t', quote = FALSE, row.names = FALSE)
