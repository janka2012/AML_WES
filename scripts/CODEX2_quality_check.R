# based on https://github.com/yuchaojiang/CODEX2

library(CODEX2)
library(BSgenome.Hsapiens.UCSC.hg38)

# usage
# Rscript --vanilla CODEX2_quality_check.R ALA ALA_new ALA_new_new

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Preprocessing
PROJECT_DIR=file.path("/storage","brno3-cerit", "home", "xputer00", "000020-Shares", "999990-rcxag", "Projects", "Ostrava_JP", "Ostrava_data")
print(PROJECT_DIR)

TARGETS_DIR=file.path(PROJECT_DIR, "targets")
bedFile <- file.path(TARGETS_DIR, "SureSelectV5_targets_hg38.srt.no_alt.no_ovl.bed")
print(bedFile)

bamFiles_all = c()
bamdirs_all = c()
for (i in 1:length(args)) { 
	DATADIR=file.path(PROJECT_DIR, "analysis", args[i], "realign_recalibrate", "good_samples")
	bamFiles <- list.files(DATADIR, pattern = '*.bam$')
	bamdirs <- file.path(DATADIR, bamFiles)
        bamFiles_all <- c(bamFiles_all, bamFiles)
        bamdirs_all <- c(bamdirs_all, bamdirs)
}


sample_names <- sapply(bamFiles_all, function(x) strsplit(x, "_")[[1]][1], USE.NAMES=FALSE)
DIAGNOSE = paste(args, collapse="_")
bambedObj <- getbambed(bamdir = bamdirs_all, bedFile = bedFile, sampname = sample_names, projectname = paste0("CODEX2_", DIAGNOSE))

bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname
genome(ref) <- "hg38"

# Getting GC content and mappability
# obrain gc content and mappability for each exon/target
gc <- getgc(ref, genome = BSgenome.Hsapiens.UCSC.hg38)
mapp <- getmapp(ref, genome = BSgenome.Hsapiens.UCSC.hg38)
values(ref) <- cbind(values(ref), DataFrame(gc, mapp))

# Getting raw read depth
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y

coverage_file = file.path(PROJECT_DIR, "analysis", DIAGNOSE, "clonality", "copynumbers", "CODEX2", paste0(DIAGNOSE, '_coverage.csv'))

write.csv(Y, file = coverage_file, quote = FALSE)

# Quality control
qcObj <- qc(Y, sampname, ref, cov_thresh = c(10, 4000), length_thresh = c(20, 20000),
            mapp_thresh = 0.9, gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc

qc_file <- file.path(PROJECT_DIR, "analysis", DIAGNOSE, "clonality", "copynumbers", "CODEX2", paste0(DIAGNOSE, '_qcmat.csv'))
write.table(qcmat, file = qc_file, sep = '\t', quote = FALSE, row.names = FALSE)


# Estimate library size factor based on genome-wide read depth after QC
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){prod(x)^(1/length(x))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)

norm_samples_index <- grep("PB$", colnames(Y_qc))

output_file = file.path(PROJECT_DIR, "analysis", DIAGNOSE, "clonality", "copynumbers", "CODEX2", paste0(DIAGNOSE, "_workspace.RData"))

save.image(file = output_file)
