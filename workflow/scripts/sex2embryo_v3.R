args <- commandArgs(TRUE)

# logfile
log <- file(snakemake@log[[1]], open="wt")
sink(log)

# libraries
library(Matrix)
set.seed(42)

# data
#cellranger_dir <- args[1]
#matrix_dir = args[1]
#cell2embryo_output <- args[2]
#outDir <- args[3]
#ID <- args[4]

cell2embryo_output <- snakemake@input[[1]]
#matrix_dir <- snakemake@config[["matrix_dir"]]
matrix_dir <- snakemake@params[["matrix_dir"]]
#outDir <- snakemake@config[["outDir"]]
#ID <- snakemake@config[["ID"]]

# files
#out_RData <- file.path(outDir, paste(ID, "_sex2embryo.RData", sep=""))
out_RData <- snakemake@output[[1]]

# variables
genome <- "mm10"

XY <- c("ENSMUSG00000096768", "ENSMUSG00000069045", "ENSMUSG00000069049")
XX <- "ENSMUSG00000086503"

# load cellranger output and extract RNA expression information
#cellranger <- load_cellranger_matrix(cellranger_dir, genome=genome)
#expr <- as.matrix(exprs(cellranger))
#colnames(expr) <- gsub("-.", "", colnames(expr))
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1

# reduce to XX and XY determining genes
expr <- mat[c(XX,XY), ]

# load BC (cell) to cluster (embryo) matrix
load(cell2embryo_output)

if(!(length(grep("-1", colnames(expr)))>0 & length(grep("-1", names(cluster_assignment)))>0)){
    names(cluster_assignment) <- paste0(names(cluster_assignment), "-1")
}

# calculate fraction of cells per cluster expressing XX or XY genes
sex_classification <- data.frame(cluster=integer(), XX=integer(), XY=integer())
for(i in 1:max(cluster_assignment)){
    XXexpr <- mean(expr[XX, which(colnames(expr) %in% names(cluster_assignment[cluster_assignment==i]))]>0)
    XYexpr <- mean(expr[XY, which(colnames(expr) %in% names(cluster_assignment[cluster_assignment==i]))]>0)
    sex_classification[i,] <- c(i, XXexpr, XYexpr)
}
save(centers, cluster_assignment, sex_classification, XX, XY, file=out_RData)
