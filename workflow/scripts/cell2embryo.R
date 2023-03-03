#/scratch/cluster/kretzmer/SingleCell/KO_EED/SG7_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_EED SG7
#/scratch/cluster/kretzmer/SingleCell/KO_EED/SG8_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_EED SG8
#/scratch/cluster/kretzmer/SingleCell/KO_EED/SG6_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_EED SG6
#/scratch/cluster/kretzmer/SingleCell/KO_EED/agg.SG9.SG10_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_EED agg.SG9.SG10
#/scratch/cluster/kretzmer/SingleCell/KO_DNMT1/SG22_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_DNMT1 SG22
#/scratch/cluster/kretzmer/SingleCell/KO_DNMT1/SG11_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_DNMT1 SG11
#/scratch/cluster/kretzmer/SingleCell/KO_DNMT3B/SG13_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_DNMT3B SG13
#/scratch/cluster/kretzmer/SingleCell/KO_DNMT3B/SG15_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_DNMT3B SG15
#/scratch/cluster/kretzmer/SingleCell/KO_DNMT3B/agg.SG16.SG17_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_DNMT3B agg.SG16.S17
#/scratch/cluster/kretzmer/SingleCell/KO_TRIM28/SG21_SNPcount.tsv /scratch/cluster/kretzmer/SingleCell/KO_TRIM28 SG21
args <- commandArgs(TRUE)

# logfile
log <- file(snakemake@log[[1]], open="wt")
sink(log)

# libraries
require(RColorBrewer)
require(reshape2)
require(ggplot2)
require(scales)
require(sfsmisc)

set.seed(42)

# data bash
#SNPcountFile <- args[1]
#outDir <- args[2]
#ID <- args[3]
#NrOfEmbryos <- as.numeric(args[4])

# files
#pdf_file <- file.path(outDir, paste(ID, "_cell2embryo.pdf", sep=""))
#cluster_RData <- file.path(outDir, paste(ID, "_cell2embryo.RData", sep=""))

# data snakemake
SNPcountFile <- snakemake@input[[1]]
pdf_file <- snakemake@output[[1]] 
cluster_RData <- snakemake@output[[2]] 
snp_profile <- snakemake@output[[3]]
session_RData <- snakemake@output[[4]]
#NrOfEmbryos <- snakemake@params[["embryo_nr"]]
NrOfEmbryos <- NA
#ID <- snakemake@config[["ID"]]
ID <- snakemake@params[["identifier"]]

# variables
SNPcutoff <- 1000
SNPperc <- 0.2
corret_cluster_threshold <- 100
iterations <- 1000
considered_iterations <- 100
maxNA <- 2

# define colors
col_pals <- brewer.pal.info[brewer.pal.info$category=='qual',]
color <- unlist(mapply(brewer.pal, col_pals$maxcolors, rownames(col_pals)))

# open pdf for plotting statistics
pdf(pdf_file)

# read SNP count data
#data <- as.matrix(read.table(SNPcountFile, header=T, row.names=1))
#load(SNPcountFile)
data <- read.table(SNPcountFile, header=F, col.names=c("BC","ref","count"))
data <- acast(data, BC~ref, fill=0, fun.aggregate=sum)

rmX <- grep("chrX", colnames(data))
rmY <- grep("chrY", colnames(data))
if(length(c(rmX,rmY))>0){data <- data[,-c(rmX,rmY)]}

# filter for cells with >= SNPcutoff (1000) SNPs
SNPcounts <- data[which(rowSums(data)>=SNPcutoff), ]
# plot boxplot
ggplot(data=data.frame(count=rowSums(data)), aes(x="#SNPs", y=count)) + geom_boxplot() + theme_bw() + scale_y_continuous(trans=log_trans(), breaks=c(1,10,100,1000,10000))
ggplot(data=data.frame(SNPs=colMeans(SNPcounts), chr=gsub(".*G","G",colnames(SNPcounts),perl=T), g=gsub(".G.","",colnames(SNPcounts),perl=T)), aes(x=g, y=SNPs)) + geom_hline(yintercept=mean(rowMeans(SNPcounts))) + geom_bar(aes(fill=chr), position="dodge", width=0.5, color="black", stat='identity') + theme_classic()

# fill matrix with percent of CAST (G2) specific SNPs per chromosome
# remove cells with missing SNP information for a chromosome
SNPfraction <- SNPcounts[, seq(2, ncol(SNPcounts), 2)]/(SNPcounts[, seq(1, ncol(SNPcounts), 2)]+SNPcounts[,seq(2, ncol(SNPcounts), 2)])
SNPfraction <- SNPfraction[complete.cases(SNPfraction),]

# plot SNPfraction distribution per chromosome
ggplot(data=melt(SNPfraction), aes(x=value)) + geom_histogram(binwidth=0.01) + theme_bw() + facet_wrap(~Var2)


# function to calculate AIC based on kmeans clustering (based on http://sherrytowers.com/2013/10/24/k-means-clustering/)
kmeansAIC = function(fit){
  m <- ncol(fit$centers)
  n <- length(fit$cluster)
  k <- nrow(fit$centers)
  D <- fit$tot.withinss
  return(D + 2*m*k)
}

# for each k run kmeans on i random starts to find best model (k)
kmax <- 30               # the maximum number of clusters we will examine
imax <- 100              # rounds of kmean cluster per k
km_fit <- list()         # kmeans information for each k from best i/simulation (i.e. smallest totwss from each simulation)
all_km_totwss <- NULL    # totwss for each k from best i/simulation (i.e. smallest totwss from each simulation)
median_km_totwss <- NULL    # or each k median totwss (from each simulation)
median_km_AIC <- NULL
all_km_AIC <- NULL
for(k in 2:kmax){
  print(paste0("Run kmeans clustering for k=", k))
  all_i_totwss <- NULL       # all totwss from all simulations
  all_i_AIC <- NULL
  best_i_totwss <- Inf       # for each k smallest totwss (from each simulation)
  best_i_fit <- NULL         # kmeans information from best i (i.e. smallest totwss from each simulation)
  best_i_AIC <- NULL         # AIC for each k from best i (i.e. smallest totwss from each simulation)
  for(i in 1:imax){
    i_fit = kmeans(SNPfraction, k)
    i_totwss <- i_fit$tot.withinss
    i_AIC <- kmeansAIC(i_fit)
    all_i_totwss <- c(all_i_totwss, i_totwss)
    all_i_AIC <- c(all_i_AIC, i_AIC)
    if(i_totwss < best_i_totwss){
      best_i_totwss <- i_totwss
      best_i_fit <- i_fit
      best_i_AIC <- kmeansAIC(i_fit)
    }
  }
  km_fit[[k]] <- best_i_fit
  all_km_totwss <- c(all_km_totwss, best_i_totwss)
  median_km_totwss <- c(median_km_totwss, median(all_i_totwss))
  median_km_AIC <- c(median_km_AIC, median(all_i_AIC))
  all_km_AIC <- c(all_km_AIC, best_i_AIC)
}

# select the best model with at least 20% improvment in totwss
plot(seq(2, kmax), all_km_totwss, xlab="Number of clusters", ylab="totwss", pch=20, cex=2)
for(i in 1:kmax){
  # if i=kmax-1 then we have not found the solution of correct clusters yet and stop
  if(i==kmax) stop("No good model found!")
  if(min(all_km_totwss[(i+1):14])<(0.8*all_km_totwss[i])){
    next
  } else{
    nclus = i+1
    break
  }
}
points((i+1), all_km_totwss[i+1], col=2, pch=20, cex=2)
cat("Number of clusters (totwss): ", nclus, "\n")

# select the best model (taken from http://sherrytowers.com/2013/10/24/k-means-clustering/)
mult.fig(1,main="Simulated data with two clusters")
plot(seq(2,kmax),all_km_AIC,xlab="Number of clusters",ylab="AIC",pch=20,cex=2)
v <- -diff(all_km_AIC)
nv <- length(v)
fom <- v[1:(nv-1)]/v[2:nv]
nclus <- which.max(fom)+2
cat("Number of clusters (AIC): ", nclus, "\n")
points(nclus,all_km_AIC[nclus],col=2,pch=20,cex=2)

kclus <- km_fit[[nclus]]

if(!is.na(NrOfEmbryos)){
    nclus <- NrOfEmbryos
    kclus <- km_fit[[nclus]]
    cat("Predefined number of clusters: ", nclus, "\n")
}

# plot cell cluster using PCA representation
SNP_PCA <- prcomp(SNPfraction)
centers <- kclus$centers
plot(SNP_PCA$x[,c("PC1", "PC2")], col=color[kclus$cluster], main=paste0("PCA: SNP profiles of all cells in ", ID, " (Number of clusters: ", nclus, ")"))


# test for cluster stability, i.e. test if cells switch cluster when only a subset of SNPs is used for assignment
# remove those unstable (shaky) cells
# unstable means here that the cells might get closer to another cluster center (assigned to another embryo) using a subset of SNPs

# function to get squared euclidean distance from each cell to each cluster center and report index of min distance
subset_clusters <- function(x, centers) {
  distances <- sapply(seq_len(nrow(x)), function(i) apply(centers, 1, function(v) mean((x[i, ]-v)^2, na.rm=T)))
  return(max.col(-t(distances)))
}

# number of correct assignmets of cells (i.e. obtained using all SNPs)
correct_clustering <- rep(NA, nrow(SNPfraction))

# for each cell sample SNPs, calculate distance to all cluster centers, assign to "new" cluster, count clusterings to "correct" cluster
for(i in 1:nrow(SNPcounts)){
  if(i%%100 == 0){
      print(paste0("Distance calculation done with ", i, " (of ", nrow(SNPcounts), ") cells\n"))
  }

  # SNPperc percent of SNPs as subset
  subset_size <- floor(SNPperc*sum(SNPcounts[i,]))

  # simulate SNPs as name of chromosome repeated by number of SNPs per chromosome
  SNPs_per_chr <- rep(colnames(SNPcounts), SNPcounts[i,])

  # sample 100 times subset_size number of SNPs
  random_SNPcount <- t(replicate(iterations, table(sample(SNPs_per_chr)[1:subset_size])[colnames(SNPcounts)]))
  random_SNPcount[is.na(random_SNPcount)] <- 0
  colnames(random_SNPcount) <- colnames(SNPcounts)

  # fraction of CAST SNPs, preferentially use simulation with low number of NA
  random_SNPfraction <- random_SNPcount[, seq(2, ncol(random_SNPcount), 2)]/(random_SNPcount[, seq(1, ncol(random_SNPcount), 2)]+random_SNPcount[,seq(2, ncol(random_SNPcount), 2)])
  random_SNPfraction <- random_SNPfraction[order(rowSums(is.na(random_SNPfraction))),][1:considered_iterations,]

  # get "new" cluster and count number of correct assignments
  new_cluster <- subset_clusters(random_SNPfraction, centers)
  correct_clustering[i] <- sum(new_cluster  == kclus$cluster[i])
}

# keep cells with >=corret_cluster_threshold correct assignments
keep <- which(correct_clustering>=corret_cluster_threshold)

#and consider the plot with just these cells (which looks much better now)
plot(SNP_PCA$x[keep, c("PC1", "PC2")], col=color[kclus$cluster[keep]], main=paste0("PCA: SNP profiles of all stable cells in ", ID, " (Number of clusters: ", nclus, ")"))
legend("topright", as.character(data.frame(e=unique(kclus$cluster), c=color[unique(kclus$cluster)])$e), col=as.character(data.frame(e=unique(kclus$cluster), c=color[unique(kclus$cluster)])$c), pch=c(1))
#ggplot(as.data.frame(cbind(SNP_PCA$x[keep, ], kclus$cluster[keep])), aes(x=PC1, y=PC2)) + geom_point(aes(color=as.factor(V20))) + theme_classic()

#histogram of SNP number for cells which are NOT considered
#hist((1:nrow(SNPcounts))[-keep], breaks=100, xlab="Indices of removed cells", main=paste0("Histogram of indices of not considered cells (arranged according to SNP count)"), cex.main=0.9)
#ggplot(data.frame(index=(1:nrow(SNPcounts))[-keep]), aes(x=index)) + geom_histogram(color="black", fill="lightgrey", binwidth=100) + theme_bw()

# plot SNPprofile per embryo
t1 <- as.data.frame(SNPfraction)
t1$cell <- rownames(t1)
t2 <- data.frame(cell=names(kclus$cluster), embryo=kclus$cluster)
t3 <- merge(t1, t2)
t4 <- melt(t3, id=c("cell","embryo"))
#ggplot(data=t4, aes(x=value)) + geom_histogram(binwidth=0.01) + theme_bw() + facet_grid(embryo~variable)

# additional PCA plots to visualize location of instable cells
# new colors
colors_new <- color[kclus$cluster]
colors_new[intersect(which(rowSums(SNPcounts) < median(rowSums(SNPcounts)[keep])), seq(1, nrow(SNPcounts))[-keep])] <- "red"
colors_new[intersect(which(rowSums(SNPcounts) >= median(rowSums(SNPcounts)[keep])), seq(1, nrow(SNPcounts))[-keep])] <- "black"

plot(SNP_PCA$x[keep, c("PC1", "PC2")], col=colors_new[keep], main=paste0("PCA: SNP profiles of all cells in ", ID, " (Number of clusters: ", nclus, ")"))
points(SNP_PCA$x[-keep, c("PC1", "PC2")], col=colors_new[-keep])

# print barcode/cell to cluster/embryo assignment
cluster_assignment <- kclus$cluster[keep]
names(cluster_assignment) <- rownames(SNPcounts)[keep]
save(cluster_assignment, centers, file=cluster_RData)
dev.off()

#pdf(file.path(outDir, paste(ID, "_cell2embryo_SNPprofile.pdf", sep="")), width=28, height=14)
pdf(snp_profile, width=28, height=14)
ggplot(data=t4, aes(x=value)) + geom_histogram(binwidth=0.01) + theme_bw() + facet_grid(embryo~variable, scales="free_y")
dev.off()

#save.image(file.path(outDir, paste(ID, "_cell2embryo_session.RData", sep="")))
save.image(session_RData)
