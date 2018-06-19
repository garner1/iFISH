library("TDA")
library("scatterplot3d")

# SET ylim, hdim and title of the output file to switch between hdim 0 and 1

hdim = 1
files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr1", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr1_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr10", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr10_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr11", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr11_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr13", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr13_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr15", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr15_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr17", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr17_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr19", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr19_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr2", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr2_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr21", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr21_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr22", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr22_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr5", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr5_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr6", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr6_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr7", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr7_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr8", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr8_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chr9", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chr9_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################

files <- list.files(path="/home/garner1/Work/dataset/iFISH/G1_cells_only/numbClusters_1/chrX", pattern=NULL, full.names=T, recursive=FALSE)
ind = 0
l = c()
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  l[ind] <- max(dist(X,method = "euclidean")) #find max distance in cloud
}
maxscale = max(l)/2 # set maxscale to half of max distance
tseq <- seq(0,maxscale,length=1*maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
p = 1
ind = 0
for (x in files) {
  ind = ind+1
  X <- as.matrix(read.table(x, sep=" "))
  # cat(ind, " of ", length(files)," and max dist equal to ", max(dist(X,method = "euclidean")),"\n")
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}
png("chrX_d1.png")
bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = FALSE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
dev.off()
################################################################################
