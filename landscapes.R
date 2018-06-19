library("TDA")
library("scatterplot3d")

files <- list.files(path="data/UserDots_IEG364_004", pattern="*.csv", full.names=T, recursive=FALSE)
files <- list.files(path="data/UserDots_iEG408_003_d0_a594", pattern="*.csv", full.names=T, recursive=FALSE)
files <- list.files(path="data/UserDots_iEG408_003_d0_cy5", pattern="*.csv", full.names=T, recursive=FALSE)

maxscale = 5000
tseq <- seq(0,maxscale,length=maxscale)
Land <- matrix(0, nrow=length(files), ncol = length(tseq))
hdim = 1
ind = 0
KK = 1
for (x in files) {
  ind = ind+1
  cat(ind, " of ", length(files),"\n")
  X <- as.matrix(read.table(x, sep=","))
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=5000,library = "GUDHI",printProgress = TRUE)
  Land[ind, ] <- landscape(Diag[["diagram"]],dimension = hdim, KK = KK, tseq)
}

bootLand <- multipBootstrap(Land, B = 100, alpha = 0.05, parallel = TRUE, printProgress = TRUE)
plot(tseq, bootLand[["mean"]], main="Mean Landscape")
polygon(c(tseq,rev(tseq)),c(bootLand[["band"]][, 1],rev(bootLand[["band"]][, 2])),col = "pink")
lines(tseq,bootLand[["mean"]],lwd = 2, col = 2)
#####################################################
files <- list.files(path="data/UserDots_IEG364_004", pattern="*.csv", full.names=T, recursive=FALSE)
files <- list.files(path="data/UserDots_iEG408_003_d0_a594", pattern="*.csv", full.names=T, recursive=FALSE)
files <- list.files(path="data/UserDots_iEG408_003_d0_cy5", pattern="*.csv", full.names=T, recursive=FALSE)

maxscale = 5000
tseq <- seq(0,maxscale,length=maxscale)
Sil <- matrix(0, nrow=length(files), ncol = length(tseq))
hdim = 1
ind = 0
p = 1
for (x in files) {
  ind = ind+1
  cat(ind, " of ", length(files),"\n")
  X <- as.matrix(read.table(x, sep=","))
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=5000,library = "GUDHI",printProgress = TRUE)
  Sil[ind, ] <-silhouette(Diag[["diagram"]],dimension = hdim, p = p, tseq)
}

bootLand <- multipBootstrap(Sil, B = 100, alpha = 0.05, parallel = TRUE, printProgress = TRUE)
plot(tseq, bootLand[["mean"]], main="Mean Silhouette")
polygon(c(tseq,rev(tseq)),c(bootLand[["band"]][, 1],rev(bootLand[["band"]][, 2])),col = "pink")
lines(tseq,bootLand[["mean"]],lwd = 2, col = 2)
