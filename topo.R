library("TDA")
library("scatterplot3d")

X <- as.matrix(read.table('/home/garner1/Work/dataset/iFISH/for_silvano_NEW/chr1/iEG452_20171204_007_calc-002.NM-a594-1', sep=" "))

files <- list.files(path="/home/garner1/Work/dataset/iFISH/for_silvano_NEW/chr1", pattern=NULL, full.names=T, recursive=FALSE)
files <- list.files(path="/home/garner1/Work/dataset/iFISH/for_silvano_NEW/chr21", pattern=NULL, full.names=T, recursive=FALSE)

maxscale = 20
tseq <- seq(0,maxscale,length=1*maxscale)
for (file in files[1:5]){
  X <- as.matrix(read.table(file, sep=" "))
  filename<- paste(file, ".png", sep="")
  
  png(filename=filename)
  Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)

  # plot(Diag[["diagram"]], barcode=TRUE)
  # plot(Diag[["diagram"]], rotated=TRUE)
  Sil[ind, ] <- silhouette(Diag[["diagram"]],dimension = 0, p = 1, tseq)
  halfLife <- Sil[ind, ]
  plot(tseq, halfLife, type="l",main = file)
  dev.off()
  
  # Land1 <- landscape(Diag[["diagram"]],dimension = 0, KK = 1, tseq)
  # Land2 <- landscape(Diag[["diagram"]],dimension = 0, KK = 2, tseq)
  # Land3 <- landscape(Diag[["diagram"]],dimension = 0, KK = 3, tseq)
  # plot(tseq,Land1,type="l",main=file)
  # lines(tseq,Land2,col="green")
  # lines(tseq,Land3,col="red")
}
################
files <- list.files(path="/home/garner1/Work/dataset/iFISH/for_silvano_NEW/chr21", pattern=NULL, full.names=T, recursive=FALSE)
files <- list.files(path="/home/garner1/Work/dataset/iFISH/for_silvano_NEW/chr1_allele1", pattern=NULL, full.names=T, recursive=FALSE)

for (file in files[1:10]){
  X <- as.matrix(read.table(file, sep=" "))
  filename<- paste(file, ".png", sep="")
  # png(filename=filename)
  scatterplot3d(X, main = file) 
  # dev.off()
}
  
# s3d <- scatterplot3d( X1, color = "blue", xlim = c(min(X1[,1],X2[,1]),max(X1[,1],X2[,1])), ylim = c(min(X1[,2],X2[,2]),max(X1[,2],X2[,2])), 
#                       zlim = c(min(X1[,3],X2[,3]),max(X1[,3],X2[,3])) )
# s3d$points3d(X2, col="red")

################
files <- list.files(path="/home/garner1/Work/dataset/iFISH/data/iEG467/chr6", pattern="*.csv", full.names=T, recursive=FALSE)
maxscale = 30
tseq <- seq(0,maxscale,length=maxscale)
hdim = 0
lapply(files, function(x) {
  X <- as.matrix(read.table(x, sep=" "))
  if (nrow(X) >= 10) {
    # Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
    # tseq <- seq(0,maxscale=maxscale,length=maxscale)
    # Land1 <- landscape(Diag[["diagram"]],dimension = hdim, KK = 1, tseq)
    # Land2 <- landscape(Diag[["diagram"]],dimension = hdim, KK = 2, tseq)
    # Land3 <- landscape(Diag[["diagram"]],dimension = hdim, KK = 3, tseq)
    # plot(tseq,Land1,type="l",main = x)
    # lines(tseq,Land2,col="green")
    # lines(tseq,Land3,col="red")
    # plot(Diag[["diagram"]], barcode=TRUE, main = x)
    scatterplot3d(X,main = x)
  }
})
########################
maxscale = 100
tseq <- seq(0,maxscale,length=maxscale)
X <- as.matrix(read.table('/home/garner1/Work/dataset/iFISH/data/iEG467/chr2/iEG467_20171226_002_calc-001.NM_a594_1.csv', sep=" "))
Diag <- ripsDiag(X=X,maxdimension=1, maxscale=maxscale,library = "GUDHI",printProgress = TRUE)
hdim=0
Land1 <- landscape(Diag[["diagram"]],dimension = hdim, KK = 1, tseq)
Land2 <- landscape(Diag[["diagram"]],dimension = hdim, KK = 2, tseq)
Land3 <- landscape(Diag[["diagram"]],dimension = hdim, KK = 3, tseq)
plot(tseq,Land1,type="l",main = x)
lines(tseq,Land2,col="green")
lines(tseq,Land3,col="red")
Sil <- silhouette(Diag[["diagram"]],dimension = 0, p = 1, tseq)
plot(tseq, Sil, type="l",main = x)
# plot(Diag[["diagram"]], barcode=TRUE, main = x)
scatterplot3d(X)
