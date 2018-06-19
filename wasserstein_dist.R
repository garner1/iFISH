library("TDA")
library("scatterplot3d")

files <- list.files(path="data/UserDots_IEG364_004", pattern="*.csv", full.names=T, recursive=FALSE)

lista = c()
for (file in files){
  X <- as.matrix(read.table(file, sep=","))
  if (nrow(X) >= 20){
    lista <- c(lista, file)
  }
}
print(length(lista))

distances_chr1_dim0 = c()
hdim = 0
for (file1 in lista[1:(length(lista))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  pos = match(file1,lista)
  if (pos < length(lista)){
    for (file2 in lista[(pos+1):length(lista)]){
      X2 <- as.matrix(read.table(file2, sep=","))
      Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
      count2 = nrow(X2)
      M <- max(count1,count2)
      distances_chr1_dim0 <- c(distances_chr1_dim0, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
    }
  }
}
distances_chr1_dim1 = c()
hdim = 1
for (file1 in lista[1:(length(lista))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  pos = match(file1,lista)
  if (pos < length(lista)){
    for (file2 in lista[(pos+1):length(lista)]){
      X2 <- as.matrix(read.table(file2, sep=","))
      Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
      count2 = nrow(X2)
      M <- max(count1,count2)
      distances_chr1_dim1 <- c(distances_chr1_dim1, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
    }
  }
}
##################
files <- list.files(path="data/UserDots_iEG408_003_d0_a594", pattern="*.csv", full.names=T, recursive=FALSE)

lista = c()
for (file in files){
  X <- as.matrix(read.table(file, sep=","))
  if (nrow(X) >= 20){
    lista <- c(lista, file)
  }
}
print(length(lista))

distances_chr20_dim0 = c()
hdim = 0
for (file1 in lista[1:(length(lista))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  pos = match(file1,lista)
  if (pos < length(lista)){
    for (file2 in lista[(pos+1):length(lista)]){
      X2 <- as.matrix(read.table(file2, sep=","))
      Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
      count2 = nrow(X2)
      M <- max(count1,count2)
      distances_chr20_dim0 <- c(distances_chr20_dim0, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
    }
  }
}
distances_chr20_dim1 = c()
hdim = 1
for (file1 in lista[1:(length(lista))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  pos = match(file1,lista)
  if (pos < length(lista)){
    for (file2 in lista[(pos+1):length(lista)]){
      X2 <- as.matrix(read.table(file2, sep=","))
      Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
      count2 = nrow(X2)
      M <- max(count1,count2)
      distances_chr20_dim1 <- c(distances_chr20_dim1, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
    }
  }
}
##################
files <- list.files(path="data/UserDots_iEG408_003_d0_cy5", pattern="*.csv", full.names=T, recursive=FALSE)

lista = c()
for (file in files){
  X <- as.matrix(read.table(file, sep=","))
  if (nrow(X) >= 20){
    lista <- c(lista, file)
  }
}
print(length(lista))

distances_chr2_dim0 = c()
hdim = 0
for (file1 in lista[1:(length(lista))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  pos = match(file1,lista)
  if (pos < length(lista)){
    for (file2 in lista[(pos+1):length(lista)]){
      X2 <- as.matrix(read.table(file2, sep=","))
      Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
      count2 = nrow(X2)
      M <- max(count1,count2)
      distances_chr2_dim0 <- c(distances_chr2_dim0, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
    }
  }
}
distances_chr2_dim1 = c()
hdim = 1
for (file1 in lista[1:(length(lista))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  pos = match(file1,lista)
  if (pos < length(lista)){
    for (file2 in lista[(pos+1):length(lista)]){
      X2 <- as.matrix(read.table(file2, sep=","))
      Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
      count2 = nrow(X2)
      M <- max(count1,count2)
      distances_chr2_dim1 <- c(distances_chr2_dim1, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
    }
  }
}
##################
# library(vioplot)
# vioplot(distances_chr1,distances_chr2,distances_chr20,names = c("chr1","chr2","chr20") )
# title("Violin plot for homological dim = 1")

#################
# h <- hist(distances,20,main = "chr2, homological dim = 1")
# print(h$breaks)

