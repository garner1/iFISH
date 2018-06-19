library("TDA")
library("scatterplot3d")

files_chr1 <- list.files(path="data/UserDots_IEG364_004", pattern="*.csv", full.names=T, recursive=FALSE)
files_chr2 <- list.files(path="data/UserDots_iEG408_003_d0_cy5", pattern="*.csv", full.names=T, recursive=FALSE)
files_chr20 <- list.files(path="data/UserDots_iEG408_003_d0_a594", pattern="*.csv", full.names=T, recursive=FALSE)

lista_1 = c()
for (file in files_chr1){
  X <- as.matrix(read.table(file, sep=","))
  if (nrow(X) >= 20){
    lista_1 <- c(lista_1, file)
  }
}

lista_2 = c()
for (file in files_chr2){
  X <- as.matrix(read.table(file, sep=","))
  if (nrow(X) >= 20){
    lista_2 <- c(lista_2, file)
  }
}

lista_20 = c()
for (file in files_chr20){
  X <- as.matrix(read.table(file, sep=","))
  if (nrow(X) >= 20){
    lista_20 <- c(lista_20, file)
  }
}

distances_chr1_chr2_dim0 = c()
hdim = 0
for (file1 in lista_1[1:(length(lista_1))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  for (file2 in lista_2[1:length(lista_2)]){
    X2 <- as.matrix(read.table(file2, sep=","))
    Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
    count2 = nrow(X2)
    M <- max(count1,count2)
    distances_chr1_chr2_dim0 <- c(distances_chr1_chr2_dim0, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
  }
}
################
distances_chr1_chr2_dim1 = c()
hdim = 1
for (file1 in lista_1[1:(length(lista_1))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  for (file2 in lista_2[1:length(lista_2)]){
    X2 <- as.matrix(read.table(file2, sep=","))
    Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
    count2 = nrow(X2)
    M <- max(count1,count2)
    distances_chr1_chr2_dim1 <- c(distances_chr1_chr2_dim1, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
  }
}
##################
distances_chr1_chr20_dim0 = c()
hdim = 0
for (file1 in lista_1[1:(length(lista_1))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  for (file2 in lista_20[1:length(lista_20)]){
    X2 <- as.matrix(read.table(file2, sep=","))
    Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
    count2 = nrow(X2)
    M <- max(count1,count2)
    distances_chr1_chr20_dim0 <- c(distances_chr1_chr20_dim0, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
  }
}

##################
distances_chr1_chr20_dim1 = c()
hdim = 1
for (file1 in lista_1[1:(length(lista_1))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  for (file2 in lista_20[1:length(lista_20)]){
    X2 <- as.matrix(read.table(file2, sep=","))
    Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
    count2 = nrow(X2)
    M <- max(count1,count2)
    distances_chr1_chr20_dim1 <- c(distances_chr1_chr20_dim1, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
  }
}
##################
distances_chr2_chr20_dim0 = c()
hdim = 0
for (file1 in lista_2[1:(length(lista_2))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  for (file2 in lista_20[1:length(lista_20)]){
    X2 <- as.matrix(read.table(file2, sep=","))
    Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
    count2 = nrow(X2)
    M <- max(count1,count2)
    distances_chr2_chr20_dim0 <- c(distances_chr2_chr20_dim0, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
  }
}

##################
distances_chr2_chr20_dim1 = c()
hdim = 1
for (file1 in lista_2[1:(length(lista_2))]){
  X1 <- as.matrix(read.table(file1, sep=","))
  Diag1 <- ripsDiag(X=X1,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
  count1 = nrow(X1)
  for (file2 in lista_20[1:length(lista_20)]){
    X2 <- as.matrix(read.table(file2, sep=","))
    Diag2 <- ripsDiag(X=X2,maxdimension=2, maxscale=5000,library = "GUDHI",printProgress = FALSE)
    count2 = nrow(X2)
    M <- max(count1,count2)
    distances_chr2_chr20_dim1 <- c(distances_chr2_chr20_dim1, wasserstein(Diag1[["diagram"]],Diag2[["diagram"]], p = 2, dimension = hdim)/M)
  }
}

##################
library(vioplot)

vioplot(distances_chr1_dim0,distances_chr2_dim0,distances_chr20_dim0,distances_chr1_chr2_dim0,distances_chr1_chr20_dim0,distances_chr2_chr20_dim0,
        names = c("chr1","chr2","chr20","chr1-2","chr1-20","chr2-20") )
title("intra- and inter-chromosomal violin plot for homological dim = 0")

vioplot(distances_chr1_dim1,distances_chr2_dim1,distances_chr20_dim1,distances_chr1_chr2_dim1,distances_chr1_chr20_dim1,distances_chr2_chr20_dim1,
        names = c("chr1","chr2","chr20","chr1-2","chr1-20","chr2-20") )
title("intra- and inter-chromosomal violin plot for homological dim = 1")

#################
# h <- hist(distances_chr1_chr2_dim0,20,main = "chr1 and chr2, homological dim = 0")
# h <- hist(distances_chr1_chr2_dim1,20,main = "chr1 and chr2, homological dim = 1")
# h <- hist(distances_chr1_chr20_dim0,20,main = "chr1 and chr20, homological dim = 0")
# h <- hist(distances_chr1_chr20_dim1,20,main = "chr1 and chr20, homological dim = 1")
# h <- hist(distances_chr2_chr20_dim0,20,main = "chr2 and chr20, homological dim = 0")
# h <- hist(distances_chr2_chr20_dim1,20,main = "chr2 and chr20, homological dim = 1")
# 
# print(h$breaks)
 
