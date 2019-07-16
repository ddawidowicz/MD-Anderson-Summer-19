library(sigclust2)
library(mclust)
library(fpc)
library(factoextra)
library(gridExtra)
library(ggplot2)

#Simulate data
set.seed(1507)
n1 <- 60; n2 <- 40; n3 <- 50; n <- n1 + n2 + n3
p <- 500
data <- (matrix(rnorm(n*p), nrow=n, ncol=p))
data[, 1] <- data[, 1] + c(rep(2, n1), rep(-2, n2), rep(0, n3))
data[, 2] <- data[, 2] + c(rep(0, n1+n2), rep(sqrt(3)*3, n3))
cluster_index = c(rep(1,n1), rep(2,n2), rep(3,n3)) #Real cluster index
df = as.data.frame(data)

x = c(2, 10, 25, 50, 100, 200, 300, 500)
rand = matrix(nrow = 6, ncol = length(x))
colnames(rand) = as.character(x)
rownames(rand) = c ("Ward.D2", "single", "average", "complete", "Mclust", "kmeans")

for (i in 1:length(x)){
  #SigClust2
  shc_ward <- shc(data[,1:x[i]], metric="euclidean", linkage="ward.D2", alpha = 0.05)
  shc_single <- shc(data[,1:x[i]], metric="euclidean", linkage="single", alpha = 0.05)
  shc_avg <- shc(data[,1:x[i]], metric="euclidean", linkage="average", alpha = 0.05)
  shc_complete <- shc(data[,1:x[i]], metric="euclidean", linkage="complete", alpha = 0.05)
  
  ward_index = cutree(shc_ward$hc_dat, k = 3)
  rand[1,i] = adjustedRandIndex(ward_index, cluster_index)
  
  single_index = cutree(shc_single$hc_dat, k = 3)
  rand[2,i] = adjustedRandIndex(single_index, cluster_index)
  
  avg_index = cutree(shc_avg$hc_dat, k = 3)
  rand[3,i] = adjustedRandIndex(avg_index, cluster_index)
  
  complete_index = cutree(shc_complete$hc_dat, k = 3)
  rand[4,i] = adjustedRandIndex(complete_index, cluster_index)
  
  #Mclust
  Mfit = Mclust(df[,1:x[i]])
  rand[5,i] = adjustedRandIndex(Mfit$classification, cluster_index)
  
  #Kmeans
  kfit = kmeans(data[,1:x[i]], centers = 3)
  rand[6,i] = adjustedRandIndex(kfit$cluster, cluster_index)
}

#plot
plot(x, rand[1,], main = "Clustering Preformance on Noisy Data", 
     type="o",  col="red", pch="*", lwd=2, ylim = c(-0.1,1.5),
     xlab = "# fake features", ylab = "Adjusted Rand index")

points(x, rand[2,], col="darkorange", pch="*")
lines(x, rand[2,], col="darkorange", lwd = 2)

points(x, rand[3,], col="gold", pch="*")
lines(x, rand[3,], col="gold", lwd = 2)

points(x, rand[4,], col="forestgreen", pch="*")
lines(x, rand[4,], col="forestgreen", lwd = 2)

points(x, rand[5,], col="blue", pch="*")
lines(x, rand[5,], col="blue", lwd = 2)

points(x, rand[6,], col="darkorchid", pch="*")
lines(x, rand[6,], col="darkorchid", lwd = 2)

legend(350,1.5, legend=rownames(rand), col=c("red","darkorange","gold","forestgreen", "blue", "darkorchid"),
       lwd = 2, ncol=2)

