library(sigclust2)
library(mclust)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

#Simulate data
set.seed(1507)
n1 <- 60; n2 <- 40; n3 <- 50; n <- n1 + n2 + n3
p <- 500
data <- (matrix(rnorm(n*p), nrow=n, ncol=p))
data[, 1] <- data[, 1] + c(rep(2, n1), rep(-2, n2), rep(0, n3))
data[, 2] <- data[, 2] + c(rep(0, n1+n2), rep(sqrt(3)*3, n3))
cluster_index = c(rep(1,n1), rep(2,n2), rep(3,n3)) #Real cluster index
dataf = as.data.frame(data)

p = c(2, 10, 25, 50, 100, 200, 300, 500)
df = as.data.frame(matrix(data = p, ncol = 7, nrow = length(p),  byrow = FALSE))
colnames(df) = c ("p", "Ward.D2", "single", "average", "complete", "Mclust", "kmeans")

for (i in 1:length(x)){
  #SigClust2
  shc_ward <- shc(data[,1:x[i]], metric="euclidean", linkage="ward.D2", alpha = 0.05)
  shc_single <- shc(data[,1:x[i]], metric="euclidean", linkage="single", alpha = 0.05)
  shc_avg <- shc(data[,1:x[i]], metric="euclidean", linkage="average", alpha = 0.05)
  shc_complete <- shc(data[,1:x[i]], metric="euclidean", linkage="complete", alpha = 0.05)
  
  ward_index = cutree(shc_ward$hc_dat, k = 3)
  df$Ward.D2[i] = adjustedRandIndex(ward_index, cluster_index)
  
  single_index = cutree(shc_single$hc_dat, k = 3)
  df$single[i] = adjustedRandIndex(single_index, cluster_index)
  
  avg_index = cutree(shc_avg$hc_dat, k = 3)
  df$average[i] = adjustedRandIndex(avg_index, cluster_index)
  
  complete_index = cutree(shc_complete$hc_dat, k = 3)
  df$complete[i] = adjustedRandIndex(complete_index, cluster_index)
  
  #Mclust
  Mfit = Mclust(dataf[,1:x[i]])
  df$Mclust[i] = adjustedRandIndex(Mfit$classification, cluster_index)
  
  #Kmeans
  kfit = kmeans(data[,1:x[i]], centers = 3)
  df$kmeans[i] = adjustedRandIndex(kfit$cluster, cluster_index)
}

#plot
df_long = gather(df, key = "method", value = "rand", -p)
g = ggplot(df_long, aes(x = p, y = rand)) + 
  geom_line(aes(col = method), size = 1) +
  geom_point(aes(col = method)) +
  labs(title = "Cluster Performance on Noisy Data", x = "# Samples", y = "Adjusted Rand Index") +
  scale_color_brewer(palette = "Set1")
plot(g)

