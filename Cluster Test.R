library(sigclust2)
library(mclust)
library(fpc)
library(factoextra)
library(gridExtra)

#Simulate data
set.seed(1508)
n1 <- 60; n2 <- 40; n3 <- 50; n <- n1 + n2 + n3
p <- 100
data <- (matrix(rnorm(n*p), nrow=n, ncol=p))
data[, 1] <- data[, 1] + c(rep(2, n1), rep(-2, n2), rep(0, n3))
data[, 2] <- data[, 2] + c(rep(0, n1+n2), rep(sqrt(3)*3, n3))
data <- data[, 1:2]
df = as.data.frame(data)

cluster_index = c(rep(1,n1), rep(2,n2), rep(3,n3)) #Real cluster index

#SigClust2
shc_ward <- shc(data, metric="euclidean", linkage="ward.D2", alpha = 0.05)
shc_single <- shc(data, metric="euclidean", linkage="single", alpha = 0.05)
shc_avg <- shc(data, metric="euclidean", linkage="average", alpha = 0.05)
shc_complete <- shc(data, metric="euclidean", linkage="complete", alpha = 0.05)

ward_index = cutree(shc_ward$hc_dat, k = 3)
a = adjustedRandIndex(ward_index, cluster_index)

single_index = cutree(shc_single$hc_dat, k = 3)
b = adjustedRandIndex(single_index, cluster_index)

avg_index = cutree(shc_avg$hc_dat, k = 3)
c = adjustedRandIndex(avg_index, cluster_index)

complete_index = cutree(shc_complete$hc_dat, k = 3)
d = adjustedRandIndex(complete_index, cluster_index)

#plot
ward_plot = plot(shc_ward)
single_plot = plot(shc_single)
avg_plot = plot(shc_avg)
complete_plot = plot(shc_complete, main = "complete")
grid.arrange(ward_plot, single_plot, avg_plot, complete_plot, nrow = 2, top = "SHC Methods",
             bottom = "[1,1] = Ward.D2   [1,2] = single    [2,1] = average    [2,2] = complete")

#Mclust
Mfit = Mclust(df)
e = adjustedRandIndex(Mfit$classification, cluster_index)

#Fviz
wss = fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)+
  labs(subtitle = "Elbow method")

sil = fviz_nbclust(df, kmeans, method = "silhouette") +
  labs(subtitle = "Silhouette method")

gap = fviz_nbclust(df, kmeans, method = "gap_stat") +
  labs(subtitle = "gap method")

grid.arrange(wss, sil, gap, nrow = 2)

#Kmeans
kfit = kmeans(data[,1:2], centers = 3)
f = adjustedRandIndex(kfit$cluster, cluster_index)

#Compare Cluster Algorithms
results = matrix(c(a, b, c, d, e, f), ncol = 1, byrow = TRUE)
colnames(results) = "C. Rand"
rownames(results) = c("ward.D2", "single", "average", "complete", "Mclust", "Kmeans")
results = as.table(results)
results