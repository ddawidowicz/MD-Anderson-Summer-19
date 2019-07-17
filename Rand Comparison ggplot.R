library(sigclust2)
library(mclust)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

z = 25 #iterations in inner loop
p = c(2, 10, 25, 50, 100, 200, 300, 500) #no. samples to test
df = as.data.frame(matrix(data = rep(p, each = z), nrow = length(p)*z, ncol = 7, byrow = FALSE))
colnames(df) = c ("p", "Ward.D2", "single", "average", "complete", "Mclust", "kmeans")
mean_df = as.data.frame(matrix(data = p, ncol = 7, nrow = length(p),  byrow = FALSE))
colnames(mean_df) = c ("p", "Ward.D2", "single", "average", "complete", "Mclust", "kmeans")

h = 1 #outer counting variable
k = 1 #inner counting variable

start.time = Sys.time()
for(i in 1:length(p)){
  for (j in 1:z){
    #Simulate data
    n1 <- 60; n2 <- 40; n3 <- 50; n <- n1 + n2 + n3
    data <- (matrix(rnorm(n*500), nrow=n, ncol=500))
    data[, 1] <- data[, 1] + c(rep(2, n1), rep(-2, n2), rep(0, n3))
    data[, 2] <- data[, 2] + c(rep(0, n1+n2), rep(sqrt(3)*3, n3))
    cluster_index = c(rep(1,n1), rep(2,n2), rep(3,n3)) #Real cluster index
    dataf = as.data.frame(data)
    
    #SigClust2
    shc_ward <- shc(data[,1:p[i]], metric="euclidean", linkage="ward.D2", alpha = 0.05)
    shc_single <- shc(data[,1:p[i]], metric="euclidean", linkage="single", alpha = 0.05)
    shc_average <- shc(data[,1:p[i]], metric="euclidean", linkage="average", alpha = 0.05)
    shc_complete <- shc(data[,1:p[i]], metric="euclidean", linkage="complete", alpha = 0.05)
  
    ward_index = cutree(shc_ward$hc_dat, k = 3)
    single_index = cutree(shc_single$hc_dat, k = 3)
    average_index = cutree(shc_average$hc_dat, k = 3)
    complete_index = cutree(shc_complete$hc_dat, k = 3)
    Mfit = Mclust(dataf[,1:p[i]])
    kfit = kmeans(data[,1:p[i]], centers = 3)
    
    df$Ward.D2[k] = adjustedRandIndex(ward_index, cluster_index)
    df$single[k] = adjustedRandIndex(single_index, cluster_index)
    df$average[k] = adjustedRandIndex(average_index, cluster_index)
    df$complete[k] = adjustedRandIndex(complete_index, cluster_index)
    df$Mclust[k] = adjustedRandIndex(Mfit$classification, cluster_index)
    df$kmeans[k] = adjustedRandIndex(kfit$cluster, cluster_index)
    
    k = k+1 #increment counting variable
  }
  mean_df$Ward.D2[i] = mean(df$Ward.D2[h:(k-1)])
  mean_df$single[i] = mean(df$single[h:(k-1)])
  mean_df$average[i] = mean(df$average[h:(k-1)])
  mean_df$complete[i] = mean(df$complete[h:(k-1)])
  mean_df$Mclust[i] = mean(df$Mclust[h:(k-1)])
  mean_df$kmeans[i] = mean(df$kmeans[h:(k-1)])
  
  h = h+z #increment outer counting variable
}
end.time = Sys.time()
end.time - start.time
#plot
mean_df_long = gather(mean_df, key = "method", value = "rand", -p)
df_long = gather(df, key = "method", value = "rand", -p)
g = ggplot(df_long, aes(x = p, y = rand)) + 
  geom_point(aes(col = method), size = 2) +
  geom_line(data = mean_df_long, aes(x = p, y = rand, col = method), size = 1) +
  labs(title = "Cluster Performance on Noisy Data", x = "# Fake Features", y = "Adjusted Rand Index") +
  scale_color_brewer(palette = "Set1")
plot(g)

