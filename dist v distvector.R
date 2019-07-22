library(sigclust2)
library(mclust)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

distvector = function(X){
  d = as.matrix(dist(X))
  n = dim(d)[1]
  E = matrix(nrow = n, ncol = n)
  
  for (i in 1:n){
    for(j in 1:n){
      temp = 0
      for (t in 1:n){
        if(t != i & t != j){
          temp = temp + (d[i,t] - d[j,t])^2
        }
        E[i,j] = sqrt(temp)
      }
    }
  }
  E = as.dist(E)
  return(E)
}

z = 10 #iterations in inner loop
p = c(2, 10, 50, 100, 250, 500) #no. samples to test
df = as.data.frame(matrix(data = rep(p, each = z), nrow = length(p)*z, ncol = 3, byrow = FALSE))
colnames(df) = c("p", "dist", "distvector")
mean_df = as.data.frame(matrix(data = p, ncol = 3, nrow = length(p),  byrow = FALSE))
colnames(mean_df) = c("p", "dist", "distvector")

dist_rand = vector(length = z)
distvector_rand =vector(length = z)

h = 1 #outer counting variable
k = 1 #inner counting variable

for(i in 1:length(p)){
  for (j in 1:z){
    #Simulate data
    n1 <- 60; n2 <- 40; n3 <- 50; n <- n1 + n2 + n3
    data <- (matrix(rnorm(n*500), nrow=n, ncol=500))
    data[, 1] <- data[, 1] + c(rep(2, n1), rep(-2, n2), rep(0, n3))
    data[, 2] <- data[, 2] + c(rep(0, n1+n2), rep(sqrt(3)*3, n3))
    cluster_index = c(rep(1,n1), rep(2,n2), rep(3,n3)) #Real cluster index
    dataf = as.data.frame(data)
    
    #cluster
    dist_fit =  hclust(dist(data[,1:p[i]]), method = "ward.D2")
    distvector_fit = hclust(distvector(data[,1:p[i]]), method = "ward.D2")
    
    #index
    dist_index = cutree(dist_fit, 3)
    distvector_index = cutree(distvector_fit, 3)
    
    df$dist[k] = adjustedRandIndex(dist_index, cluster_index)
    df$distvector[k] = adjustedRandIndex(distvector_index, cluster_index)
    k = k+1 #increment counting variable
  }
  mean_df$dist[i] = mean(df$dist[h:(k-1)])
  mean_df$distvector[i] = mean(df$distvector[h:(k-1)])

  h = h+z #increment outer counting variable
}

#plot
mean_df_long = gather(mean_df, key = "method", value = "rand", -p)
df_long = gather(df, key = "method", value = "rand", -p)
g = ggplot(df_long, aes(x = p, y = rand)) + 
  geom_point(aes(col = method), size = 2) +
  geom_line(data = mean_df_long, aes(x = p, y = rand, col = method), size = 1) +
  labs(title = "Cluster Performance on Noisy Data", x = "# Fake Features", y = "Adjusted Rand Index") +
  scale_color_brewer(palette = "Set1")
plot(g)

