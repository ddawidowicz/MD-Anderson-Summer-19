#Author: John Steinman
#Goal: Apply random forrest to UNC337 data to determine which genes are most important for
#      breast cancer subtype. Then use classifier to provide confidence for unsupervised clustering.
#Version: Random Forest is fitted to entrie data and then used to predict entire data for probs.
#Train Time: ~10 min
# SETUP#######################################################################
setwd("~/R/win-library/3.6/UNCSupervised") 
rm(list = ls())
library(caret)
library(nnet)
library(class)
library(randomForest)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# PREPROCESSING ##############################################################
#Load gene expression data
dat = read.table("UNCBreastRna.cdt", header = TRUE, sep = "", fill = TRUE)

#Store Pam_50
pam_50 = as.factor(as.matrix(dat[4,5:338])); rownames(pam_50) = NULL; colnames(pam_50) = NULL;
summary(pam_50)
pam_50 = t(as.data.frame(pam_50))
pam_index = as.factor(pam_50)
  
#Extract gene data
gene_dat = dat[7:nrow(dat), 5:338] 
  
#Convert to matrix
gene_dat = matrix(data = as.numeric(as.matrix(gene_dat)), 
                  nrow = nrow(gene_dat), ncol = ncol(gene_dat))
rownames(gene_dat) = dat[c(7:nrow(dat)),1]
colnames(gene_dat) = colnames(dat)[5:338]
  
#Scale Data
gene_dat = scale(gene_dat)  #center and scale
  
# UNSUPERVISED CLUSTERING ####################################################
cluster_results = hclust(dist(t(gene_dat)), method = "ward.D2")
cluster_ID = as.factor(cutree(cluster_results, 5))

# SUPERVISED PREP ############################################################
#Generate Fake Data
fake_dat = matrix(data = rnorm(20*ncol(gene_dat)), nrow = 20, ncol = ncol(gene_dat))
prefix = "Fake row"
suffix = 1:20
fake_row_names = paste(prefix, suffix, sep = " ")
rownames(fake_dat) = fake_row_names
colnames(fake_dat) = colnames(gene_dat)

#Attach Fake Data
new_dat = rbind(gene_dat, fake_dat)

#Transpose for Supervised Algorithm
t_dat = t(new_dat)

# SUPERVISED MODEL ###########################################################
#set control
myControl <- trainControl(
  method = "cv", ## cross validation
  number = 10,   ## 10-fold
  summaryFunction = multiClassSummary,
  verboseIter = FALSE)

{ #Train Model
  train.start = Sys.time()  
  model <- train(x = t_dat, y = cluster_ID,
                 method = "rf",
                 tuneLength = 15,
                 ntree = 500,
                 trControl = trainControl(method = "cv", number = 5))
  train.end = Sys.time()
  train.end - train.start  #train time
}

print(model)
plot(model)  #plots mtry vs accuracy

#Calculate variable importance
imp = varImp(model)
plot(imp)
imp_index = order(imp$importance, decreasing = TRUE)
imp_genes = as.matrix(imp$importance[imp_index,])
rownames(imp_genes) = rownames(new_dat)[imp_index]
colnames(imp_genes) = "Importance"
  
#Select dat with most important genes
fake_index = vector(length = 20)
for(i in 1:20){
  fake_index[i] = which(rownames(imp_genes) == paste(prefix, i, sep = " "))
}
fake_index = sort(fake_index, decreasing = FALSE)

select_dat = new_dat[imp_index,]
select_dat1 = select_dat[1:(fake_index[1]-1),] #keeps all genes before first fake feature
select_dat50 = select_dat[1:50,]

probability = t(model$finalModel$votes) #stores votes as confidence measure

# HEATMAPS ################################################################
gene_color = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) #colors for gene expression
label_color = c("red" ,"orange", "yellow", "green", "blue")
color_fun1 = colorRamp2(c(0, 1), c("white", "red"))
color_fun2 = colorRamp2(c(0, 1), c("white", "orange"))
color_fun3 = colorRamp2(c(0, 1), c("white", "gold"))
color_fun4 = colorRamp2(c(0, 1), c("white", "forestgreen"))
color_fun5 = colorRamp2(c(0, 1), c("white", "dodgerblue")) 
w = 30  #heatmap width

#Pam 50
pam_names = c("Normal", "Basal", "Her2", "LumB", "LumA")     
pam_colors = structure(1:5, names = pam_names)  #colors for Pam 50
ht_pam50 = Heatmap((pam_50), name = "Pam 50",col = pam_colors,
                   width = unit(w, "cm"), height = unit(0.5, "cm"))

#Traditional Heat Map
traditional_ht = Heatmap(gene_dat, name = "Gene expression", column_title = "Traditional Heat map showing gene expression" ,
                         col = gene_color, 
                         clustering_distance_columns = "euclidean",  clustering_method_columns = "ward.D2",
                         clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
                         show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                         width = unit(w, "cm"), height = unit(20, "cm"))
traditional_labels = Heatmap(t(as.character(cluster_ID)), name = "Cluster ID",
                             col = label_color, width = w, height = unit(0.5, "cm"),
                             show_column_dend = FALSE, show_row_dend = FALSE,
                             show_column_names = FALSE, show_row_names = FALSE)
draw(traditional_labels)
traditional_ht_list = traditional_ht %v% traditional_labels %v% ht_pam50
draw(traditional_ht_list)

#New Heat Map
ht = Heatmap(select_dat50, name = "Gene Expression", column_title = "Heatmap with 50 Most Important Genes" ,
             col = gene_color, 
             clustering_distance_columns = "euclidean",  clustering_method_columns = "ward.D2",
             clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
             show_row_dend = FALSE, show_row_names = TRUE, show_column_names = FALSE,
             row_names_gp = gpar(col = c(rep("red", 10), rep("black", 40))),
             width = unit(w, "cm"), height = unit(20, "cm"))
ht2 = Heatmap(t(probability[1,]), name = "Cluster 1", col = color_fun1,
              cluster_columns = FALSE, show_column_names = FALSE,
              width = unit(w, "cm"), height = unit(0.5, "cm"))
ht3 = Heatmap(t(probability[2,]), name = "Cluster 2", col = color_fun2,
              cluster_columns = FALSE, show_column_names = FALSE,
              width = unit(w, "cm"), height = unit(0.5, "cm"))
ht4 = Heatmap(t(probability[3,]), name = "Cluster 3", col = color_fun3,
              cluster_columns = FALSE, show_column_names = FALSE,
              width = unit(w, "cm"), height = unit(0.5, "cm"))
ht5 = Heatmap(t(probability[4,]), name = "Cluster 4", col = color_fun4,
              cluster_columns = FALSE, show_column_names = FALSE,
              width = unit(w, "cm"), height = unit(0.5, "cm"))
ht6 = Heatmap(t(probability[5,]), name = "Cluster 5", col = color_fun5,
              cluster_columns = FALSE, show_column_names = FALSE,
              width = unit(w, "cm"), height = unit(0.5, "cm"))

#Vertical concatentation--Note: columns will be ordered automatically
ht_list = ht %v% ht2 %v% ht3 %v% ht4 %v% ht5 %v% ht6 %v% ht_pam50
ht_list_draw = draw(ht_list)



