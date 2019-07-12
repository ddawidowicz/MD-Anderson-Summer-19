##############################################################################
#MD Anderson Summer Project
#Recreate Figure 3 from "Statistical Significance for Hierarchical Clustering"
##############################################################################

#Load Packages 
suppressPackageStartupMessages(library("sigclust2"))
suppressPackageStartupMessages(library("sigclust2"))
library(sigclust2)
library(ComplexHeatmap)
library(BiocManager)
library(circlize)

#READ DATA #################################################################
dat<-read.table("UNCBreastRna.cdt", header = TRUE, sep = "", fill = TRUE)
dat<-dat[,1:(ncol(dat)-1)]

#Remove Headers
sub_dat<-dat[7:nrow(dat), 5:338]
sub_dat = matrix(as.numeric(as.matrix((sub_dat))), nrow = nrow(sub_dat), ncol = ncol(sub_dat))
rownames(sub_dat)<-dat[7:nrow(dat),1]
colnames(sub_dat)<-colnames(dat)[5:338]

summary(as.numeric(as.character(sub_dat)))  #5 number summary

#Store Pam_50
pam_50 = as.factor(as.matrix(dat[4,5:338])); rownames(pam_50) = NULL; colnames(pam_50) = NULL;
summary(pam_50)
pam_50 = t(as.data.frame(pam_50))

#CLEAN DATA #################################################################
dat_std = scale(sub_dat) #center and scale samples

#HEATMAPS ####################################################################
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) #colors for gene expression
pam_names = c("Normal", "Basal", "Her2", "LumB", "LumA")     
pam_colors = structure(1:5, names = pam_names)               #colors for Pam 50

ht= Heatmap(dat_std, name = "Gene Expression", col = col_fun, 
            clustering_distance_columns = "euclidean",  clustering_method_columns = "ward.D2",
            clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2",
            show_row_dend = FALSE, show_row_names = FALSE, show_column_names = FALSE,
            width = unit(14, "cm"), height = unit(10, "cm"))

ht_draw = draw(ht)

pam_50_ordered = t(pam_50[column_order(ht_draw)])
ht2 = Heatmap(pam_50_ordered, name = "Pam 50",col = pam_colors,
              width = unit(14, "cm"), height = unit(1, "cm"))
ht_list = ht %v% ht2
ht_list_draw = draw(ht_list)

#CLUSTER  ###################################################################
t = t(dat_std) #transpose data; shc performs clustering on rows
shc_result = shc(t, metric = "euclidean", linkage = "ward.D2")
plot(shc_result)

#Time Function Analysis ###################################################
test = c(50,100,250,500)
time_results=vector(mode = "numeric", length = length(test))
for (i in 1:length(test)){
  start_time = Sys.time()
  shc_results = shc(dat_std[1:test[i],1:(test[i]/5)])
  plot(shc_results)
  end_time = Sys.time()
  
  time_results[i] = end_time - start_time
}
plot(test,time_results)
