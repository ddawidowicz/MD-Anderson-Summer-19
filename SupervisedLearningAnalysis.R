#Set Working Directory
setwd("/Users/620jo/OneDrive/Documents/R/win-library/3.6/BreastCancerMachineLearning")

#libraries
library(class)
library(caret)
library(randomForest)
library(BBmisc)
rm(list = ls())

#Import Data
data = read.csv(file = "breast-cancer-wisconsin.data", header = FALSE, sep = ",")
col_titles = c("ID", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "y")
colnames(data) = col_titles

#Convert Y to binomial form
for(i in 1:length(data$y)){
  if(data$y[i] == 2){ 
    data$y[i] = "B"} #Benign Tumor
  else{
    data$y[i] = "M"  #Malignant Tumor
  }
}
data$y = factor(data$y)

#Split data
tr <- sample(nrow(data), round(nrow(data) * 0.6))
train = data[tr, ]
test = data[-tr, ]

#Fit Models
myControl <- trainControl(
  method = "cv", ## cross validation
  number = 10,   ## 10-fold
  summaryFunction = twoClassSummary, ## NEW
  classProbs = TRUE, 
  verboseIter = FALSE
)

#K nearest Neighbors
knn_model <- train(y ~ ., test,
                   method = "knn", ## to use glm's logistic regression
                   trControl = myControl)
confusionMatrix(knn_model)
print(knn_model)
plot(knn_model)

#Generalized linear model
glm_model <- train(y ~ ., test,
               method = "glm", ## to use glm's logistic regression
               trControl = myControl)
print(glm_model)

p = predict(glm_model, test, type = "prob")
confusionMatrix(glm_model)

#Bayes Generalized linear model
bayesglm_model <- train(y ~ ., test,
                   method = "bayesglm", ## to use glm's logistic regression
                   trControl = myControl)

# Random Forest
rf_model <- train(y ~ ., test,
                   method = "rf", ## to use glm's logistic regression
                   trControl = myControl)
confusionMatrix(rf_model)
print(rf_model)

