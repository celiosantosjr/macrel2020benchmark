#!/usr/bin env

##########################################################################
# Hemolytic model implemented in FACS testing
# Author: Célio D. Santos Jr.
##########################################################################

##########################################################################
# Required libraries
##########################################################################
if(!require(randomForest)){
  install.packages("randomForest")
  library(randomForest)
}

if(!require(caret)){
  install.packages("caret", dependencies = TRUE)
  library(caret)
}

if(!require(dplyr)){
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if(!require(data.table)){
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

##########################################################################
# CMDs
##########################################################################
#  Setting random seed
set.seed(95014)

# Taking argument
model <- readRDS("FACS/rf_dataset1.rds")

# testing model
testing <- fread(file = "HemoPI-1_Benchmark.tsv", sep="\t")
testing <- testing[,3:25]
testing$group <- as.factor(testing$group)

p <- predict(model, testing)
confusionMatrix(p, testing$group, positive = "HEMO")
varImpPlot(model)
