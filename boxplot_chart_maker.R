#!/usr/bin env

##########################################################################
# FACS - Fast AMP Clustering System
# Box plot maker from figure 2
# Author: Célio D. Santos Jr.
##########################################################################

##########################################################################
# Required libraries
##########################################################################
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")
library(ggpubr)
library(dplyr)
##########################################################################
# CMDs
##########################################################################
# Importing data
file <- read.table("scores_comparison", sep="\t", header=TRUE)
Spurious <- file$Spurious
Non_Spurious <- file$Non.spurious
##########################################################################
# Testing data
res <- t.test(Non_Spurious, Spurious, alternative = "two.sided", var.equal = FALSE)
##########################################################################
# Drawing plot
svg("boxplot.svg")
boxplot(Spurious, Non_Spurious, names=c("Spurious", "Non-Spurious"), col="gray", ylab="FACS AMP probability", ylim=c(0.5,1), outpch = 8, outcex = 1)
segments(1, 0.955, 2.0, 0.955, col = "black")
text(x = 1.5, y = 0.97, paste0("T-test, n.s. (p = ", round(res$p.value, 2),")"))
dev.off()

