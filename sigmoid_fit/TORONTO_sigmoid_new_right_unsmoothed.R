# May 5th 2019
# Can the cortical intensity profile be approximated with a sigmoid function?
### Load required libraries
##
library(reshape2)
library(plyr)
library(ggplot2)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd('/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/')
#args <- paste("TO_", args, sep="")
print(args)
gf <- data.frame(matrix(NA, nrow =1, ncol = 9))
gf$right_50 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_midright.txt", sep = "")
gf$right_25 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_25right.txt", sep = "")
gf$right_18_75 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_18_75right.txt", sep = "")
gf$right_12_5 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_12_5right.txt", sep = "")
gf$right_6_25 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_6_25right.txt", sep = "")

gf$right_whiteboundary <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_boundary_surface_right.txt", sep = "")

gf$right_wm50 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_mid_right_WM.txt", sep = "")
gf$right_wm25 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_25_right_WM.txt", sep = "")
gf$right_wm18_75 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_18_75_right_WM.txt", sep = "")
gf$right_wm12_5 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_12_5_right_WM.txt", sep = "")
gf$right_wm6_25 <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/",args,"_6_25_right_WM.txt", sep = "")

gray_right50 <- vertexTable(gf$right_50)
gray_right25 <- vertexTable(gf$right_25)
gray_right18_75 <- vertexTable(gf$right_18_75)
gray_right12_5 <- vertexTable(gf$right_12_5)
gray_right6_25 <- vertexTable(gf$right_6_25)
right_whiteboundary <- vertexTable(gf$right_whiteboundary)
white_right25 <-  vertexTable(gf$right_wm25)
white_right18_75 <-  vertexTable(gf$right_wm18_75)
white_right12_5 <-  vertexTable(gf$right_wm12_5)
white_right6_25 <-  vertexTable(gf$right_wm6_25)
sf <- rbind(t(white_right25), t(white_right18_75), t(white_right12_5), t(white_right6_25), t(right_whiteboundary), t(gray_right6_25), t(gray_right12_5), t(gray_right18_75), t(gray_right25), t(gray_right50))
df <- cbind(sf, c(-25, -18.75, -12.5, -6.25, 0, 6.25, 12.5, 18.75, 25, 50))

dff <- as.data.frame(df)
model_a <- list()
model_k <- list()
model_c<- list()
model_d <- list()
model_c_adj <- list()
model_ratio <- list()
model_max <-list()
model_min <- list()

for (i in 1:40963){ 
  xvalues <- dff[,40963]
  yvalues <- dff[,i]
  nls_fit_rest <- nls(formula = yvalues ~ a+exp(k)+(-exp(k))/(1 + exp(-c*(xvalues-d))),  start = list(a =min(yvalues), k = 5, c = 0.1, d = 0), algorithm="port", lower=c(min(yvalues), 0, 0, -50), upper=c(2000,100, 1, 50),control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/5096, printEval = FALSE, warnOnly = TRUE))
  coeffs <- coef(nls_fit_rest)
  model_c[i] <- coeffs[3]
  model_ratio[i] <- yvalues[1]/yvalues[9]
}

model_c_log<- log(as.numeric(model_c)+0.1)

 #write lists to .txt files
write.table(as.numeric(model_c), paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/sigmoid_fit/TORONTO/",args,"_model_c_right.txt", sep=""), sep = "", col.names = FALSE, row.names = FALSE)
write.table(as.numeric(model_ratio), paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/sigmoid_fit/TORONTO/",args,"_model_ratio_right.txt", sep=""), sep = "", col.names = FALSE, row.names = FALSE)
write.table(as.numeric(model_c_log), paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/sigmoid_fit/TORONTO/",args,"_model_c_right_log.txt", sep=""), sep = "", col.names = FALSE, row.names = FALSE)


