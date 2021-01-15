# Fit a sigmoid function to T1w intensities sampled across the cortex and into the white matter.
# Run separately on each subject.

### Load required libraries

library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 

args = commandArgs(trailingOnly=TRUE)
directory <- "/data/chamal/projects/emilyO/Tissue_contrast/derivatives/sampled_intensities/TORONTO/" #replace
resultsdir <- "/data/chamal/projects/emilyO/Tissue_contrast/results/sigmoid_fit/TORONTO/" #replace

setwd(directory)

# right hemisphere
df <- data.frame(matrix(NA, nrow =1, ncol = 9))
df$right_50 <- paste(directory,args,"_midright.txt", sep = "")
df$right_25 <- paste(directory,args,"_25right.txt", sep = "")
df$right_18_75 <- paste(directory,args,"_18_75right.txt", sep = "")
df$right_12_5 <- paste(directory,args,"_12_5right.txt", sep = "")
df$right_6_25 <- paste(directory,args,"_6_25right.txt", sep = "")
df$right_whiteboundary <- paste(directory,args,"_boundary_surface_right.txt", sep = "")
df$right_wm25 <- paste(directory,args,"_25_right_WM.txt", sep = "")
df$right_wm18_75 <- paste(directory,args,"_18_75_right_WM.txt", sep = "")
df$right_wm12_5 <- paste(directory,args,"_12_5_right_WM.txt", sep = "")
df$right_wm6_25 <- paste(directory,args,"_6_25_right_WM.txt", sep = "")

gray_right50 <- vertexTable(df$right_50)
gray_right25 <- vertexTable(df$right_25)
gray_right18_75 <- vertexTable(df$right_18_75)
gray_right12_5 <- vertexTable(df$right_12_5)
gray_right6_25 <- vertexTable(df$right_6_25)
right_whiteboundary <- vertexTable(df$right_whiteboundary)
white_right25 <-  vertexTable(df$right_wm25)
white_right18_75 <-  vertexTable(df$right_wm18_75)
white_right12_5 <-  vertexTable(df$right_wm12_5)
white_right6_25 <-  vertexTable(df$right_wm6_25)

sf <- rbind(t(white_right25), t(white_right18_75), t(white_right12_5), t(white_right6_25), t(right_whiteboundary), t(gray_right6_25), t(gray_right12_5), t(gray_right18_75), t(gray_right25), t(gray_right50))
ef <- cbind(sf, c(-25, -18.75, -12.5, -6.25, 0, 6.25, 12.5, 18.75, 25, 50))

dff <- as.data.frame(ef)
model_a <- list()
model_k <- list()
model_c<- list()
model_d <- list()
model_c_adj <- list()
model_ratio <- list()
model_max <-list()
model_min <- list()

for (i in 1:40962){ 
  xvalues <- dff[,40963]
  yvalues <- dff[,i]
  nls_fit_rest <- nls(formula = yvalues ~ a+exp(k)+(-exp(k))/(1 + exp(-c*(xvalues-d))),  start = list(a =min(yvalues), k = 5, c = 0.5, d = 0), algorithm="port", lower=c(min(yvalues), 0, 0, -50), upper=c(2000,100, 1, 50),control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/5096, printEval = FALSE, warnOnly = TRUE))
  coeffs <- coef(nls_fit_rest)
  model_c[i] <- coeffs[3]
}

model_c_log<- log(as.numeric(model_c)+0.1)

#write results to .txt files
write.table(as.numeric(model_c_log), paste(resultsdir,args,"_model_c_right_log.txt", sep=""), sep = "", col.names = FALSE, row.names = FALSE)


# left hemisphere
df <- data.frame(matrix(NA, nrow =1, ncol = 9))
df$left_50 <- paste(directory,args,"_midleft.txt", sep = "")
df$left_25 <- paste(directory,args,"_25left.txt", sep = "")
df$left_18_75 <- paste(directory,args,"_18_75left.txt", sep = "")
df$left_12_5 <- paste(directory,args,"_12_5left.txt", sep = "")
df$left_6_25 <- paste(directory,args,"_6_25left.txt", sep = "")
df$left_whiteboundary <- paste(directory,args,"_boundary_surface_left.txt", sep = "")
df$left_wm25 <- paste(directory,args,"_25_left_WM.txt", sep = "")
df$left_wm18_75 <- paste(directory,args,"_18_75_left_WM.txt", sep = "")
df$left_wm12_5 <- paste(directory,args,"_12_5_left_WM.txt", sep = "")
df$left_wm6_25 <- paste(directory,args,"_6_25_left_WM.txt", sep = "")

gray_left50 <- vertexTable(df$left_50)
gray_left25 <- vertexTable(df$left_25)
gray_left18_75 <- vertexTable(df$left_18_75)
gray_left12_5 <- vertexTable(df$left_12_5)
gray_left6_25 <- vertexTable(df$left_6_25)
left_whiteboundary <- vertexTable(df$left_whiteboundary)
white_left25 <-  vertexTable(df$left_wm25)
white_left18_75 <-  vertexTable(df$left_wm18_75)
white_left12_5 <-  vertexTable(df$left_wm12_5)
white_left6_25 <-  vertexTable(df$left_wm6_25)

sf <- rbind(t(white_left25), t(white_left18_75), t(white_left12_5), t(white_left6_25), t(left_whiteboundary), t(gray_left6_25), t(gray_left12_5), t(gray_left18_75), t(gray_left25), t(gray_left50))
ef <- cbind(sf, c(-25, -18.75, -12.5, -6.25, 0, 6.25, 12.5, 18.75, 25, 50))

dff <- as.data.frame(ef)
model_a <- list()
model_k <- list()
model_c<- list()
model_d <- list()
model_c_adj <- list()
model_ratio <- list()
model_max <-list()
model_min <- list()

for (i in 1:40962){ 
  xvalues <- dff[,40963]
  yvalues <- dff[,i]
  nls_fit_rest <- nls(formula = yvalues ~ a+exp(k)+(-exp(k))/(1 + exp(-c*(xvalues-d))),  start = list(a =min(yvalues), k = 5, c = 0.5, d = 0), algorithm="port", lower=c(min(yvalues), 0, 0, -50), upper=c(2000,100, 1, 50),control = nls.control(maxiter = 50, tol = 1e-05, minFactor = 1/5096, printEval = FALSE, warnOnly = TRUE))
  coeffs <- coef(nls_fit_rest)
  model_c[i] <- coeffs[3]
}

model_c_log<- log(as.numeric(model_c)+0.1)

#write results to .txt files
write.table(as.numeric(model_c_log), paste(resultsdir,args,"_model_c_left_log.txt", sep=""), sep = "", col.names = FALSE, row.names = FALSE)

