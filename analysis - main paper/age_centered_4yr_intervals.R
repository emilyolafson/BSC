### Load required libraries
library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 

### Set default options for font & number of digits
options(digits=10)

#############
#source('/data/chamal/projects/matt/R_scripts/vertexLmer.R')
#optional: exclude civet QC and motion fails
gf <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_alldemog_Emily.csv")
motionQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_demog_motionQC_Emily.csv")
civetQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/Total_CIVET_QC_Emily.csv")
mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')
gf <- merge(gf ,civetQC, by = "Subject_ID")
gf <- merge(gf,motionQC, by = "Subject_ID")
gf <- subset(gf,FINAL_MOTION_QC < 2)
gf <- subset(gf, FINAL_QC > 0)
write.table(gf, file="gf.csv", sep=",", row.names=F)


gf$Age_4<-gf$Age-4
gf$Age_8<-gf$Age-8
gf$Age_12<- gf$Age-12
gf$Age_16<-gf$Age-16
gf$Age_20<- gf$Age-20
gf$Age_24<-gf$Age-24
gf$Age_28<-gf$Age-28
gf$Age_32<-gf$Age-32

##need to load these in. having paths in csv file does not work. 
gf$left_model_c <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID ,"_model_c_resid_left.txt", sep ="")
gf$right_model_c <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID , "_model_c_resid_right.txt", sep ="")

gf$DX<-relevel(gf$DX, ref="Control")
gf$Sex_edit<-relevel(gf$Sex_edit, ref="Male")

#-------------------------------------------------------------------------------------------------

results_left_model_c_4 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_4 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_4'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_4=rbind.fill(results_left_model_c_4, site_results)
}



meta_analysis_left_model_c_4 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_4$vertex))){
  vertex=i
  data=subset(results_left_model_c_4, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_4[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_4)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")



##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_4 <- as.data.frame(meta_analysis_left_model_c_4)
###FDR correction
meta_analysis_left_model_c_4$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_4$dx_p, "fdr")
meta_analysis_left_model_c_4$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_4$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_4, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_4-DX.vertstats", headers=TRUE)


results_right_model_c_4 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_4 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_4'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_4=rbind.fill(results_right_model_c_4, site_results)
}



meta_analysis_right_model_c_4 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_4$vertex))){
  vertex=i
  data=subset(results_right_model_c_4, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_4[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_4)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_4 <- as.data.frame(meta_analysis_right_model_c_4)
###FDR correction
meta_analysis_right_model_c_4$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_4$dx_p, "fdr")
meta_analysis_right_model_c_4$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_4$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_4, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_4-DX.vertstats", headers=TRUE)



results_right_model_c_8 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_8+ Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_8'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_8=rbind.fill(results_right_model_c_8, site_results)
}



meta_analysis_right_model_c_8 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_8$vertex))){
  vertex=i
  data=subset(results_right_model_c_8, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_8[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_8)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_8 <- as.data.frame(meta_analysis_right_model_c_8)
###FDR correction
meta_analysis_right_model_c_8$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_8$dx_p, "fdr")
meta_analysis_right_model_c_8$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_8$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_8, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_8-DX.vertstats", headers=TRUE)

results_left_model_c_8 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_8+ Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_8'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_8=rbind.fill(results_left_model_c_8, site_results)
}



meta_analysis_left_model_c_8 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_8$vertex))){
  vertex=i
  data=subset(results_left_model_c_8, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_8[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_8)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_8 <- as.data.frame(meta_analysis_left_model_c_8)
###FDR correction
meta_analysis_left_model_c_8$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_8$dx_p, "fdr")
meta_analysis_left_model_c_8$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_8$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_8, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_8-DX.vertstats", headers=TRUE)




results_left_model_c_12 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- mni.vertex.statistics(data, 'y ~ DX*Age_12 + Sex_edit', vertexdata)
  site_results <- as.data.frame(site_results)
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$tstatistic.DXASD.Age_14 * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$tstatistic.DXASD * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_12=rbind.fill(results_left_model_c_12, site_results)
}



meta_analysis_left_model_c_12 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_12$vertex))){
  vertex=i
  data=subset(results_left_model_c_12, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_12[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_12)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_12 <- as.data.frame(meta_analysis_left_model_c_12)
###FDR correction
meta_analysis_left_model_c_12$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_12$dx_p, "fdr")
meta_analysis_left_model_c_12$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_12$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_12, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_12-DX.vertstats", headers=TRUE)







results_right_model_c_12 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_12+ Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_14'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_12=rbind.fill(results_right_model_c_12, site_results)
}


meta_analysis_right_model_c_12 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_12$vertex))){
  vertex=i
  data=subset(results_right_model_c_12, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_12[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_12)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")



##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_12 <- as.data.frame(meta_analysis_right_model_c_12)
###FDR correction
meta_analysis_right_model_c_12$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_12$dx_p, "fdr")
meta_analysis_right_model_c_12$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_12$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_12, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_12-DX.vertstats", headers=TRUE)






results_left_model_c_16 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_16 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_18'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_16=rbind.fill(results_left_model_c_16, site_results)
}



meta_analysis_left_model_c_16 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_16$vertex))){
  vertex=i
  data=subset(results_left_model_c_16, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_16[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_16)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_16 <- as.data.frame(meta_analysis_left_model_c_16)
###FDR correction
meta_analysis_left_model_c_16$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_16$dx_p, "fdr")
meta_analysis_left_model_c_16$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_16$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_16, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_16-DX.vertstats", headers=TRUE)







results_right_model_c_16 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_16 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_18'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_16=rbind.fill(results_right_model_c_16, site_results)
}



meta_analysis_right_model_c_16 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_16$vertex))){
  vertex=i
  data=subset(results_right_model_c_16, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_16[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_16)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")


##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_16 <- as.data.frame(meta_analysis_right_model_c_16)
###FDR correction
meta_analysis_right_model_c_16$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_16$dx_p, "fdr")
meta_analysis_right_model_c_16$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_16$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_16, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_16-DX.vertstats", headers=TRUE)





results_left_model_c_20 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_20+ Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_20'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_20=rbind.fill(results_left_model_c_20, site_results)
}



meta_analysis_left_model_c_20 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_20$vertex))){
  vertex=i
  data=subset(results_left_model_c_20, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_20[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_20)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_20 <- as.data.frame(meta_analysis_left_model_c_20)
###FDR correction
meta_analysis_left_model_c_20$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_20$dx_p, "fdr")
meta_analysis_left_model_c_20$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_20$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_20, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_20-DX.vertstats", headers=TRUE)







results_right_model_c_20 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_20 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_20'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_20=rbind.fill(results_right_model_c_20, site_results)
}



meta_analysis_right_model_c_20 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_20$vertex))){
  vertex=i
  data=subset(results_right_model_c_20, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_20[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_20)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_20 <- as.data.frame(meta_analysis_right_model_c_20)
###FDR correction
meta_analysis_right_model_c_20$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_20$dx_p, "fdr")
meta_analysis_right_model_c_20$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_20$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_20, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_20-DX.vertstats", headers=TRUE)






results_left_model_c_24 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_24 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_24'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_24=rbind.fill(results_left_model_c_24, site_results)
}


meta_analysis_left_model_c_24 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_24$vertex))){
  vertex=i
  data=subset(results_left_model_c_24, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_24[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_24)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_24 <- as.data.frame(meta_analysis_left_model_c_24)
###FDR correction
meta_analysis_left_model_c_24$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_24$dx_p, "fdr")
meta_analysis_left_model_c_24$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_24$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_24, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_24-DX.vertstats", headers=TRUE)







results_right_model_c_24 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_24 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_24'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_24=rbind.fill(results_right_model_c_24, site_results)
}


meta_analysis_right_model_c_24 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_24$vertex))){
  vertex=i
  data=subset(results_right_model_c_24, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_24[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_24)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_24 <- as.data.frame(meta_analysis_right_model_c_24)
###FDR correction
meta_analysis_right_model_c_24$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_24$dx_p, "fdr")
meta_analysis_right_model_c_24$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_24$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_24, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_24-DX.vertstats", headers=TRUE)






results_left_model_c_28 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_28 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_28'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_28=rbind.fill(results_left_model_c_28, site_results)
}


meta_analysis_left_model_c_28 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_28$vertex))){
  vertex=i
  data=subset(results_left_model_c_28, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_28[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_28)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_28 <- as.data.frame(meta_analysis_left_model_c_28)
###FDR correction
meta_analysis_left_model_c_28$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_28$dx_p, "fdr")
meta_analysis_left_model_c_28$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_28$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_28, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_28-DX.vertstats", headers=TRUE)



results_right_model_c_28 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_28 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_28'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_28=rbind.fill(results_right_model_c_28, site_results)
}


meta_analysis_right_model_c_28 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_28$vertex))){
  vertex=i
  data=subset(results_right_model_c_28, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_28[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_28)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_28 <- as.data.frame(meta_analysis_right_model_c_28)
###FDR correction
meta_analysis_right_model_c_28$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_28$dx_p, "fdr")
meta_analysis_right_model_c_28$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_28$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_28, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_28-DX.vertstats", headers=TRUE)






results_left_model_c_32 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(left_model_c ~ DX*Age_32 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_32'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_left_model_c_32=rbind.fill(results_left_model_c_32, site_results)
}


meta_analysis_left_model_c_32 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_left_model_c_32$vertex))){
  vertex=i
  data=subset(results_left_model_c_32, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_left_model_c_32[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_left_model_c_32)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_left_model_c_32 <- as.data.frame(meta_analysis_left_model_c_32)
###FDR correction
meta_analysis_left_model_c_32$dx_p_adjust <- p.adjust(meta_analysis_left_model_c_32$dx_p, "fdr")
meta_analysis_left_model_c_32$interac_p_adjust <- p.adjust(meta_analysis_left_model_c_32$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_left_model_c_32, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_left_model_c_age_32-DX.vertstats", headers=TRUE)



results_right_model_c_32 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_model_c)
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(right_model_c ~ DX*Age_32 + Sex_edit, data)
  
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  #add df to site_results
  site_results <- as.data.frame(site_results)
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###(Changed formula for calculating Cohen’s d to equation 10 in Nakagawa & Cuthill)
  
  ###Cohen’s d for interaction effect
  site_results$Cohens_d_interac <- { site_results$'tvalue-DXASD:Age_32'* (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_interac <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) *  { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_interac)**2)/8))} )
  
  ###Cohen’s d for dx main effect
  site_results$Cohens_d_dx <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se_dx <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d_dx)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  #site_results$Project <-as.factor(unique(data$Project))
  
  results_right_model_c_32=rbind.fill(results_right_model_c_32, site_results)
}


meta_analysis_right_model_c_32 <- matrix(data=NA, nrow=40962, ncol=19)
for (i in 1:max(as.numeric(results_right_model_c_32$vertex))){
  vertex=i
  data=subset(results_right_model_c_32, vertex==i)
  
  ## Interaction (interac) 
  data$yi_interac<-as.numeric(as.character(data$Cohens_d_interac))
  data$sei_interac <- as.numeric(as.character(data$Cohens_d_se_interac))
  interac_meta=rma(yi=yi_interac, sei=sei_interac, method="DL", data=data)
  
  interac_het=cbind(interac_meta$tau2, interac_meta$I2, interac_meta$QE, interac_meta$QEp)
  
  #DX 
  data$yi_dx<-as.numeric(as.character(data$Cohens_d_dx))
  data$sei_dx <- as.numeric(as.character(data$Cohens_d_se_dx))
  dx_meta=rma(yi=yi_dx, sei=sei_dx, method="DL", data=data)
  
  dx_het=cbind(dx_meta$tau2, dx_meta$I2, dx_meta$QE, dx_meta$QEp)
  
  ###Add interac (now interac and dx) 
  output=cbind(vertex,dx_meta$pval,dx_meta$b,dx_meta$se,dx_meta$ci.lb, dx_meta$ci.ub, dx_het, interac_meta$pval, interac_meta$b, interac_meta$se, interac_meta$ci.lb, interac_meta$ci.ub, interac_het)
  meta_analysis_right_model_c_32[i,] <- output
}

###Edited to account for above interac + dx
colnames(meta_analysis_right_model_c_32)=c("vertex","dx_p","dx_beta","dx_se","dx_ci_lb", "dx_ci_ub", "dx_tau2", "dx_I2", "dx_QE", "dx_QEp","interac_p","interac_beta","interac_se","interac_ci_lb", "interac_ci_ub", "interac_tau2", "interac_I2", "interac_QE", "interac_QEp")




##Convert meta-analysis results to data.frame for easy manipulation
meta_analysis_right_model_c_32 <- as.data.frame(meta_analysis_right_model_c_32)
###FDR correction
meta_analysis_right_model_c_32$dx_p_adjust <- p.adjust(meta_analysis_right_model_c_32$dx_p, "fdr")
meta_analysis_right_model_c_32$interac_p_adjust <- p.adjust(meta_analysis_right_model_c_32$interac_p, "fdr")
###Write out files for visualization
mni.write.vertex.stats(meta_analysis_right_model_c_32, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/age_centered/meta_analysis_right_model_c_age_32-DX.vertstats", headers=TRUE)





savehistory(file = ".Rhistory")


