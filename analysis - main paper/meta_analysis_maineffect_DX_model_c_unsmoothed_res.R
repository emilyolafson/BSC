### Load required libraries
library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 

setwd("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/")

#optional: exclude civet QC and motion fails
gf <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_alldemog_Emily.csv")
motionQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_demog_motionQC_Emily.csv")
civetQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/Total_CIVET_QC_Emily.csv")
mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')

gf$model_c_unsmoothed_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID ,"_model_c_resid_left.txt", sep ="")
gf$model_c_unsmoothed_right <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID , "_model_c_resid_right.txt", sep ="")
gf$ratio_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/ratio/unsmoothed/",gf$Subject_ID ,"_model_ratio_left_20mm_rsl.txt", sep ="")
gf$ratio_right <-  paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/ratio/unsmooothed/",gf$Subject_ID ,"_model_ratio_right_20mm_rsl.txt", sep ="")

gf <- merge(gf ,civetQC, by = "Subject_ID")
gf <- merge(gf,motionQC, by = "Subject_ID")
gf <- subset(gf,FINAL_MOTION_QC < 2)
gf <- subset(gf, FINAL_QC > 0)
setwd("/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/")
write.table(gf, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/Final_gf.csv")

gf<-subset(gf, Site_combined !='IP')
gf<-subset(gf, FIQ !='ND')
gf<-subset(gf, FIQ !='NA')
gf<-subset(gf, Age !='ND')
gf<-subset(gf, Age !='NA')

#IMPORTANT

gf$DX<-relevel(gf$DX, ref="Control")
gf$Sex_edit<-relevel(gf$Sex_edit, ref="Male")

#summary of sex/DX distribution
sum(gf$Sex_edit == "Female" & gf$DX == "ASD")
sum(gf$Sex_edit == "Female" & gf$DX == "Control")
sum(gf$Sex_edit == "Male" & gf$DX == "ASD")
sum(gf$Sex_edit == "Male" & gf$DX == "Control")

#-------------------------------------------------------------------------------------------------

##Run analysis per site, store in results_*, assuming all “Sites” in data “gf” are appropriately separate.
results_model_c_unsmoothed_left <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  if (dim(data)[1] == 0){
    next()
  }
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(data$model_c_unsmoothed_left~ DX + Age + Sex_edit + as.numeric(FIQ), data = data)
  #degrees of freedom does not get included if you immediately turn site_results into a dataframe. must extract it first.
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  site_results <- as.data.frame(site_results)
  #add df to site_results
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###Cohen’s d for DX effect
  site_results$Cohens_d <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  site_results$Project <-as.factor(unique(data$Project))
  
  results_model_c_unsmoothed_left=rbind.fill(results_model_c_unsmoothed_left, site_results)
}


meta_analysis_model_c_unsmoothed_left <- matrix(data=NA, nrow=40962, ncol=10)
for (i in 1:max(as.numeric(results_model_c_unsmoothed_left$vertex))){
  vertex=i
  data=subset(results_model_c_unsmoothed_left, vertex==i)
  
  data$yi_DX<-as.numeric(as.character(data$Cohens_d))
  data$sei_DX <- as.numeric(as.character(data$Cohens_d_se))
  DX_meta=rma(yi=yi_DX, sei=sei_DX, method="DL", data=data)
  
  DX_het=cbind(DX_meta$tau2, DX_meta$I2, DX_meta$QE, DX_meta$QEp)
  
  output=cbind(vertex,DX_meta$pval,DX_meta$b,DX_meta$se,DX_meta$ci.lb, DX_meta$ci.ub, DX_het)
  meta_analysis_model_c_unsmoothed_left[i,] <- output
  
}
colnames(meta_analysis_model_c_unsmoothed_left)=c("vertex","DX_p","DX_beta","DX_se","DX_ci_lb", "DX_ci_ub", "DX_tau2", "DX_I2", "DX_QE", "DX_QEp")


results_model_c_unsmoothed_right <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  if (dim(data)[1] == 0){
    next()
  }
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(data$model_c_unsmoothed_right~ DX + Age + Sex_edit + as.numeric(FIQ), data = data)
  #degrees of freedom does not get included if you immediately turn site_results into a dataframe. must extract it first.
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  site_results <- as.data.frame(site_results)
  #add df to site_results
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###Cohen’s d for DX effect
  site_results$Cohens_d <- { site_results$'tvalue-DXASD' * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  site_results$Project <-as.factor(unique(data$Project))
  
  results_model_c_unsmoothed_right=rbind.fill(results_model_c_unsmoothed_right, site_results)
}

###Within-disorder meta-analysis. Uses the metafor package & rma function to run random-effects meta-analysis over the entire cortex (nrow=40962 vertices). Cohen's d and corresponding standard errors as inputs—here, within disorder only.

meta_analysis_model_c_unsmoothed_right <- matrix(data=NA, nrow=40962, ncol=10)
for (i in 1:max(as.numeric(results_model_c_unsmoothed_right$vertex))){
  vertex=i
  data=subset(results_model_c_unsmoothed_right, vertex==i)
  
  data$yi_DX<-as.numeric(as.character(data$Cohens_d))
  data$sei_DX <- as.numeric(as.character(data$Cohens_d_se))
  DX_meta=rma(yi=yi_DX, sei=sei_DX, method="DL", data=data)
  
  DX_het=cbind(DX_meta$tau2, DX_meta$I2, DX_meta$QE, DX_meta$QEp)
  
  output=cbind(vertex,DX_meta$pval,DX_meta$b,DX_meta$se,DX_meta$ci.lb, DX_meta$ci.ub, DX_het)
  meta_analysis_model_c_unsmoothed_right[i,] <- output
  
}
colnames(meta_analysis_model_c_unsmoothed_right)=c("vertex","DX_p","DX_beta","DX_se","DX_ci_lb", "DX_ci_ub", "DX_tau2", "DX_I2", "DX_QE", "DX_QEp")


###Convert meta-analysis results to data.frame for easy manipulation
#eta_analysis_model_c_unsmoothed_left_1 <- as.data.frame(meta_analysis_model_c_unsmoothed_left)
meta_analysis_model_c_unsmoothed_right_1 <- as.data.frame(meta_analysis_model_c_unsmoothed_right)

#mask beta-coefficient values so that they are 0 around the subcortical structures.
#meta_analysis_model_c_unsmoothed_left_beta_masked <- meta_analysis_model_c_unsmoothed_left_1$DX_beta*mask
meta_analysis_model_c_unsmoothed_right_beta_masked <- meta_analysis_model_c_unsmoothed_right_1$DX_beta*mask

#write out masked beta coefficients for visualization.
#write.table(meta_analysis_model_c_unsmoothed_left_beta_masked,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_unsmoothed_maineffect_DX_left_beta_masked.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(meta_analysis_model_c_unsmoothed_right_beta_masked,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_unsmoothed_maineffect_DX_right_beta_masked.txt", sep = "", row.names= FALSE, col.names =FALSE)

#p values that are masked must be set to 1000 instead of 0, so that when stats are displayed on the brain, subcortical structures don't appear to be below 0.01.
meta_analysis_model_c_unsmoothed_right_p_masked <- meta_analysis_model_c_unsmoothed_right_1$DX_p*mask
#meta_analysis_model_c_unsmoothed_left_p_masked <- meta_analysis_model_c_unsmoothed_left_1$DX_p*mask
meta_analysis_model_c_unsmoothed_right_p_masked[meta_analysis_model_c_unsmoothed_right_p_masked == 0] <-1000
#meta_analysis_model_c_unsmoothed_left_p_masked[meta_analysis_model_c_unsmoothed_left_p_masked == 0] <-1000

#write out masked, FDR-adjusted p-values for visualization.
#write.table(meta_analysis_model_c_unsmoothed_left_p_masked,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effectmodel_c_unsmoothed_maineffect_DX_left_p_masked.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(meta_analysis_model_c_unsmoothed_right_p_masked,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_unsmoothed_maineffect_DX_right_p_masked.txt", sep = "", row.names= FALSE, col.names =FALSE)

###FDR correction
#meta_analysis_model_c_unsmoothed_left_1$DX_p_adjust <- p.adjust(meta_analysis_model_c_unsmoothed_left_1$DX_p, "fdr")
meta_analysis_model_c_unsmoothed_right_1$DX_p_adjust <- p.adjust(meta_analysis_model_c_unsmoothed_right_1$DX_p, "fdr")

#FDR-adjusted p values that are masked must be set to 1000 instead of 0, so that when stats are displayed on the brain, subcortical structures don't appear to be below 0.01.
meta_analysis_model_c_unsmoothed_right_p_adjust_masked <-meta_analysis_model_c_unsmoothed_right_1$DX_p_adjust*mask
#meta_analysis_model_c_unsmoothed_left_p_adjust_masked <-meta_analysis_model_c_unsmoothed_left_1$DX_p_adjust*mask
meta_analysis_model_c_unsmoothed_right_p_adjust_masked[meta_analysis_model_c_unsmoothed_right_p_adjust_masked == 0] <-1000
#meta_analysis_model_c_unsmoothed_left_p_adjust_masked[meta_analysis_model_c_unsmoothed_left_p_adjust_masked == 0] <-1000

#write out masked, FDR-adjusted p-values for visualization.
#write.table(meta_analysis_model_c_unsmoothed_left_p_adjust_masked,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_unsmoothed_maineffect_DX_left_p_adjust_masked.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(meta_analysis_model_c_unsmoothed_right_p_adjust_masked,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_unsmoothed_maineffect_DX_right_p_adjust_masked.txt", sep = "", row.names= FALSE, col.names =FALSE)


###Write out vertstats files for visualization
#mni.write.vertex.stats(meta_analysis_model_c_unsmoothed_left_1, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_left_unsmoothed_maineffect_DX.vertstats", headers=TRUE)
mni.write.vertex.stats(meta_analysis_model_c_unsmoothed_right_1, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/model_c_right_unsmoothed_maineffect_DX.vertstats", headers=TRUE)


save.image(file = "/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/meta_analysis_maineffect_DX_model_c.RData")
