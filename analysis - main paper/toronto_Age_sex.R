### Load required libraries
library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 
#optional: exclude civet QC and motion fails
gf <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_alldemog_Emily.csv")
motionQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_demog_motionQC_Emily.csv")
civetQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/Total_CIVET_QC_Emily.csv")
mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')
gf <- merge(gf ,civetQC, by = "Subject_ID")
gf <- merge(gf,motionQC, by = "Subject_ID")
gf <- subset(gf,FINAL_MOTION_QC < 2)
gf <- subset(gf, FINAL_QC > 0)
gf <- subset(gf, Site == "TORONTO")
gf <- subset(gf, DX == "Control")

write.csv(left, '/data/chamal/projects/emilyO/Tissue_contrast/analysis/avg_BSC/left.csv')

gf$model_c_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID ,"_model_c_resid_left.txt", sep ="")
gf$model_c_right <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID , "_model_c_resid_right.txt", sep ="")

##Run analysis per site, store in results_*, assuming all “Sites” in data “gf” are appropriately separate.
results_model_c_unsmoothed_left <- data.frame()
i=1

  data=gf
  if (dim(data)[1] == 0){
    next()
  }
  N= nrow(data)
  N1=summary(data$DX)[1]  ##Control N
  N2= N - N1  ##ASD N 
  
  #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
  site_results <- vertexLm(data$model_c_unsmoothed_left~ DX*Age+ as.numeric(FIQ) + Sex_edit, data = data)
  #degrees of freedom does not get included if you immediately turn site_results into a dataframe. must extract it first.
  df <- attr(site_results, 'df')[[2]][1]
  #turn into a df to append data
  site_results <- as.data.frame(site_results)
  #add df to site_results
  site_results$df <- df
  
  #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
  ###Cohen’s d for DX effect
  site_results$Cohens_d <- { site_results$'tvalue-DXASD:Age'  * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
  site_results$Cohens_d_se <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d)**2)/8))} )
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  site_results$Project <-as.factor(unique(data$Project))
  
  results_model_c_unsmoothed_left=rbind.fill(results_model_c_unsmoothed_left, site_results)

results_model_c_unsmoothed_right <- data.frame()
i=1
  data=gf
    if (dim(data)[1] == 0){
      next()
    }
    N= nrow(data)
    N1=summary(data$DX)[1]  ##Control N
    N2= N - N1  ##ASD N 
    
    #Statistical model to be used here. Uses mni.vertex.statistics from mni.cortical.statistics package
    site_results <- vertexLm(data$model_c_unsmoothed_right~  DX*Age+ as.numeric(FIQ) +Age + Sex_edit, data = data)
    #degrees of freedom does not get included if you immediately turn site_results into a dataframe. must extract it first.
    df <- attr(site_results, 'df')[[2]][1]
    #turn into a df to append data
    site_results <- as.data.frame(site_results)
    #add df to site_results
    site_results$df <- df
    
    #Convert t-statistic of interest to Cohen's d values-- need to change tstatistic.DXASD to suit your own data
    ###Cohen’s d for DX effect
    site_results$Cohens_d <- {  site_results$'tvalue-DXASD:Age'  * (N1 + N2) } / { sqrt(N1 * N2) * sqrt(site_results$df) }
    site_results$Cohens_d_se <- sqrt( ((N1 + N2 - 1)/(N1 + N2 - 3)) * { (4/(N1 + N2)) * (1 + (((site_results$Cohens_d)**2)/8))} )
    
    site_results$vertex<-as.factor(row.names(site_results))
    site_results$Site_combined <-as.factor(site)
    site_results$Project <-as.factor(unique(data$Project))
    
    results_model_c_unsmoothed_right=rbind.fill(results_model_c_unsmoothed_right, site_results)
    
