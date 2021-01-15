
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

#change per site.
args = commandArgs(trailingOnly=TRUE)
print(args)
gf <- merge(gf ,civetQC, by = "Subject_ID")
gf <- merge(gf,motionQC, by = "Subject_ID")
gf <- subset(gf,FINAL_MOTION_QC < 2)
gf <- subset(gf, FINAL_QC > 0)

#add all files: model_c unsmoothed, ratio values, and curvature values. all are smoothed and resampled to the common mesh.
gf$model_c_unsmoothed_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/model_c/unsmoothed/",gf$Subject_ID ,"_model_c_left_log_20mm_rsl.txt", sep ="")
gf$model_c_unsmoothed_right <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/model_c/unsmoothed/",gf$Subject_ID , "_model_c_right_log_20mm_rsl.txt", sep ="")
gf$ratio_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/ratio/unsmoothed/",gf$Subject_ID ,"_model_ratio_left_20mm_rsl.txt", sep ="")
gf$ratio_right <-  paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/ratio/unsmooothed/",gf$Subject_ID ,"_model_ratio_right_20mm_rsl.txt", sep ="")
gf$curvature_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/raw_data/mean_curvature/",gf$Subject_ID ,"_native_mc_rsl_28.28mm_white_left.txt", sep ="")
gf$curvature_right <-paste("/data/chamal/projects/emilyO/Tissue_contrast/raw_data/mean_curvature/",gf$Subject_ID ,"_native_mc_rsl_28.28mm_white_right.txt", sep ="")

site <- subset(gf, Site_combined == args)


vertexwise_response_variable_left <- vertexTable(site$model_c_unsmoothed_left)
vertexwise_response_variable_right <- vertexTable(site$model_c_unsmoothed_right)
vertexwise_covariate_left <- vertexTable(site$curvature_left)
vertexwise_covariate_right <- vertexTable(site$curvature_right)

residualized_vertexwise_data_left <-matrix(nrow=40962, ncol=ncol(vertexwise_covariate_left), 0)
residualized_vertexwise_data_right <-matrix(nrow=40962, ncol=ncol(vertexwise_covariate_left), 0)

# This loop conducts a linear model at each vertex. 
# The response variable is the distribution of values at a vertex i across your subjects, and you covariate is the distribution of values at vertex i across subjects.
# Tables of this format (rows = vertices, columns = subjects) can be easily created by calling the vertexTable function (from the RMINC package) on an array of filenames.

for (i in 1:40962){
  residualized_vertexwise_data_left[i,] <- residuals(lm(vertexwise_response_variable_left[i,] ~ vertexwise_covariate_left[i,]))
  residualized_vertexwise_data_right[i,] <- residuals(lm(vertexwise_response_variable_right[i,] ~ vertexwise_covariate_right[i,]))
}


for (i in 1:ncol(residualized_vertexwise_data_left)){
  txt <- residualized_vertexwise_data_left[,i]
  filename <- site$Subject_ID[i] # subject ID's go here
  filename <- paste(filename, "_model_c_resid_left.txt", sep = "") 
  path <- "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/"
  filename <- paste(path, filename, sep = "")
  write.table(txt,filename,sep="\t",row.names=FALSE, col.names =  FALSE)
}
for (i in 1:ncol(residualized_vertexwise_data_right)){
  txt <- residualized_vertexwise_data_right[,i]
  filename <- site$Subject_ID[i] # subject ID's go here
  filename <- paste(filename, "_model_c_resid_right.txt", sep = "") 
  path <- "/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/"
  filename <- paste(path, filename, sep = "")
  write.table(txt,filename,sep="\t",row.names=FALSE, col.names =  FALSE)
}
