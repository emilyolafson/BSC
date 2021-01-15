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
gf <- subset(gf, Project == "UKAIMS")
gf <- subset(gf, gf$Age > 22)
gf <- subset(gf, gf$Age < 35)

gf <- subset(gf, DX == "Control")

write.csv(left, '/data/chamal/projects/emilyO/Tissue_contrast/analysis/avg_BSC/left.csv')

gf$model_c_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID ,"_model_c_resid_left.txt", sep ="")
gf$model_c_right <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID , "_model_c_resid_right.txt", sep ="")

left <-vertexTable(gf$model_c_left)
right <- vertexTable(gf$model_c_right)

left_avg <- rowMeans(left, na.rm=T)
right_avg <- rowMeans(right, na.rm=T)

left_avg_masked[left_avg*mask==0] <- -1000 
right_avg_masked[right_avg*mask==0] <- -1000 

write.table(left_avg_masked, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/avg_BSC/left_HCP_age_range_UKAIMS_residualized_masked.txt", sep = "", row.names = F, col.names = F)
write.table(right_avg_masked, "/data/chamal/projects/emilyO/Tissue_contrast/analysis/avg_BSC/right_HCP_age_range_UKAIMS_residualized_masked.txt", sep = "", row.names = F, col.names = F)
