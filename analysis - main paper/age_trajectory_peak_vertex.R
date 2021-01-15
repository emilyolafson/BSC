#### Load required libraries
library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 
library(ggplot2)

setwd("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/")

#optional: exclude civet QC and motion fails
gf <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_alldemog_Emily.csv")
motionQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_demog_motionQC_Emily.csv")
civetQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/Total_CIVET_QC_Emily.csv")
mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')


gf$model_c_unsmoothed_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID ,"_model_c_resid_left.txt", sep ="")
gf$model_c_unsmoothed_right <- paste("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID , "_model_c_resid_right.txt", sep ="")

gf <- merge(gf ,civetQC, by = "Subject_ID")
gf <- merge(gf,motionQC, by = "Subject_ID")
gf <- subset(gf,FINAL_MOTION_QC < 2)
gf <- subset(gf, FINAL_QC > 0)


#right hemisphere: large differnces in superior temporal cortex, vertex #11637
model_c <- vertexTable(gf$model_c_unsmoothed_right)
bsc <- model_c[22534,]
age <- gf$Age
dx <- gf$DX

df<- as.data.frame(cbind(bsc, age, dx))

bmp(file="/data/chamal/projects/emilyO/Tissue_contrast/figures/paper/correct_sizing/age_trajectory.bmp",width=4.75, height=3, units="in", res=300)
ggplot(mapping=aes(x = age, y =bsc), data = df)+geom_jitter(aes(x = age, y = bsc, fill=dx, color=dx))+geom_smooth(aes(x = age, y = bsc, group=dx, color=dx), method="lm")+ theme(legend.position="none")+theme_classic(base_size = 10)+labs(x='Age', y = "Boundary Sharpness Coefficient at a peak \n vertex in the Right Superior Temporal lobe")
dev.off()
