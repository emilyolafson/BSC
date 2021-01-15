### Load required libraries
library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 
library(ggplot2)

### Set default options for font & number of digits
options(digits=10)

#############

#setwd("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/final_pipeline_outputs/residualized/model_c/")

#optional: exclude civet QC and motion fails
gf <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_alldemog_Emily.csv")
motionQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_demog_motionQC_Emily.csv")
civetQC <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/Total_CIVET_QC_Emily.csv")
mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')

gf$model_c_unsmoothed_left <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID ,"_model_c_resid_left.txt", sep ="")
gf$model_c_unsmoothed_right <- paste("/data/chamal/projects/emilyO/Tissue_contrast/derivatives/final_pipeline_outputs/residualized/model_c/",gf$Subject_ID , "_model_c_resid_right.txt", sep ="")

gf <- merge(gf ,civetQC, by = "Subject_ID")
gf <- merge(gf,motionQC, by = "Subject_ID")
gf <- subset(gf,FINAL_MOTION_QC < 2)
gf <- subset(gf, FINAL_QC > 0)

TO <- subset(gf, Site_combined == "TORONTO")
SDSU <- subset(gf, Site_combined == "SDSU")
OHSU <- subset(gf, Site_combined == "OHSU")
UM <- subset(gf, Site_combined == "UM")
KKI <- subset(gf, Site_combined == "KKI")
NYU <- subset(gf, Site_combined == "NYU")
MAX_MUN <- subset(gf, Site_combined == "MAX_MUN")
IP <- subset(gf, Site_combined == "IP")
IoP <- subset(gf, Site_combined == "IoP")
CAM <- subset(gf, Site_combined == "CAM")
CAMBRIDGE <- subset(gf, Site_combined == "Cambridge")

##need to load these in. having paths in csv file does not work. 
gf$DX<-relevel(gf$DX, ref="Control")
gf$Sex_edit<-relevel(gf$Sex_edit, ref="Male")

#-------------------------------------------------------------------------------------------------

#Toronto 

#left
model1_TO_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , TO)
vertexFDR(model1_TO_L)
model2_TO_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit + as.numeric(FIQ) +Age, TO)
vertexFDR(model2_TO_L)
model3_TO_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit + as.numeric(FIQ) +Age, TO)
vertexFDR(model3_TO_L)
AICc_TO_left<-compare_models(model1_TO_L, model2_TO_L, model3_TO_L,metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_TO_L <- as.data.frame(AICc_TO_left)

AIC_comp_TO_L$vertex <- seq.int(nrow(AIC_comp_TO_L))

whichbest <- matrix(data=NA, nrow=40962, ncol=1)

for (i in 1:max(as.numeric(AIC_comp_TO_L$vertex))){
  whichbest[i] = which.min(AIC_comp_TO_L[i,c(1,2,3)])
}
whichbest_TO_L=whichbest

#write.table(whichbest, file="whichbest_TO_L_AIC.csv", sep=",", row.names=F)
whichbest_TO_L<-as.data.frame(whichbest_TO_L)
table(whichbest_TO_L)
#this will tell you at how many vertices each one is better

#####
#right
model1_TO_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age ,  TO)
vertexFDR(model1_TO_R)
model2_TO_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, TO)
vertexFDR(model2_TO_R)
model3_TO_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, TO)
vertexFDR(model3_TO_R)

AICc_TO_right<-compare_models(model1_TO_R, model2_TO_R,model3_TO_R, metric = AICc)


#AIC: which model is best at each vertex (which AIC is lowest) :

AIC_comp_TO_R <- as.data.frame(AICc_TO_right)

AIC_comp_TO_R$vertex <- seq.int(nrow(AIC_comp_TO_R))

whichbest <- matrix(data=NA, nrow=40962, ncol=1)

for (i in 1:max(as.numeric(AIC_comp_TO_R$vertex))){
  whichbest[i] = which.min(AIC_comp_TO_R[i,c(1,2,3)])
}
whichbest_TO_R=whichbest

#write.table(whichbest, file="whichbest_TO_R_AIC.csv", sep=",", row.names=F)

whichbest_TO_R<-as.data.frame(whichbest_TO_R)
table(whichbest_TO_R)
#whichbest_TO_R
#1     2    


#KKI
KKI = subset(gf, Site_combined=="KKI")

model1_KKI_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , KKI)
vertexFDR(model1_KKI_L)
model2_KKI_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, KKI)
vertexFDR(model2_KKI_L)
model3_KKI_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, KKI)
vertexFDR(model3_KKI_L)
AICc_KKI_left<-compare_models(model1_KKI_L, model2_KKI_L,model3_KKI_L, metric = AICc)


#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_KKI_L <- as.data.frame(AICc_KKI_left)

AIC_comp_KKI_L$vertex <- seq.int(nrow(AIC_comp_KKI_L))

whichbest <- matrix(data=NA, nrow=40962, ncol=1)

for (i in 1:max(as.numeric(AIC_comp_KKI_L$vertex))){
  whichbest[i] = which.min(AIC_comp_KKI_L[i,c(1,2,3)])
}
whichbest_KKI_L=whichbest

#write.table(whichbest, file="whichbest_KKI_L_AIC.csv", sep=",", row.names=F)

whichbest_KKI_L<-as.data.frame(whichbest_KKI_L)
table(whichbest_KKI_L)


model1_KKI_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , KKI)
vertexFDR(model1_KKI_R)
model2_KKI_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, KKI)
vertexFDR(model2_KKI_R)
model3_KKI_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, KKI)
vertexFDR(model3_KKI_R)
AICc_KKI_right<-compare_models(model1_KKI_R, model2_KKI_R, model3_KKI_R,metric = AICc)


#AIC: which model is best at each vertex (which AIC is lowest) :

AIC_comp_KKI_R <- as.data.frame(AICc_KKI_right)

AIC_comp_KKI_R$vertex <- seq.int(nrow(AIC_comp_KKI_R))

whichbest <- matrix(data=NA, nrow=40962, ncol=1)

for (i in 1:max(as.numeric(AIC_comp_KKI_R$vertex))){
  whichbest[i] = which.min(AIC_comp_KKI_R[i,c(1,2,3)])
}
whichbest_KKI_R=whichbest

#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_KKI_R_AIC.csv", sep=",", row.names=F)

whichbest_KKI_R<-as.data.frame(whichbest_KKI_R)
table(whichbest_KKI_R)




#NYU
NYU=subset(gf, Site_combined=="NYU")

model1_NYU_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , NYU)
vertexFDR(model1_NYU_L)
model2_NYU_L <- vertexLm(model_c_unsmoothed_left ~DX*Sex_edit + as.numeric(FIQ) + Age, NYU)
vertexFDR(model2_NYU_L)
model3_NYU_L <- vertexLm(model_c_unsmoothed_left ~DX + Sex_edit + as.numeric(FIQ) + Age, NYU)
vertexFDR(model3_NYU_L)
AICc_NYU_left<-compare_models(model1_NYU_L, model2_NYU_L, model3_NYU_L,metric = AICc)


#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_NYU_L <- as.data.frame(AICc_NYU_left)

AIC_comp_NYU_L$vertex <- seq.int(nrow(AIC_comp_NYU_L))

whichbest <- matrix(data=NA, nrow=40962, ncol=1)

for (i in 1:max(as.numeric(AIC_comp_NYU_L$vertex))){
  whichbest[i] = which.min(AIC_comp_NYU_L[i,c(1,2,3)])
}
whichbest_NYU_L=whichbest

#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_NYU_L_AIC.csv", sep=",", row.names=F)

whichbest_NYU_L<-as.data.frame(whichbest_NYU_L)
table(whichbest_NYU_L)

#####



model1_NYU_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , NYU)
vertexFDR(model1_NYU_R)
model2_NYU_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, NYU)
vertexFDR(model2_NYU_R)
model3_NYU_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, NYU)
vertexFDR(model3_NYU_R)
AICc_NYU_right<-compare_models(model1_NYU_R, model2_NYU_R, model3_NYU_R,metric = AICc)


#AIC: which model is best at each vertex (which AIC is lowest) :

AIC_comp_NYU_R <- as.data.frame(AICc_NYU_right)

AIC_comp_NYU_R$vertex <- seq.int(nrow(AIC_comp_NYU_R))

whichbest <- matrix(data=NA, nrow=40962, ncol=1)

for (i in 1:max(as.numeric(AIC_comp_NYU_R$vertex))){
  whichbest[i] = which.min(AIC_comp_NYU_R[i,c(1,2,3)])
}
whichbest_NYU_R=whichbest

#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_NYU_R_AIC.csv", sep=",", row.names=F)

whichbest_NYU_R<-as.data.frame(whichbest_NYU_R)
table(whichbest_NYU_R)



#OHSU
OHSU=subset(gf, Site_combined=="OHSU")

model1_OHSU_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , OHSU)
vertexFDR(model1_OHSU_L)
model2_OHSU_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, OHSU)
vertexFDR(model2_OHSU_L)
model3_OHSU_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, OHSU)
vertexFDR(model3_OHSU_L)
AICc_OHSU_left<-compare_models(model1_OHSU_L, model2_OHSU_L,model3_OHSU_L, metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_OHSU_L <- as.data.frame(AICc_OHSU_left)
AIC_comp_OHSU_L$vertex <- seq.int(nrow(AIC_comp_OHSU_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_OHSU_L$vertex))){
  whichbest[i] = which.min(AIC_comp_OHSU_L[i,c(1,2,3)])
}
whichbest_OHSU_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_OHSU_L_AIC.csv", sep=",", row.names=F)
whichbest_OHSU_L<-as.data.frame(whichbest_OHSU_L)
table(whichbest_OHSU_L)
 
#####
model1_OHSU_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , OHSU)
vertexFDR(model1_OHSU_R)
model2_OHSU_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, OHSU)
vertexFDR(model2_OHSU_R)
model3_OHSU_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, OHSU)
vertexFDR(model3_OHSU_R)

AICc_OHSU_right<-compare_models(model1_OHSU_R, model2_OHSU_R, model3_OHSU_R,metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_OHSU_R <- as.data.frame(AICc_OHSU_right)
AIC_comp_OHSU_R$vertex <- seq.int(nrow(AIC_comp_OHSU_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_OHSU_R$vertex))){
  whichbest[i] = which.min(AIC_comp_OHSU_R[i,c(1,2,3)])
}
whichbest_OHSU_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_OHSU_R_AIC.csv", sep=",", row.names=F)
whichbest_OHSU_R<-as.data.frame(whichbest_OHSU_R)
table(whichbest_OHSU_R)

#SDSU
SDSU=subset(gf, Site_combined=="SDSU")

model1_SDSU_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , SDSU)
vertexFDR(model1_SDSU_L)
model2_SDSU_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, SDSU)
vertexFDR(model2_SDSU_L)
model3_SDSU_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, SDSU)
vertexFDR(model3_SDSU_L)
AICc_SDSU_left<-compare_models(model1_SDSU_L, model2_SDSU_L,model3_SDSU_L, metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_SDSU_L <- as.data.frame(AICc_SDSU_left)
AIC_comp_SDSU_L$vertex <- seq.int(nrow(AIC_comp_SDSU_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_SDSU_L$vertex))){
  whichbest[i] = which.min(AIC_comp_SDSU_L[i,c(1,2,3)])
}
whichbest_SDSU_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_SDSU_L_AIC.csv", sep=",", row.names=F)
whichbest_SDSU_L<-as.data.frame(whichbest_SDSU_L)
table(whichbest_SDSU_L)

#####
model1_SDSU_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , SDSU)
vertexFDR(model1_SDSU_R)
model2_SDSU_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, SDSU)
vertexFDR(model2_SDSU_R)
model3_SDSU_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, SDSU)
vertexFDR(model3_SDSU_R)
AICc_SDSU_right<-compare_models(model1_SDSU_R, model2_SDSU_R, model3_SDSU_R,metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_SDSU_R <- as.data.frame(AICc_SDSU_right)
AIC_comp_SDSU_R$vertex <- seq.int(nrow(AIC_comp_SDSU_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_SDSU_R$vertex))){
  whichbest[i] = which.min(AIC_comp_SDSU_R[i,c(1,2,3)])
}
whichbest_SDSU_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_SDSU_R_AIC.csv", sep=",", row.names=F)
whichbest_SDSU_R<-as.data.frame(whichbest_SDSU_R)
table(whichbest_SDSU_R)


#UM
UM=subset(gf, Site_combined=="UM")

model1_UM_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , UM)
vertexFDR(model1_UM_L)
model2_UM_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, UM)
vertexFDR(model2_UM_L)
model3_UM_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, UM)
vertexFDR(model3_UM_L)
AICc_UM_left<-compare_models(model1_UM_L, model2_UM_L, model3_UM_L,metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_UM_L <- as.data.frame(AICc_UM_left)
AIC_comp_UM_L$vertex <- seq.int(nrow(AIC_comp_UM_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_UM_L$vertex))){
  whichbest[i] = which.min(AIC_comp_UM_L[i,c(1,2,3)])
}
whichbest_UM_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_UM_L_AIC.csv", sep=",", row.names=F)
whichbest_UM_L<-as.data.frame(whichbest_UM_L)
table(whichbest_UM_L)

#####
model1_UM_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , UM)
vertexFDR(model1_UM_R)
model2_UM_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, UM)
vertexFDR(model2_UM_R)
model3_UM_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, UM)
vertexFDR(model3_UM_R)
AICc_UM_right<-compare_models(model1_UM_R, model2_UM_R, model3_UM_R,metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_UM_R <- as.data.frame(AICc_UM_right)
AIC_comp_UM_R$vertex <- seq.int(nrow(AIC_comp_UM_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_UM_R$vertex))){
  whichbest[i] = which.min(AIC_comp_UM_R[i,c(1,2,3)])
}
whichbest_UM_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_UM_R_AIC.csv", sep=",", row.names=F)
whichbest_UM_R<-as.data.frame(whichbest_UM_R)
table(whichbest_UM_R)



#IP
IP=subset(gf, Site_combined=="IP")

model1_IP_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , IP)
vertexFDR(model1_IP_L)
model2_IP_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, IP)
vertexFDR(model2_IP_L)
model3_IP_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, IP)
vertexFDR(model3_IP_L)
AICc_IP_left<-compare_models(model1_IP_L, model2_IP_L, model3_IP_L,metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_IP_L <- as.data.frame(AICc_IP_left)
AIC_comp_IP_L$vertex <- seq.int(nrow(AIC_comp_IP_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_IP_L$vertex))){
  whichbest[i] = which.min(AIC_comp_IP_L[i,c(1,2,3)])
}
whichbest_IP_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_IP_L_AIC.csv", sep=",", row.names=F)
whichbest_IP_L<-as.data.frame(whichbest_IP_L)
table(whichbest_IP_L)

#####
model1_IP_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , IP)
vertexFDR(model1_IP_R)
model2_IP_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, IP)
vertexFDR(model2_IP_R)
model3_IP_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, IP)
vertexFDR(model3_IP_R)
AICc_IP_right<-compare_models(model1_IP_R, model2_IP_R, model3_IP_R, metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_IP_R <- as.data.frame(AICc_IP_right)
AIC_comp_IP_R$vertex <- seq.int(nrow(AIC_comp_IP_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_IP_R$vertex))){
  whichbest[i] = which.min(AIC_comp_IP_R[i,c(1,2,3)])
}
whichbest_IP_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_IP_R_AIC.csv", sep=",", row.names=F)
whichbest_IP_R<-as.data.frame(whichbest_IP_R)
table(whichbest_IP_R)

#IoP
IoP=subset(gf, Site_combined=="IoP")

model1_IoP_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , IoP)
vertexFDR(model1_IoP_L)
model2_IoP_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, IoP)
vertexFDR(model2_IoP_L)
model3_IoP_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, IoP)
vertexFDR(model3_IoP_L)
AICc_IoP_left<-compare_models(model1_IoP_L, model2_IoP_L,model3_IoP_L, metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_IoP_L <- as.data.frame(AICc_IoP_left)
AIC_comp_IoP_L$vertex <- seq.int(nrow(AIC_comp_IoP_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_IoP_L$vertex))){
  whichbest[i] = which.min(AIC_comp_IoP_L[i,c(1,2,3)])
}
whichbest_IoP_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_IoP_L_AIC.csv", sep=",", row.names=F)
whichbest_IoP_L<-as.data.frame(whichbest_IoP_L)
table(whichbest_IoP_L)

#####
model1_IoP_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , IoP)
vertexFDR(model1_IoP_R)
model2_IoP_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, IoP)
vertexFDR(model2_IoP_R)
model3_IoP_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, IoP)
vertexFDR(model3_IoP_R)
AICc_IoP_right<-compare_models(model1_IoP_R, model2_IoP_R, model3_IoP_R, metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_IoP_R <- as.data.frame(AICc_IoP_right)
AIC_comp_IoP_R$vertex <- seq.int(nrow(AIC_comp_IoP_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_IoP_R$vertex))){
  whichbest[i] = which.min(AIC_comp_IoP_R[i,c(1,2,3)])
}
whichbest_IoP_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_IoP_R_AIC.csv", sep=",", row.names=F)
whichbest_IoP_R<-as.data.frame(whichbest_IoP_R)
table(whichbest_IoP_R)


#CAM
CAM=subset(gf, Site_combined=="CAM")

model1_CAM_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , CAM)
vertexFDR(model1_CAM_L)
model2_CAM_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, CAM)
vertexFDR(model2_CAM_L)
model3_CAM_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, CAM)
vertexFDR(model3_CAM_L)
AICc_CAM_left<-compare_models(model1_CAM_L, model2_CAM_L, model3_CAM_L, metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_CAM_L <- as.data.frame(AICc_CAM_left)
AIC_comp_CAM_L$vertex <- seq.int(nrow(AIC_comp_CAM_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_CAM_L$vertex))){
  whichbest[i] = which.min(AIC_comp_CAM_L[i,c(1,2,3)])
}
whichbest_CAM_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_CAM_L_AIC.csv", sep=",", row.names=F)
whichbest_CAM_L<-as.data.frame(whichbest_CAM_L)
table(whichbest_CAM_L)

#####
model1_CAM_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , CAM)
vertexFDR(model1_CAM_R)
model2_CAM_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, CAM)
vertexFDR(model2_CAM_R)
model3_CAM_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, CAM)
vertexFDR(model3_CAM_R)
AICc_CAM_right<-compare_models(model1_CAM_R, model2_CAM_R, model3_CAM_R, metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_CAM_R <- as.data.frame(AICc_CAM_right)
AIC_comp_CAM_R$vertex <- seq.int(nrow(AIC_comp_CAM_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_CAM_R$vertex))){
  whichbest[i] = which.min(AIC_comp_CAM_R[i,c(1,2,3)])
}
whichbest_CAM_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_CAM_R_AIC.csv", sep=",", row.names=F)
whichbest_CAM_R<-as.data.frame(whichbest_CAM_R)
table(whichbest_CAM_R)

#CAMBRIDGE
model1_CAMBRIDGE_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , CAMBRIDGE)
vertexFDR(model1_CAMBRIDGE_L)
model2_CAMBRIDGE_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, CAMBRIDGE)
vertexFDR(model2_CAMBRIDGE_L)
model3_CAMBRIDGE_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, CAMBRIDGE)
vertexFDR(model3_CAMBRIDGE_L)
AICc_CAMBRIDGE_left<-compare_models(model1_CAMBRIDGE_L, model2_CAMBRIDGE_L, model3_CAMBRIDGE_L, metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_CAMBRIDGE_L <- as.data.frame(AICc_CAMBRIDGE_left)
AIC_comp_CAMBRIDGE_L$vertex <- seq.int(nrow(AIC_comp_CAMBRIDGE_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_CAMBRIDGE_L$vertex))){
  whichbest[i] = which.min(AIC_comp_CAMBRIDGE_L[i,c(1,2,3)])
}
whichbest_CAMBRIDGE_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_CAMBRIDGE_L_AIC.csv", sep=",", row.names=F)
whichbest_CAMBRIDGE_L<-as.data.frame(whichbest_CAMBRIDGE_L)
table(whichbest_CAMBRIDGE_L)

#####
model1_CAMBRIDGE_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , CAMBRIDGE)
vertexFDR(model1_CAMBRIDGE_R)
model2_CAMBRIDGE_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, CAMBRIDGE)
vertexFDR(model2_CAMBRIDGE_R)
model3_CAMBRIDGE_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, CAMBRIDGE)
vertexFDR(model3_CAMBRIDGE_R)
AICc_CAMBRIDGE_right<-compare_models(model1_CAMBRIDGE_R, model2_CAMBRIDGE_R, model3_CAMBRIDGE_R, metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_CAMBRIDGE_R <- as.data.frame(AICc_CAMBRIDGE_right)
AIC_comp_CAMBRIDGE_R$vertex <- seq.int(nrow(AIC_comp_CAMBRIDGE_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_CAMBRIDGE_R$vertex))){
  whichbest[i] = which.min(AIC_comp_CAMBRIDGE_R[i,c(1,2,3)])
}
whichbest_CAMBRIDGE_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_CAMBRIDGE_R_AIC.csv", sep=",", row.names=F)
whichbest_CAMBRIDGE_R<-as.data.frame(whichbest_CAMBRIDGE_R)
table(whichbest_CAMBRIDGE_R)



#MAX_MUN
model1_MAX_MUN_L <- vertexLm(model_c_unsmoothed_left ~ DX +as.numeric(FIQ) +Age , MAX_MUN)
vertexFDR(model1_MAX_MUN_L)
model2_MAX_MUN_L <- vertexLm(model_c_unsmoothed_left ~ DX*Sex_edit +as.numeric(FIQ)  + Age, MAX_MUN)
vertexFDR(model2_MAX_MUN_L)
model3_MAX_MUN_L <- vertexLm(model_c_unsmoothed_left ~ DX+Sex_edit +as.numeric(FIQ)  + Age, MAX_MUN)
vertexFDR(model3_MAX_MUN_L)
AICc_MAX_MUN_left<-compare_models(model1_MAX_MUN_L, model2_MAX_MUN_L, model3_MAX_MUN_L, metric = AICc)
#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_MAX_MUN_L <- as.data.frame(AICc_MAX_MUN_left)
AIC_comp_MAX_MUN_L$vertex <- seq.int(nrow(AIC_comp_MAX_MUN_L))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_MAX_MUN_L$vertex))){
  whichbest[i] = which.min(AIC_comp_MAX_MUN_L[i,c(1,2,3)])
}
whichbest_MAX_MUN_L=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_MAX_MUN_L_AIC.csv", sep=",", row.names=F)
whichbest_MAX_MUN_L<-as.data.frame(whichbest_MAX_MUN_L)
table(whichbest_MAX_MUN_L)

#####
model1_MAX_MUN_R <- vertexLm(model_c_unsmoothed_right ~ DX +as.numeric(FIQ) +Age , MAX_MUN)
vertexFDR(model1_MAX_MUN_R)
model2_MAX_MUN_R <- vertexLm(model_c_unsmoothed_right ~ DX*Sex_edit +as.numeric(FIQ)  + Age, MAX_MUN)
vertexFDR(model2_MAX_MUN_R)
model3_MAX_MUN_R <- vertexLm(model_c_unsmoothed_right ~ DX+Sex_edit +as.numeric(FIQ)  + Age, MAX_MUN)
vertexFDR(model3_MAX_MUN_R)
AICc_MAX_MUN_right<-compare_models(model1_MAX_MUN_R, model2_MAX_MUN_R, model3_MAX_MUN_R, metric = AICc)

#AIC: which model is best at each vertex (which AIC is lowest) :
AIC_comp_MAX_MUN_R <- as.data.frame(AICc_MAX_MUN_right)
AIC_comp_MAX_MUN_R$vertex <- seq.int(nrow(AIC_comp_MAX_MUN_R))
whichbest <- matrix(data=NA, nrow=40962, ncol=1)
for (i in 1:max(as.numeric(AIC_comp_MAX_MUN_R$vertex))){
  whichbest[i] = which.min(AIC_comp_MAX_MUN_R[i,c(1,2,3)])
}
whichbest_MAX_MUN_R=whichbest
#write.table(whichbest, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/AIC/whichbest_MAX_MUN_R_AIC.csv", sep=",", row.names=F)
whichbest_MAX_MUN_R<-as.data.frame(whichbest_MAX_MUN_R)
table(whichbest_MAX_MUN_R)


##rename column headers of "which best" files, then merge together (cbind) to get one big csv of which model is better at each vertex, for each site
#example: names(DX_L) <- c("dx_p")

#Sites: CAM, Cambridge, IoP, IP, KKI, MAX_MUN, NIMH, NYU, OHSU, SDSU, TORONTO, UM

names(whichbest_KKI_L) <- c("KKI")
names(whichbest_NYU_L) <- c("NYU")
names(whichbest_OHSU_L) <- c("OHSU")
names(whichbest_TO_L) <- c("TO")
names(whichbest_SDSU_L) <- c("SDSU")
names(whichbest_CAM_L) <- c("CAM")
names(whichbest_CAMBRIDGE_L) <- c("CAMBRIDGE")
names(whichbest_UM_L) <- c("UM")
names(whichbest_MAX_MUN_L) <- c("MAX_MUN")
names(whichbest_IP_L) <- c("IP")
names(whichbest_IoP_L) <- c("IoP")


whichbest_L_combined=cbind(whichbest_KKI_L,whichbest_NYU_L,whichbest_OHSU_L, whichbest_TO_L, whichbest_SDSU_L, whichbest_CAM_L, whichbest_CAMBRIDGE_L ,whichbest_UM_L, whichbest_MAX_MUN_L,  whichbest_IP_L, whichbest_IoP_L)

#R
names(whichbest_KKI_R) <- c("KKI")
names(whichbest_NYU_R) <- c("NYU")
names(whichbest_OHSU_R) <- c("OHSU")
names(whichbest_TO_R) <- c("TO")
names(whichbest_SDSU_R) <- c("SDSU")
names(whichbest_CAM_R) <- c("CAM")
names(whichbest_CAMBRIDGE_R) <- c("CAMBRIDGE")
names(whichbest_UM_R) <- c("UM")
names(whichbest_MAX_MUN_R) <- c("MAX_MUN")
names(whichbest_IP_R) <- c("IP")
names(whichbest_IoP_R) <- c("IoP")

whichbest_R_combined=cbind(whichbest_KKI_R,whichbest_NYU_R,whichbest_OHSU_R, whichbest_TO_R, whichbest_SDSU_R, whichbest_CAM_R, whichbest_CAMBRIDGE_R ,whichbest_UM_R, whichbest_MAX_MUN_R,  whichbest_IP_R, whichbest_IoP_R)

write.table(whichbest_L_combined, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/revisions/AIC_sex/sex_whichbest_L_combined_model_c.csv", sep=",", row.names=F)
write.table(whichbest_R_combined, file="/data/chamal/projects/emilyO/Tissue_contrast/analysis/revisions/AIC_sex/sex_whichbest_R_combined_model_c.csv", sep=",", row.names=F)

##


