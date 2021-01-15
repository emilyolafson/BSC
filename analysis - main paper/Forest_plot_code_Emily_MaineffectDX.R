#August8
install.packages(ggplot2)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(plyr)
library(metafor)


options(digits=10)

###LH Forest plots
load("/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/meta_analysis_maineffect_DX_model_c.RData")
final.data = results_model_c_unsmoothed_left

#L superior frontal lobe (main figure)
final.data_1 = subset(final.data, vertex == 20440)
data = final.data_1
data$yi_DX<-as.numeric(as.character(data$Cohens_d))
data$sei_DX <- as.numeric(as.character(data$Cohens_d_se))
DX_meta=rma(yi=yi_DX, sei=sei_DX, method="DL", data=data, slab=data$Site_combined)
#tiff(filename=paste(subregion, ".tiff", sep=""), width=6000, height=5500, units="px", 	compression="lzw", pointsize=20, bg= "white", res=300, antialias = "none")
bmp(file="/data/chamal/projects/emilyO/Tissue_contrast/figures/paper/correct_sizing/forestplot.bmp",width=4, height=3.75, units="in", res=300)
forplot <- forest(DX_meta, showweight = T, addfit= T, cex = 0.65, lwd=1)
dev.off()

#L superior temporal gyrus (supplementary figure)
final.data_1 = subset(final.data, vertex == 11281)
data = final.data_1
data$yi_DX<-as.numeric(as.character(data$Cohens_d))
data$sei_DX <- as.numeric(as.character(data$Cohens_d_se))
DX_meta=rma(yi=yi_DX, sei=sei_DX, method="DL", data=data, slab=data$Site_combined)
#tiff(filename=paste(subregion, ".tiff", sep=""), width=6000, height=5500, units="px", 	compression="lzw", pointsize=20, bg= "white", res=300, antialias = "none")
forplot <- forest(DX_meta, showweight = T, addfit= T, cex = 1.5, lwd=1)

#L posterior postcentral gyrus (supplementary figure)
final.data_1 = subset(final.data, vertex == 21196)
data = final.data_1
data$yi_DX<-as.numeric(as.character(data$Cohens_d))
data$sei_DX <- as.numeric(as.character(data$Cohens_d_se))
DX_meta=rma(yi=yi_DX, sei=sei_DX, method="DL", data=data, slab=data$Site_combined)
#tiff(filename=paste(subregion, ".tiff", sep=""), width=6000, height=5500, units="px", 	compression="lzw", pointsize=20, bg= "white", res=300, antialias = "none")
forplot <- forest(DX_meta, showweight = T, addfit= T, cex = 1.5, lwd=1)

###RH Forest plots
load("/data/chamal/projects/emilyO/Tissue_contrast/analysis/main_effect/meta_analysis_maineffect_DX_model_c.RData")
final.data = results_model_c_unsmoothed_right
#R posterior postcentral gyrus (supplementary figure)
final.data_1 = subset(final.data, vertex == 36467)
data = final.data_1
data$yi_DX<-as.numeric(as.character(data$Cohens_d))
data$sei_DX <- as.numeric(as.character(data$Cohens_d_se))
DX_meta=rma(yi=yi_DX, sei=sei_DX, method="DL", data=data, slab=data$Site_combined)
#tiff(filename=paste(subregion, ".tiff", sep=""), width=6000, height=5500, units="px", 	compression="lzw", pointsize=20, bg= "white", res=300, antialias = "none")
forplot <- forest(DX_meta, showweight = T, addfit= T, cex = 1.5, lwd=1)

