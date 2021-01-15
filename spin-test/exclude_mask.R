### Load required libraries
library(reshape2)
library(plyr)
library(RMINC)
library(metafor)
library(mni.cortical.statistics) 

column_L_ct <- mni.read.vertstats.column('/data/chamal/projects/saashi/ASD_heterogeneity/analysis/Rstats_new_Apr26-17/analysis1_DX_all/DX_all_FINAL/meta_analysis_left_thickness_DX_June29.vertstats', column.name='DX_p_adjust')
column_R_ct <- mni.read.vertstats.column('/data/chamal/projects/saashi/ASD_heterogeneity/analysis/Rstats_new_Apr26-17/analysis1_DX_all/DX_all_FINAL/meta_analysis_right_thickness_DX_June29.vertstats', column.name='DX_p_adjust')

column_L_bsc <- mni.read.vertstats.column('/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/main_effect/unsmoothed/residualized/model_c_left_unsmoothed_maineffect_DX.vertstats', column.name = 'DX_p_adjust')
column_R_bsc <- mni.read.vertstats.column('/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/main_effect/unsmoothed/residualized/model_c_right_unsmoothed_maineffect_DX.vertstats', column.name = 'DX_p_adjust')

mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')
column_L_bsc_m <-column_L_bsc*mask
column_L_bsc_m[column_L_bsc_m==0] <- NA

column_R_bsc_m <-column_R_bsc*mask
column_R_bsc_m[column_R_bsc_m==0] <- NA

column_L_ct_m <-column_L_ct*mask
column_L_ct_m[column_L_ct_m==0] <- NA

column_R_ct_m <-column_R_ct*mask
column_R_ct_m[column_R_ct_m==0] <- NA

write.table(column_L_bsc_m,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/BSC_data/BSC_left_p_adjust.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(column_R_bsc_m,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/BSC_data/BSC_right_p_adjust.txt", sep = "", row.names = FALSE, col.names = FALSE)

write.table(column_L_ct_m,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/CT_left_p_adjust.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(column_R_ct_m,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/CT_right_p_adjust.txt", sep = "", row.names = FALSE, col.names = FALSE)

#nomasking
write.table(column_L_bsc,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/BSC_data/BSC_left_p_adjust_nomask.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(column_R_bsc,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/BSC_data/BSC_right_p_adjust_nomask.txt", sep = "", row.names = FALSE, col.names = FALSE)

write.table(column_L_ct,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/CT_left_p_adjust_nomask.txt", sep = "", row.names = FALSE, col.names = FALSE)
write.table(column_R_ct,"/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/CT_right_p_adjust_nomask.txt", sep = "", row.names = FALSE, col.names = FALSE)
