# get betas from cortical thickness and BSC vertstats
library(mni.cortical.statistics)
#NEVERMIND this still doesn't work. vertstats suck.
left_beta <- mni.read.vertstats.column('/data/chamal/projects/saashi/ASD_heterogeneity/analysis/Rstats_new_Apr26-17/analysis1_DX_all/DX_all_FINAL/meta_analysis_left_thickness_DX_June29.vertstats', 'DX_beta')
right_beta <- mni.read.vertstats.column('/data/chamal/projects/saashi/ASD_heterogeneity/analysis/Rstats_new_Apr26-17/analysis1_DX_all/DX_all_FINAL/meta_analysis_right_thickness_DX_June29.vertstats', 'DX_beta')

write.csv(left_beta, '/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/left_CT_beta.csv')
write.csv(right_beta, '/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/right_CT_beta.csv')



left_beta <- mni.read.vertstats.column('/data/chamal/projects/saashi/ASD_heterogeneity/analysis/Rstats_new_Apr26-17/analysis1_DX_all/DX_all_FINAL/meta_analysis_left_thickness_DX_June29.vertstats', 'DX_beta')
right_beta <- mni.read.vertstats.column('/data/chamal/projects/saashi/ASD_heterogeneity/analysis/Rstats_new_Apr26-17/analysis1_DX_all/DX_all_FINAL/meta_analysis_right_thickness_DX_June29.vertstats', 'DX_beta')

write.csv(left_beta, '/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/left_CT_beta.csv')
write.csv(right_beta, '/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/CT_data/right_CT_beta.csv')


