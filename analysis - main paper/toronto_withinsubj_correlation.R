#Toronto within-subject analysis: no resampling.

library(RMINC)
mask <- read.table('/opt/quarantine/resources/CIVET/CIVET-CC-mask.txt')

left_curv <- dir('/data/chamal/projects/emilyO/Tissue_contrast/raw_data/TORONTO/mean_curvature', pattern ="_mc_left.txt", full.names = TRUE)
right_curv <-dir('/data/chamal/projects/emilyO/Tissue_contrast/raw_data/TORONTO/mean_curvature', pattern ="_mc_right.txt", full.names = TRUE)

left_ratio <- dir('/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/sigmoid_fit/TORONTO', pattern ="_model_ratio_left.txt", full.names = TRUE)
right_ratio <-dir('/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/sigmoid_fit/TORONTO', pattern ="_model_ratio_right.txt", full.names = TRUE)


left_curvature <- read.table("/data/chamal/projects/emilyO/Tissue_contrast/raw_data/TORONTO/mean_curvature/d8_0001_04_mc_smoothed_left.txt")
left_ratio<- read.table("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/sigmoid_fit/TORONTO/TO_d8_0001_04_model_ratio_left_smoothed.txt")

left_ratio <- left_ratio[1:40962,1]

left_curvature[left_curvature*mask ==0 ]<-NA
left_ratio[left_ratio*mask ==0 ]<-NA
  

bmp(file="/data/chamal/projects/emilyO/Tissue_contrast/figures/paper/correct_sizing/withinsubject_correlation.bmp",width=3.75, height=3, units="in", res=300)
ggplot() + geom_point(mapping = aes(x = left_curvature, y =left_ratio), alpha =0.05) + labs(x = "Mean Curvature", y = "Tissue Intensity Ratio") +theme_bw(base_size = 10)
dev.off()

