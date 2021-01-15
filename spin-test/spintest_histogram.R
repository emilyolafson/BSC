#spin test
library(ggplot2)
nullrho <- read.table('/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/nullrho_betas.txt')
realrho <- read.table('/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/realrho_betas.txt')

nullrho_qs <- read.table('/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/nullrho_combinedFDR_500pm.txt')
realrho_qs <- read.table('/data/chamal/projects/emilyO/Tissue_contrast/analysis/spin-test/realrho_combinedFDR_500pm.txt')


bmp(file="/data/chamal/projects/emilyO/Tissue_contrast/figures/paper/correct_sizing/spintest_histogram.bmp",width=3.75, height=3, units="in", res=300)
ggplot()+geom_histogram(mapping = aes(x = nullrho$V1), fill='gray', color='black', binwidth =0.01) + geom_segment(aes(x = realrho, y = 0,xend=realrho, yend=30),lineend= "round",color='red')+ theme_classic(base_size = 10)+labs(x = 'Pearson correlation coefficient', y='No. of permutations', title ='Null correlation histogram from 1000 \n permutations')
dev.off()




bmp(file="/data/chamal/projects/emilyO/Tissue_contrast/figures/paper/correct_sizing/spintest_histogram_qvals.bmp",width=7.5, height=7.5, units="in", res=300)
ggplot()+geom_histogram(mapping = aes(x = nullrho_qs$V1), fill='gray', color='black', binwidth =0.01) + geom_segment(aes(x = realrho_qs, y = 0,xend=realrho_qs, yend=20),lineend= "round",color='red')+ theme_classic(base_size = 10)+labs(x = 'Pearson correlation coefficient', y='No. of permutations', title ='Null correlation histogram from 500 \n permutations')
dev.off()


