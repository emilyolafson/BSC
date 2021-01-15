library(ggplot2)
library(RMINC)
library(lme4)
library(doBy)
library(metafor)
library(mni.cortical.statistics) 

gf <- read.csv('/data/chamal/projects/saashi/Sex_diffs2/ASD_combined_total.csv')

QC_motion <- read.csv('/data/chamal/projects/saashi/Sex_diffs2/MOTION_QC/final_QC_tomerge/FINAL_MOTION_QC.csv')
QC_CIVET <- read.csv('/data/chamal/projects/saashi/CIVET_VERIFIES/comparison_CSVs/comparison/Total_CIVET_QC.csv')
gf <- merge(gf, QC_motion, by="Subject_ID")
gf <-merge(gf, QC_CIVET, by="Subject_ID")
gf = subset(gf, FINAL_MOTION_QC <2)
gf = subset(gf, FINAL_QC>0)
gf = subset(gf, Site_combined!="CMU" & Site_combined!="GU" & Site_combined!="OLIN" & Site_combined!="PITT" & Site_combined!="UCD" & Site_combined!="UCLA" & Site_combined!="IU")

#ADOS2
gf_ADOS2 <- subset(gf, gf$ADOS_2_TOTAL!='NA')
gf_ADOS2 <- subset(gf_ADOS2, ADOS_2_TOTAL>0)
gf_ADOS2 <-subset(gf_ADOS2, DX=="ASD" )

gf_ADOS2$Severity <- gf_ADOS2$ADOS_2_TOTAL

gf=gf_ADOS2

#module 3 = 143, module 4 = 15 - so run only only those with module 3 
gf_ADOS2_3 = subset(gf, ADOS_MODULE==3)

gf= gf_ADOS2_3

##need to load these in. having paths in csv file does not work. 
gf$left_thickness <- paste("/data/chamal/projects/saashi/Sex_diffs2/CIVET_output_links_updated/",gf$Subject_ID, "_native_rms_rsl_tlink_28.28mm_left.txt", sep="")
gf$right_thickness <- paste("/data/chamal/projects/saashi/Sex_diffs2/CIVET_output_links_updated/",gf$Subject_ID,"_native_rms_rsl_tlink_28.28mm_right.txt", sep="")
gf$left_surface <- paste("/data/chamal/projects/saashi/Sex_diffs2/CIVET_output_links_updated/",gf$Subject_ID,"_mid_surface_rsl_left_native_area_56.57mm.txt", sep="")
gf$right_surface <- paste("/data/chamal/projects/saashi/Sex_diffs2/CIVET_output_links_updated/",gf$Subject_ID,"_mid_surface_rsl_right_native_area_56.57mm.txt", sep="")

gf$Sex_edit<-relevel(gf$Sex_edit, ref="Male")

#------------------------------------------------------------------------------------------


##Run vertexLm per site, store in results_*
#For ADOS 2 Total Score 

results_left_thickness_1 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$left_thickness)
  
  #Statistical model to be used here. Uses vertexLm
  site_results <- vertexLm(left_thickness ~ Severity + Age + Sex_edit, data)
  site_results <- as.data.frame(site_results)
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  
  results_left_thickness_1=rbind(results_left_thickness_1, site_results)
}

results_right_thickness_1 <- data.frame()
for (i in 1: length(unique(gf$Site_combined))){
  site= as.character(unique(gf$Site_combined))[i]
  data=subset(gf, Site_combined==site)
  
  vertexdata <- vertexTable(data$right_thickness)
  
  #Statistical model to be used here. Uses vertexLm
  site_results <- vertexLm(right_thickness ~ Severity + Age + Sex_edit, data)
  site_results <- as.data.frame(site_results)
  
  site_results$vertex<-as.factor(row.names(site_results))
  site_results$Site_combined <-as.factor(site)
  
  results_right_thickness_1=rbind(results_right_thickness_1, site_results)
}




#meta analysis part - escalc to get effect sizes (partial correlation) and then rma 

#Left thickness
mysummary = summaryBy(left_thickness ~ Site_combined, data=gf, FUN=length)
colnames(mysummary) = c("Site_combined", "length")

stats = merge(results_left_thickness_1, mysummary, by="Site_combined")

#Setup a cluster and send variables to it
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(metafor))
clusterExport(cl, "stats")

#use escalc to calculate the effect size to be used (here, the partial correlation coefficient (Fisher’s r-to-z transformed ) - using ZPCOR)
yivi = parSapply(cl, 1:40962 , function(x) escalc(measure="SPCOR", ti=stats[stats$vertex==x,]$`tvalue-Severity`, ni=stats[stats$vertex==x,]$length, r2i=stats[stats$vertex==x,]$`R-squared`, mi=rep(3, length(stats[stats$vertex==x,]$length)) ), simplify=F)

stopCluster(cl)
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(metafor))
clusterExport(cl, "yivi")

rmaresults =  parSapply(cl, 1:40962 , function(x) rma(yi,vi, data = yivi[[x]], method="DL"), simplify=F)

stopCluster(cl)

pvalues = sapply(1:40962, function(x) rmaresults[[x]]$pval, simplify="array")
padjusted = p.adjust(pvalues, "fdr")



#right thickness
mysummary_r = summaryBy(right_thickness ~ Site_combined, data=gf, FUN=length)
colnames(mysummary_r) = c("Site_combined", "length")

stats_r = merge(results_right_thickness_1, mysummary_r, by="Site_combined")

#Setup a cluster and send variables to it
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(metafor))
clusterExport(cl, "stats_r")

#use escalc to calculate the effect size to be used (here, the partial correlation coefficient (Fisher’s r-to-z transformed ) - using ZPCOR)
yivi_r = parSapply(cl, 1:40962 , function(x) escalc(measure="SPCOR", ti=stats_r[stats_r$vertex==x,]$`tvalue-Severity`, ni=stats_r[stats_r$vertex==x,]$length, r2i=stats_r[stats_r$vertex==x,]$`R-squared`, mi=rep(3, length(stats_r[stats_r$vertex==x,]$length)) ), simplify=F)

stopCluster(cl)
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(metafor))
clusterExport(cl, "yivi_r")

rmaresults_r =  parSapply(cl, 1:40962 , function(x) rma(yi,vi, data = yivi_r[[x]], method="DL"), simplify=F)

stopCluster(cl)

pvalues_r = sapply(1:40962, function(x) rmaresults_r[[x]]$pval, simplify="array")
padjusted_r = p.adjust(pvalues_r, "fdr")

save.image()



mni.write.vertex.stats(padjusted, "meta_analysis_left_thickness_SeverityADOS_2_total.vertstats", headers=TRUE)
mni.write.vertex.stats(pvalues, "meta_analysis_left_thickness_SeverityADOS_2_total_UNCORRECTED.vertstats", headers=TRUE)

mni.write.vertex.stats(padjusted_r, "meta_analysis_right_thickness_SeverityADOS_2_total.vertstats", headers=TRUE)
mni.write.vertex.stats(pvalues_r, "meta_analysis_right_thickness_SeverityADOS_2_total_UNCORRECTED.vertstats", headers=TRUE)


