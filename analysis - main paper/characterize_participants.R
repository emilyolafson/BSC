#characterize study participants
library(plyr)
initial <- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/spreadsheets/ASD_alldemog_Emily.csv")
final<- read.csv("/data/chamal/projects/emilyO/Tissue_contrast/analysis/cortical_gradient/meta-analysis/main_effect/unsmoothed/residualized/Final_gf.csv")

sum(final$Group ==  "ASD_Female")
sum(final$Group ==  "ASD_Male")
sum(final$Group ==  "Control_Female")
sum(final$Group ==  "Control_Male")

sum(initial$Group ==  "ASD_Female")
sum(initial$Group ==  "ASD_Male")
sum(initial$Group ==  "Control_Female")
sum(initial$Group ==  "Control_Male")


# F tests for Age, FIQ
abide <- subset(final$Age, final$Project=="ABIDE")
abide2 <-subset(final$Age, final$Project=="ABIDE_II")
cam <- subset(final$Age, final$Project=="Cambridge_family")
TO <- subset(final$Age, final$Project=="TORONTO_sickkids")
ukaims <-subset(final$Age, final$Project=="UKAIMS")

Data <- data.frame( Y=c(abide, abide2, cam, TO, ukaims),site =factor(rep(c("site1", "site2", "site3", "site4", "site5"), times=c(length(abide), length(abide2), length(cam), length(TO), length(ukaims)))), check.rows = F)

a <- aov(Y ~ site, data = Data)
summary(a)

abide <- subset(initial$Age, initial$Project=="ABIDE")
abide2 <-subset(initial$Age, initial$Project=="ABIDE_II")
cam <- subset(initial$Age, initial$Project=="Cambridge_family")
TO <- subset(initial$Age, initial$Project=="TORONTO_sickkids")
ukaims <-subset(initial$Age, initial$Project=="UKAIMS")

Data <- data.frame( Y=c(abide, abide2, cam, TO, ukaims),site =factor(rep(c("site1", "site2", "site3", "site4", "site5"), times=c(length(abide), length(abide2), length(cam), length(TO), length(ukaims)))), check.rows = F)

a <- aov(Y ~ site, data = Data)
summary(a)

abide <- subset(final$FIQ, final$Project=="ABIDE")
abide2 <-subset(final$FIQ, final$Project=="ABIDE_II")
cam <- subset(final$FIQ, final$Project=="Cambridge_family")
TO <- subset(final$FIQ, final$Project=="TORONTO_sickkids")
ukaims <-subset(final$FIQ, final$Project=="UKAIMS")

Data <- data.frame( Y=c(abide, abide2, cam, TO, ukaims),site =factor(rep(c("site1", "site2", "site3", "site4", "site5"), times=c(length(abide), length(abide2), length(cam), length(TO), length(ukaims)))), check.rows = F)

a <- aov(Y ~ site, data = Data)
summary(a)

abide <- subset(initial$FIQ, initial$Project=="ABIDE")
abide2 <-subset(initial$FIQ, initial$Project=="ABIDE_II")
cam <- subset(initial$FIQ, initial$Project=="Cambridge_family")
TO <- subset(initial$FIQ, initial$Project=="TORONTO_sickkids")
ukaims <-subset(initial$FIQ, initial$Project=="UKAIMS")

Data <- data.frame( Y=c(abide, abide2, cam, TO, ukaims),site =factor(rep(c("site1", "site2", "site3", "site4", "site5"), times=c(length(abide), length(abide2), length(cam), length(TO), length(ukaims)))), check.rows = F)

a <- aov(Y ~ site, data = Data)
summary(a)


abide <- subset(final$Group, final$Project=="ABIDE")
abide2 <-subset(final$Group, final$Project=="ABIDE_II")
cam <- subset(final$Group, final$Project=="Cambridge_family")
TO <- subset(final$Group, final$Project=="TORONTO_sickkids")
ukaims <-subset(final$Group, final$Project=="UKAIMS")

abide <- c(sum(abide == "Control_Male"), sum(abide=="Control_Female"), sum(abide=="ASD_Male"), sum(abide=="ASD_Female"))
abide2 <- c(sum(abide2 == "Control_Male"), sum(abide2=="Control_Female"), sum(abide2=="ASD_Male"), sum(abide2=="ASD_Female"))
cam <- c(sum(cam == "Control_Male"), sum(cam=="Control_Female"), sum(cam=="ASD_Male"), sum(cam=="ASD_Female"))
TO <- c(sum(TO == "Control_Male"), sum(TO=="Control_Female"), sum(TO=="ASD_Male"), sum(TO=="ASD_Female"))
ukaims <- c(sum(ukaims == "Control_Male"), sum(ukaims=="Control_Female"), sum(ukaims=="ASD_Male"), sum(ukaims=="ASD_Female"))


table <- cbind(as.matrix(abide), as.matrix(abide2), as.matrix(cam), as.matrix(TO), as.matrix(ukaims))

chisq.test(as.data.frame(table))




abide <- subset(final$Group, initial$Project=="ABIDE")
abide2 <-subset(initial$Group, initial$Project=="ABIDE_II")
cam <- subset(initial$Group, initial$Project=="Cambridge_family")
TO <- subset(initial$Group, initial$Project=="TORONTO_sickkids")
ukaims <-subset(initial$Group, initial$Project=="UKAIMS")

abide <- c(sum(abide == "Control_Male"), sum(abide=="Control_Female"), sum(abide=="ASD_Male"), sum(abide=="ASD_Female"))
abide2 <- c(sum(abide2 == "Control_Male"), sum(abide2=="Control_Female"), sum(abide2=="ASD_Male"), sum(abide2=="ASD_Female"))
cam <- c(sum(cam == "Control_Male"), sum(cam=="Control_Female"), sum(cam=="ASD_Male"), sum(cam=="ASD_Female"))
TO <- c(sum(TO == "Control_Male"), sum(TO=="Control_Female"), sum(TO=="ASD_Male"), sum(TO=="ASD_Female"))
ukaims <- c(sum(ukaims == "Control_Male"), sum(ukaims=="Control_Female"), sum(ukaims=="ASD_Male"), sum(ukaims=="ASD_Female"))


table <- cbind(as.matrix(abide), as.matrix(abide2), as.matrix(cam), as.matrix(TO), as.matrix(ukaims))

chisq.test(as.data.frame(table))



















cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
