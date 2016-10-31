
setwd("~/Repos/Master_Estadistica/Tema_6")
cet<-read.table("data/Cetaceans.txt",header=T)
head(cet)
str(cet)
# DolphinID es la identificación de cada animal
cet$DolphinID<-factor(cet$DolphinID)
#eliminar aquellos especímenes donde no se determinó el sexo
cet<-cet[-which(cet$Sex==0),]
cet$Sex<-factor(cet$Sex)
str(cet)
library("nlme")
m1<-gls(Age~1+Stain*Location, method="REML", data=cet)
# Otro que considera las playas.
m2<-lme(Age~1+Stain*Location, method="REML", data=cet, random=~1|Species/DolphinID)
# Otro que considera dentro de las playas, la altura de la marea.
m3<-lme(Age~1+Stain*Location, method="REML", data=cet, random=~Stain|Species/DolphinID)
AIC(m1,m2,m3)
anova(m1,m2,m3)
