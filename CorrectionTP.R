library(skimr)
# Importation des données
dt = read.csv("Data/Inhibiteurs_GPR54.dat (copie).txt", sep = "\t")

# Description des données
summary(dt)
skimr::skim(dt)
sum(is.na(dt))

#Nettoyage
dt[which(dt=="", arr.ind = T)]= NA
dt1=dt
dt1$N<-NULL
dt1$X<-NULL
dt1$StructureChemDraw<-NULL
dt1 = dt1[-which(is.na(dt1$C.ter)),]

#Check NA des seq
sum(is.na(dt1[,5:16]))

mydist = dist(dt1[,5:16], method = "manhattan")
# MAS QUE FAAAAAAAAAAAAAAAIIIIIIRE
