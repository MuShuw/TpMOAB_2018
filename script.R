setwd("~/Bureau/projet_moab/")

# ---- Nettoyage ----
dt = read.csv("Inhibiteurs_GPR54.dat (copie).txt", sep="\t", dec = ".")
dt = dt[,-c(19,18,4)]
dt = dt[-64,]
dt[c(78, 51, 35, 5), 4] = NA
write.csv(dt, file = "Inhibiteurs_GPR54.dat(copie-Bertrand).txt", row.names = FALSE)



# ---- Visualisation des donnees ----
dt[,2] = as.factor(dt[,2])

# probleme conversion SD en numerice
test = as.numeric(dt[-which(is.na(dt[,4])),4])
dim(dt)
summary(dt)

# Suppression ligne 64 car pas de val de mean donc aucun moyen d'évaluer sa liaison au recepteur
dt = dt[-64,]

par(mfrow = c(1,2))
plot(dt[,3], main = "MEANS")
plot(dt[,4], main = "SD")
boxplot(dt[,3], main = "MEANS")
boxplot(dt[,4], main = "SD")
# on ne supprime pas les valeures "extremes" car ce n'en est pas, c'est juste que si pas liaison au recepteur -> 
# proche de 0 et des qu'il y a une activite valeur tres superieur a 0

length(which(dt$Means >= 1))
dt[which(dt$Means >= 1), 3]
# donc 33 composé potentiellement capable de se lier au recepteur

# autre moyen de visualisation
library(FactoMineR)
dt.acp = PCA(dt[,3:4])
dev.off()



# ---- Classification des donnees ----
# Partitionnement k-means
I.intra = vector(length=15)
for (i in 1:15) {
  I.intra[i] = kmeans(dt[,3], centers=i)$tot.withinss
}
plot(I.intra, type="b")
dt.kmeans = kmeans(dt[,3], centers=5, nstart = 3)
plot(dt$Means, col= c(1,2,3,4,5,6)[as.factor(dt.kmeans$cluster)])
# on voit la séparation du jeu en groupes qui correspondent à des palier d'activité ici 6 palier, faudra 
# visualiser (rapport 2page) les ressemblance dans ces deux groupes et différence entre les groupes au niveau des séquences

#### Version en prenant Means et SD (perte composé avec liaison vu l'activit donc rentable ????) ####
# which(dt[,4] == "NA")
# dt.1 = dt[-c(77, 51, 35, 5),]
# summary(dt.1)
# 
# I.intra = vector(length=15)
# for (i in 1:15) {
#   I.intra[i] = kmeans(dt.1[,3:4], centers=i)$tot.withinss
# }
# plot(I.intra, type="b")
# 
# dt.1.kmeans = kmeans(dt.1[,3:4], centers=3, nstart = 3)
# plot(dt.1$Means, col= c(1,2,3,4,5,6)[as.factor(dt.1.kmeans$cluster)])


# Classification hiérarchique ascendante
dt.hclust = hclust(dist(dt[,3]), method = "ward.D2")
plot(dt.hclust, hang = -0.1, labels = F, main = "Cluster Dendrogram expression gene")
rect.hclust(dt.hclust, k = 2)
rect.hclust(dt.hclust, k = 3)
rect.hclust(dt.hclust, k = 5)
dt.cut = as.factor(cutree(dt.hclust, k = 5))

table(dt.cut, as.factor(dt.kmeans$cluster))
par(mfrow=c(1,2))
plot(dt$Means, col= c(1,2,3,4,5)[as.factor(dt.kmeans$cluster)], main = "Classification k-means")
plot(dt$Means, col= c(2,1,3,4,5)[dt.cut], main = "Classification hclust")
# relativement bonne corrélation entre hclust et kmeansexcepté pour deux individu (les plus proche de 0)
# qui sont mis dans le groupe inactif par hclust


# vu qu'on veut obtenir tendance sur la (ou les) position(s) essentielle(s) pour correctement se lier 
# au récepteur, on voudra préfèrera être restrictif pour être sur sûr que la modification de position
# est réellement nécessaire à la liaison au recepteur, à partir de là on "prefere" donc la classification
# de hclust
dt = cbind(dt[,1:4], dt.cut, dt[,5:16])
colnames(dt)[5] = "classe"
dt



# ---- Classification des donnees ----




















