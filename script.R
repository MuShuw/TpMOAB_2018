setwd("~/Bureau/projet_moab/")

# ---- Nettoyage ----
# Suppression pour fournir le jeu de donnée du TP, retire les colones inutiles, met de NA quand une valeur est manquante et
# suppression de la ligne 64 qui contient un peptide sans moyenne donc inutilisable
# suppression de la ligne 70 qui contient un peptide sans sequence donc inutilisable
dt = read.csv("Inhibiteurs_GPR54.dat (copie).txt", sep="\t", dec = ".")
dt = dt[,-c(19,18,4)]
dt = dt[-c(70,64),]
dt[c(78, 51, 35, 5), 4] = NA
write.csv(dt, file = "Inhibiteurs_GPR54.dat(copie-Bertrand).txt", row.names = FALSE)



# ---- Import des donnees ----
dt = read.csv("Inhibiteurs_GPR54.dat(copie-Bertrand).txt")



# ---- Visualisation des donnees ----
# 1.
dim(dt)
# 80 lignes et 16 colones
summary(dt)
# Necessite de transformer la colone N en variable qualitative
dt[,2] = as.factor(dt[,2])

# 2.
par(mfrow = c(1,2))
plot(dt[,3], main = "MEANS")
plot(dt[,4], main = "SD")

# 3.
# On observe la présence d'une base de peptides donc l'activité est quasiment nulle et quelque peptide dont 
# l'activité est très supérieure à 0. On ne supprime pas les valeures "extremes" car ce n'en est pas, c'est 
# seulement le fait que si il n'y a pas liaison au recepteur, la moyenne sera proche de 0 et des qu'il y a 
# une activite, la valeur sera tres superieur a 0



# ---- Classification des donnees ----
# Partitionnement k-means
# 1.
I.intra = vector(length=15)
for (i in 1:15) {
  I.intra[i] = kmeans(dt[,3], centers=i)$tot.withinss
}
plot(I.intra, type="b")
# On veut minimiser l'inertie intra en minimisant également le nombre de groupe, au vue de la courbe de 
# l'inertie intra, le nombre optimal de groupe semble etre 5 groupes

# 2.
dt.kmeans = kmeans(dt[,3], centers=5, nstart = 3)

# 3.
plot(dt$Means, col= c(1,2,3,4,5,6)[as.factor(dt.kmeans$cluster)])
# on voit la séparation du jeu en groupes qui correspondent à des palier d'activité ici 5 paliers
# avec, pour le groupe proche de 0, les peptides sans activite donc sans liaison

# Classification hiérarchique ascendante
# 1.
dt.dist = dist(dt[,3])

# 2.
dt.hclust = hclust(dt.dist, method = "ward.D2")

# 3.
plot(dt.hclust, hang = -0.1, labels = F, main = "Cluster Dendrogram expression gene")
rect.hclust(dt.hclust, k = 2)
rect.hclust(dt.hclust, k = 3)
rect.hclust(dt.hclust, k = 5)

# 4.
dt.cut = as.factor(cutree(dt.hclust, k = 5))


# Comparaison des deux méthodes
# 1.
table(dt.cut, as.factor(dt.kmeans$cluster))

# 2.
par(mfrow=c(1,2))
plot(dt$Means, col= c(1,2,3,4,5)[as.factor(dt.kmeans$cluster)], main = "Classification k-means")
plot(dt$Means, col= c(2,1,3,4,5)[dt.cut], main = "Classification hclust")
# On remarque une bonne corrélation entre hclust et kmeans excepté pour deux individus (les plus proche de 0)
# qui sont mis dans le groupe inactif par hclust

# Nous voyons donc que nous pouvons séparer nos peptide en differents groupes et donc observer une différence entre les peptide sans activité
# de ceux présentant une activité, donc une liaison au recepteur GPR54. Il est donc maintenant possible d'étudier la différences de séquence 
# entre le groupe basale par rapport aux groupes présentant une activitée

# 3.
summary(dt[which(dt.cut == dt.cut[1]),5:16])
summary(dt[which(dt.cut != dt.cut[1]),5:16])

# Nous pouvons eliminer des positions d'interet: AA4, AA6, AA7, AA8, C.ter car ces position ne présente soi aucune différence entre les peptides 
# liés et ceux non liés soi un profil similaire
# En revanche on remarque que pour les positions AA1, AA5, AA9, AA10 l'acide amine est a chaque fois identique pour tout les peptides lie donc pourrait des
# position conserve entre les peptides necessaire a leur liaison au recepteur et il est difficile de comparer les position N.ter, AA2, AA3 qu'il faudrait 
# visualiser au cas par cas, nous les conservons donc dans notre liste de position potentielle qui est donc composé de N.ter, AA1, AA2, AA3, AA5, AA9, AA10
# Pour voir si ces positions sont vraiment essentielles, nous allons réaliser un apprentissage se basant sur ces positions et prédsant si un peptide appartient 
# à la classe sans activite ou si elle appartient à celle des peptides lies au recepteur.



# ---- Apprentissage sur les donnees ----
library(nnet)

classe = c()
for(i in 1:length(dt.cut)){
  if(dt.cut[i] == dt.cut[1]){
    classe = c(classe, 0)
  } else {
    classe = c(classe, 1)
  }
}
classe
dt = cbind(dt, classe)


dt.1 = dt[,c(5,6,7,8,10,14, 15, 17)]
dt.1$N.term = class.ind(dt.1$N.term)
dt.1$AA1 = class.ind(dt.1$AA1)
dt.1$AA2 = class.ind(dt.1$AA2)
dt.1$AA3 = class.ind(dt.1$AA3)
dt.1$AA5 = class.ind(dt.1$AA5)
dt.1$AA9 = class.ind(dt.1$AA9)
dt.1$AA10 = class.ind(dt.1$AA10)
sortie = dt.1$classe
net.dt = nnet(dt.1[,1:7], sortie, size = 10, maxit = 1000, linout = FALSE)
predict(net.dt, dt.1[,1:7])
# marche pas






