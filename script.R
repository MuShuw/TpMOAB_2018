library(skimr)
setwd("~/Bureau/projet_moab/")

# ---- Nettoyage ----
# Suppression pour fournir le jeu de donnée du TP, retire les colones inutiles, met de NA quand une valeur est manquante et
# suppression de la ligne 64 qui contient un peptide sans moyenne donc inutilisable
# suppression de la ligne 70 qui contient un peptide sans sequence donc inutilisable
dt = read.csv("Data/Inhibiteurs_GPR54.dat (copie).txt", sep="\t", dec = ".")
dt = dt[,-c(19,18,4)]
dt = dt[-c(70,64),]
dt[c(78, 51, 35, 5), 4] = NA
write.csv(dt, file = "Inhibiteurs_GPR54.dat(copie-Bertrand).txt", row.names = FALSE)



# ---- Import des donnees ----
dt = read.csv("Inhibiteurs_GPR54.dat(copie-Bertrand).txt")



# ---- Visualisation des donnees ----
# 1.
dim(dt)
# 80 lignes et 16 colones Elyas : j'en ai 79 ??
summary(dt)
skimr::skim(dt$Means)
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
#hist(dt[which(dt.cut == dt.cut[1] & dt[,3] < 10 ),3], breaks = 50)
summary(dt[which(dt.cut != dt.cut[1]),5:16])

# Nous pouvons eliminer des positions d'interet: AA4, AA6, AA7, AA8, C.ter car ces position ne présente soi aucune différence entre les peptides
# liés et ceux non liés soi un profil similaire
# En revanche on remarque que pour les positions AA1, AA5, AA9, AA10 l'acide amine est a chaque fois identique pour tout les peptides lie donc pourrait des
# position conserve entre les peptides necessaire a leur liaison au recepteur et il est difficile de comparer les position N.ter, AA2, AA3 qu'il faudrait
# visualiser au cas par cas, nous les conservons donc dans notre liste de position potentielle qui est donc composé de N.ter, AA1, AA2, AA3, AA5, AA9, AA10
# Pour voir si ces positions sont vraiment essentielles, nous allons réaliser un apprentissage se basant sur ces positions et prédsant si un peptide appartient
# à la classe sans activite ou si elle appartient à celle des peptides lies au recepteur.

############ bazar d'observation ###############
dt_temp = apply(dt[,5:16], MARGIN = 2, as.factor)
skim(as.data.frame(dt_temp))
dt1 = as.data.frame(cbind(dt[,-c(5:16)], dt_temp))
dt1 = dt1[which(dt.cut != dt.cut[1]),]
mean_Nter = split(dt1$Means, dt1$N.term)
dev.off()
boxplot(mean_Nter, notch = FALSE, varwidth = TRUE)
mean_AA2 = split(dt1$Means, dt1$AA2)
dev.off()
boxplot(mean_AA2, notch = FALSE, varwidth = TRUE)
mean_AA3 = split(dt1$Means, dt1$AA3)
dev.off()
boxplot(mean_AA3, notch = FALSE, varwidth = TRUE)
###################################################
# ---- Apprentissage sur les donnees ----

# ---- Creation d'une sequence consensu ----

Diff_element = apply(dt[,5:16], 2, table)
Consensus = vector()
for ( i in seq_along(Diff_element)){
        tmp = as.data.frame(Diff_element[[i]])
        tmp$Var1 <- as.character(tmp$Var1)
        tmp$Freq <- as.numeric(tmp$Freq)
        print(tmp$Freq/79)
        Consensus[i] = tmp[which(tmp$Freq == max(tmp$Freq)),1]
}


# ---- Creation d'une table d'ecart avec matching matrix pour nnet ----

dt_ecart = dt

seuil = 0.2867
groupe1 = rep(0,dim(dt)[1])
groupe1[which(dt$Means > seuil)] = 1
groupe0 = rep(0,dim(dt)[1])
groupe0[which(dt$Means <= seuil)] = 1
groupe = cbind(groupe0,groupe1)
matching_matrix = as.data.frame(matrix(0, ncol = 12, nrow = dim(dt)[1]))
colnames(matching_matrix) = colnames(dt[,5:16])
for ( i in 1:dim(matching_matrix)[1]){
        for ( j in 1:12 ){
                if ( dt[i,j+4] != Consensus[j]){
                        matching_matrix[i,j] = 1
                }
        }

}

# apply(dt_ecart[,6:17], 2, var)
# dt_ecart$AA6 <- NULL
# dt_ecart$C.ter <- NULL

library(nnet)

net.dt_ecart = nnet(matching_matrix, groupe, size = 50, maxit = 1000, softmax  =  T)
x = predict(net.dt_ecart, matching_matrix)
table(round(x)[,1], groupe[,1])
table(round(x)[,2], groupe[,2])
# marche pas

#### Creation ----

# Diff_element_freq = lapply(Diff_element, function(x) sapply(x, function(x) x/79))
freq_mat = as.data.frame(matrix(0, ncol = 12, nrow = dim(dt)[1]))
colnames(matching_matrix) = colnames(dt[,5:16])
table_de_seq = apply(dt[,5:16], 2, as.character)


for ( i in 1:dim(dt_freq)[1]){
        for ( j in 1:12 ){
                tmp = as.data.frame(Diff_element[[j]])
                tmp$Var1 <- as.character(tmp$Var1)
                tmp$Freq <- as.numeric(tmp$Freq)/79

                freq_mat[i,j] = tmp[which(tmp$Var1 == table_de_seq[i,j]),2]

        }

}

Index = sample(79, 30)
Mat_train = freq_mat[-Index,]
Mat_test = freq_mat[Index,]
net.dt_freq = nnet(Mat_train, groupe[-Index,], size = 30, maxit = 5000, softmax = T)
x = predict(net.dt_freq, Mat_train)
round(x)
table(round(x)[,2], groupe[-Index,2])
