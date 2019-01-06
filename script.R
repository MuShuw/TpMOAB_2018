library(skimr)

# ---- Import des donnees ----
dt = read.csv("Inhibiteurs_GPR54.dat.txt")



# ---- Visualisation des donnees ----
# 1.
dim(dt)
# 79 lignes et 16 colones
summary(dt)
skimr::skim(dt$Means)
# Necessite de transformer la colone N en variable qualitative
dt[,2] = as.factor(dt[,2])

# 2.
par(mfrow = c(1,2))
plot(dt[,3], main = "MEANS")
plot(dt[,4], main = "SD")
dev.off()

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



# ---- Apprentissage sur les donnees ----

# -- Creation d'une sequence consensus --
# 4.1.1 Creation d'une liste contenant les tables de frequence à chaque position
Diff_element = apply(dt[,5:16], 2, table)

# Creation d'un vecteur null
Consensus = vector()

# On va entrer à chaque position de la sequence l'element majoritaire dans le jeu

# 4.1.2 On va parcourir notre liste de table
for ( i in seq_along(Diff_element)){
  # Chaque table est enregistrée dans un df temporaire
  tmp = as.data.frame(Diff_element[[i]])
  # Les elements de la première colonne sont transformer en character pour
  # tester le matching
  tmp$Var1 <- as.character(tmp$Var1)
  # Les elements de la deuxième colonne sot convertie en valeur numérique
  tmp$Freq <- as.numeric(tmp$Freq)
  # Chaque ieme element de consensus va contenir l'element dont la freq est la
  # plus importante, ie l'element majoritaire
  Consensus[i] = tmp[which(tmp$Freq == max(tmp$Freq)),1]
}


# ---- Creation d'une table contenant les séquences sous forme de character ----
# 4.1.1 Creation d'une table de sequence
table_de_seq = apply(dt[,5:16], 2, as.character)


# ---- 4.1.3 Creation d'une matrice de similarité pour l'apprentissage ----
# Creation d'une matrice de 0, avec une nombre de lignes egale au nombre de libre du jeu de donnees
# et de 12 colonnes pour les 12 elements de la séquence
matching_matrix = as.data.frame(matrix(0, ncol = 12, nrow = dim(dt)[1]))
# On va nommer les colonnes de cette matrice
colnames(matching_matrix) = colnames(dt[,5:16])

# On va parcourir les cases [i,j] de la matrice et mettre 1 dans celle-ci si le j element de la sequene
# de la ieme ligne du jeu different du consensus
# On parcours les colonnes de la matrice des differences
for ( j in 1:12 ){
  # Puis ces lignes
  for ( i in 1:dim(matching_matrix)[1]){
    # Si l'element de le sequence du jeu de donnée diffère du consensus
    if ( table_de_seq[i,j] != Consensus[j]){
      # On lui attribue la valeur 1
      matching_matrix[i,j] = 1
    }
  }
}


# 4.1.4 ---- Creation d'une matrice frequence proche d'une pssm ----

# Creation d'une matrice de frequence contenant des zero
freq_mat = as.data.frame(matrix(0, ncol = 12, nrow = dim(dt)[1]))
# On va nommer les colonnes de cette matrice
colnames(matching_matrix) = colnames(dt[,5:16])


# On va parcourir les cases [i,j] de la matrice et mettre dans celle-ci la frequence de l'element
# en position [i,j] relativement à la colonne j ( les frequence des element sont relatives aux elements
# de leur colonnes, donc la frequence en jème position dans la sequence )

# On parcours les colonnes de la matrice des frequnces
for ( j in 1:12 ){
  # Puis ses lignes
  for ( i in 1:dim(freq_mat)[1]){
    # Chaque table est enregistrée dans un df temporaire
    tmp = as.data.frame(Diff_element[[j]])
    # Les elements de la première colonne sont transformer en character pour
    # retrouver sa frequence
    tmp$Var1 <- as.character(tmp$Var1)
    # Les elements de la deuxième colonne sot convertie en valeur numérique puis divisé
    # par le nombre de lignes de la table de sequence
    tmp$Freq <- as.numeric(tmp$Freq)/dim(dt)[1]
    # On enregistre la frequence de l'element en position i,j dans la matrice de frequence
    # en position i,j
    freq_mat[i,j] = tmp[which(tmp$Var1 == table_de_seq[i,j]),2]
  }
}

# 4.2 ---- Creation de groupe pour la classification ----
# - Groupe dichotimique selon dt$Means -
# Un valeur seuil est enregistrée ( ici la médiane a été choisie, tester avec différente valeurs, ex :
# different quartiles, appartenant au cluster1 / exterieur au cluster1...[ Ici on changerai le seuil
# et la colonne a tester [ avant ajouter une colonne avec le cluster ] ])
seuil = 0.2867

# Le groupe 0 contiendra les individus dont la valeur d'activité est inférieur à la medianne
groupe0 = rep(0,dim(dt)[1])
groupe0[which(dt$Means <= seuil)] = 1
# Le groupe 1 contiendra les individus dont la valeur d'activité est strictement supérieur à la medianne
groupe1 = rep(0,dim(dt)[1])
groupe1[which(dt$Means > seuil)] = 1
# On joint ces deux groupes
groupe = cbind(groupe0,groupe1)

# 4.3.1 - Utilisation de nnet -
library(nnet)

# Observation sans split dans un premier temps$
# Cette fonction sera fournie. Elle permet de calculer les performances pour une classification
# en lui donnant une table de confusion
performanceVal <- function(ConfTable){
        # Specificité
        Spe = ConfTable[2,2]/sum(ConfTable[,2])
        # Sensibilité
        Sen = ConfTable[1,1]/sum(ConfTable[,1])
        # Taux de succes
        Acc = (ConfTable[1,1]+ConfTable[2,2]) /sum(ConfTable)
        # Taux d'erreur
        Err = (ConfTable[1,2]+ConfTable[2,1]) /sum(ConfTable)
        MyPerf = as.data.frame(cbind(Spe, Sen, Acc, Err))
        return(MyPerf)
}

# 4.3.1 -- Spitting des données --
# Vous pouvez jouer avec la seed afin d'obtenir different sampling, pour ces sampling
# vous constaterez de meilleures performances pour des sampling où le training contient une bonne
# partie des valeurs qui n'appartenaient pas au cluster 1
seed = runif(1,1,10000)
set.seed(191.0967) # 43, 44, 300, 191.0967, 114.3176
seed
Index = sample(79, 30)
Freq_train = freq_mat[-Index,]
Freq_test = freq_mat[Index,]
Match_train = matching_matrix[-Index,]
Match_test = matching_matrix[Index,]
groupe_train = groupe[-Index,]
groupe_test = groupe[Index,]

# 4.3-- training avec nnet --
# ---- Evaluation avec la Matching Matrix ----
# 43, 44, 300, 191.0967, 114.3176
net.dt_ecart = nnet(Match_train, groupe_train, size = 50, maxit = 1000, softmax  =  T)

# Mesure des perfs sur le train
x = predict(net.dt_ecart, Match_train)
ConfTab = table( round(x)[,2], groupe_train[,2])
PerfMatch_train = performanceVal(ConfTable = ConfTab)
ConfTab
PerfMatch_train

# Mesure des perfs sur le test
x = predict(net.dt_ecart, Match_test)
ConfTab = table( round(x)[,2], groupe_test[,2])
PerfMatch_test = performanceVal(ConfTable = ConfTab)
ConfTab
PerfMatch_test

# ---- Evaluation avec la Frequence  Matrix ----
# seed 300
net.dt_freq = nnet(Freq_train, groupe_train, size = 50, maxit = 1000, softmax  =  T)
# Mesure des perfs sur le train
x = predict(net.dt_freq, Freq_train)
ConfTab = table( round(x)[,2], groupe_train[,2])
PerfFreq_train = performanceVal(ConfTable = ConfTab)
PerfFreq_train
# Mesure des perfs sur le test
x = predict(net.dt_freq, Freq_test)
ConfTab = table( round(x)[,2], groupe_test[,2])
PerfMatch_test = performanceVal(ConfTable = ConfTab)
PerfMatch_test


# 4.4 ---- Element critique de la sequence ----


# 4.4. Evaluation avec la Matching Matrix
net.dt_ecart = nnet(matching_matrix, groupe, size = 50, softmax  =  T)
pred_full_match = predict(net.dt_ecart, matching_matrix)
pred_full_match = round(pred_full_match)
ConfTab = table( pred_full_match[,2], groupe[,2])
PerfMatch = performanceVal(ConfTable = ConfTab)
ConfTab
PerfMatch

# 4.4 Evaluation avec la Frequence  Matrix
net.dt_freq = nnet(freq_mat, groupe, size = 50, softmax  =  T)
pred_full_freq = predict(net.dt_freq, freq_mat)
pred_full_freq = round(pred_full_freq)
ConfTab = table(pred_full_freq[,2], groupe[,2])
PerfFreq = performanceVal(ConfTable = ConfTab)
ConfTab
PerfFreq

# On constate que l'utilisation d'une matrice de frequence permet une meilleur distinction,
# Elle est porteuse de plus d'information


# On va selectionner les sequences supérieur au seuil pour nos deux modèles complets
library(dplyr)
# On va d'abord créer des dataframe contenant les sequence predite comme sup/inf au seuil et
# les sequences réels
sequence_match_sup = dt[which(pred_full_match[,2]==1),5:16]
sequence_match_inf = dt[which(pred_full_match[,1]==1),5:16]
sequence_freq_sup = dt[which(pred_full_freq[,2]==1),5:16]
sequence_freq_inf = dt[which(pred_full_freq[,1]==1),5:16]
sequence_true_sup = dt[which(groupe1==1),5:16]
sequence_true_inf = dt[which(groupe0==1),5:16]

# Puis nous allons déterminer les intersections de ces dataframe
common_true_match_sup = inner_join(sequence_match_sup,sequence_true_sup)
common_true_freq_sup = inner_join(sequence_freq_sup,sequence_true_sup)
common_true_match_inf = inner_join(sequence_match_inf,sequence_true_inf)
common_true_freq_inf = inner_join(sequence_freq_inf,sequence_true_inf)

# On fait des resumé de ces intersections
# L'intersection est utilisé pour méler les information brut réel avec les interaction
# perçu par le modèle

skim(common_true_match_sup)
skim(common_true_match_inf)

# Premièrement nous savons que AA6 et C-ter sont des invariants

skim(common_true_match_sup)
skim(common_true_match_inf)
summary(common_true_match_sup)
summary(common_true_match_inf)
# On constate qu'une lysine acétylé en position AA4 est carachtéristique d'une valeur supérieur au seuil :
dt[which(dt$AA4 == "Lys(Ac)"),]
# L'individu 20 semble bien être une erreur de manipulation
# On constate la même chose en position AA5 tandis que Thr semble être exclusivement chez les inférieur
dt[which(dt$AA5 == "Lys(Ac)"),]
# Ces 4 individus observés semble montrer que un Acetyle en position N-ter diminue l'activitée
dt[which(dt$AA5 == "Thr"),]

# Pour les positions AA8 et AA9 il n'apparait pas de disticntions particulières cependent les fréquences
# de AA8:Y[tz]Leu et AA9:Arg(Me) sont moins importantes dans la catégorie supérieur

# En position AA1 D-Tyr par sa fréauence chez les inférieurs et les observations précédante
# semble réduire l'activité
dt[which(dt$AA1 == "D-Tyr"),]

# Les gaps en position AA2 semblent avoir un effet reduceur sur l'activitée
dt[which(dt$AA2 == "_"),]



# ---- Conclusion ----
# De ce TP vous devez retenir les deux méthodes de clustering présentées
# et l'informations quelles peuvent apporter, la façon de créer les splits et
# de les utiliser pour la création d'un modèle nnet pour de la classification,
# l'importance de la répartition des données dans ces deux sous-ensembles et
# l'importance du niveau informatif d'un jeu de donnée.
# Pour des raisons pratiques le seuil est arbitraire,
# un choix plus rigoureux passe par une compréhension biologique des valeurs de l'activitée,
# de même les matrices créées donne des performances relativement mauvaises pour de la prédiction.
# Une meilleur approche aurait était de récuperer des informations sur
# les caractèristiques physicochimique des éléments des séquences.

