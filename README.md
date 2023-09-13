# LEBIB_Ines_Projet_Clustal_M2BI


LEBIB Inès - M2 BI                                      

Projet court Clustal

13/09/2023




### **Alignement multiple heuristique par la méthode Clustal** ###


L'objectif de ce projet est de réaliser un script Python qui reprend une méthode décrite dans un article écrit par Desmond G. Higgins et Paul M. Sharp en Avril 1989 *Fast and sensitive multiple sequence alignments on a microcomputer*. Les auteurs présentent une méthodologie d'alignement de séquences multiple heuristique par la méthode Clustal. 
Durant ce projet, nous allons nous intéresser exclusivement aux séquences protéiques. 



#### **Introduction** ####

De nombreuses stratégies ont été mis en place afin de faire un alignement multiple. Nous allons nous nous intéresser à la méthode de programmation dynamique de Needleman & Wunsch (1970) et la méthode UPGMA (Sneath et Sokal, 1973) pour construire un arbre enraciné. 



#### **Informations importantes** ####

Le répertoire dans lequel est stocké le projet contient : 

- L'article scientifique écrit par Desmond G. Higgins et Paul M. Sharp sous format pdf : *Article_Clustal.pdf*
- Les séquences protéiques sous format txt : *fichier_fasta.txt*
- La matrice BLOSUM62 sous format de txt : *Blosum62.txt*
- Le script python : *LEBIB_Ines_Projet_Clustal_script.py*
- Le JupyterLab : *LEBIB_Ines_Projet_Clustal_Jupyter_Lab.ipynb*
- Le rapport du projet sous format pdf : *LEBIB_Ines_Rapport_Clustal_Projet.pdf*




#### **Packages** ####

Durant ce projet, quelques packages python sont à installer :

```
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
```



#### **Les pré-requis** ####

Afin de pouvoir utiliser au mieux le code, certaines fonctions vont être utiles.

Pour le déroulement optimal de notre script, nous aurons besoin d'avoir des séquences protéiques au format FASTA sous forme de dictionnaire. 

A partir d'un fichier texte contenant toutes les séquences FASTA à aligner, indiquez le chemin de votre fichier et exécuter cette ligne de commande afin de créer le dicionnaire :

```
fasta = lire_fasta('path/Fichier_fasta_2.txt')
print("Les séquences protéiques : " , fasta, "\n")

```

Plusieurs fonctions vont également prendre en entrée une matrice BLOSUM62, également sous forme de dictionnaire. La matrice BLOSUM62 est une matrice de similarité/substitution. Elle va permettre de définir un score de similarité ou de ressemblance entre deux acides aminés. L'intérêt de cette fonction est d'avoir un dictionnaire les valeurs de substitutions de tous les couples d'acide aminé. Cette fois encore, il faut indiquer le chemin de votre matrice BLOSUM. 

La commande pour exécuter cette fonction est : 

```
blosum62 = lire_blosum62('path/Blosum62.txt')
print("Le dictionnaire de la matrice BLOSUM62 : ", blosum62,  "\n")

```

Le dernier pré-requis est la valeur du gap. Le gap est une valeur de pénalité fixe pour chaque écart, ici la valeur du gap est -5.

```
gap = -5
print("La valeur du gap ici vaut : ", gap , "\n")
```



#### **Guide d'utilisation** ####

Maintenant que les packages et les entrées de nos fonctions sont prêtes, il est temps de rentrer dans le coeur du projet.


La première grande étape de ce projet est de calculer une matrice de score de similarité entre toutes les paires de séquences grâce à l’algorithme de Needleman & Wunsch.

Cet algorithme prend en argument deux séquences protéiques, la matrice BLOSUM62, et une pénalité de gap. En sortie, il va affichier un score de similarité entre deux séquences protéiques.

La fonction va calculer les scores de toutes les positions (i,j) d'acides aminés de nos séquences et va retourner le score de l'alignement optimal entre les deux séquences.

```
score = score_needleman_wunsch(fasta[1], fasta[2], blosum62, gap)
print("Score d'alignement:", score , "\n")
```

A présent, le but est d'exectuer une fonction qui permet de calculer une matrice de score de similarité entre toutes les paires de séquences contenues dans le fichier FASTA. 

La fonction prend en entrée le fichier FASTA et retourne la matrice de score.

```
matrice_des_scores = matrice_score(fasta)
print("La matrice des scores est : \n", matrice_des_scores, "\n")
```

Ce qui nous intéresse maintenant est de convertir la matrice des scores en matrice de distances. Cela permettra de ne pas avoir de valeurs négatives dans la matrice. 

Dans la matrice des scores, des valeurs élevées signifiaient que deux séquences étaient proches ; avec la matrice des distances, des valeurs petites indiquent des séquences proches.

La fonction prend en entrée la matrice des scores et retourne une matrice de distance.

```
matrice_des_distances = matrice_distance(matrice_des_scores)
print("La matrice des distances est : \n", matrice_des_distances, "\n")
```

Il est possible d'affichier l'arbre à partir de la matrice des distances, grâce à des packages installés plus haut : 

```
dendro = linkage(matrice_des_distances, method='average')
plt.figure(figsize=(20, 20))
dendrogram(dendro)
plt.title('Arbre UPGMA')
plt.xlabel('Seq_ID')
plt.ylabel('Distance')
plt.show()
```

La fonction qui va suivre va permettre de construire un arbre de manière incrémentiel en s'inspirant de la méthode UPGMA (Unweighted Pair Group Method with Arithmetic Mean). UPGMA est issue des méthodes de clusterisation. L’objectif est, par itérations successives, de réduire progressivement la taille de la matrice de distances en la transformant en arbre enraciné, jusqu’à ce qu’il ne reste qu’un seul cluster.

Cette fonction prend en entrée la matrice de distance et retourne l'arbre 

```
arbre_UPGMA = ConstruireArbre(matrice_des_distances)
print("L'arbre fréquentiel de nos séquences est : \n ", arbre_UPGMA, "\n")
```

La dernière étape est de faire un alignement multiple de séquences. Pour cela, il faut utiliser l'algorithme de Needleman & Wunsch. Il va permettre de construire les alignements entre les différentes séquences protéiques.
Malheureusement je n'ai pas réussis à faire cela. Mais j'ai tout de même réussis à faire un alignement entre deux séquences protéiques.
Le but est de retracer le ‘chemin’ suivant le score maximal des alignements qui donnera l’alignement optimal entre les deux séquences. 
Cette fonction prend en entrée des séquences (suivant l’ordre de ramification de l’arbre UPGMA), l'arbre UPGMA, la matrice BLOSUM62 et le gap. 

```
alignement_séquences = alignement_needleman_wunsch(fasta, arbre_UPGMA, blosum62, gap)
print("Alignement des séquences : \n ", alignement_séquences, "\n")
```
                 









