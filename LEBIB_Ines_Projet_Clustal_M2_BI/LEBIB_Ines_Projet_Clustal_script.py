import time
import numpy as np
# import matplotlib.pyplot as plt
# from scipy.cluster.hierarchy import dendrogram, linkage

tps1 = time.time()

'''
La fonction lire_fasta() est une fonction qui prend en entrée un fichier texte
contenant plusieurs fichiers au format FASTA et stocke
les séquences dans un dictionnaire
'''

def lire_fasta (fichier_fasta):
    '''
    Va retourner un dictionnaire de fichiers au format fasta
    
        Parameters :
            fichier_fasta : fichier texte
        
        Returns :
            prot_dict : un dictionnaire de séquences protéiques au format fasta

    '''
    
    compteur_seq = 0
    prot_dict = {}
    with open(fichier_fasta, "r") as fasta_file:
        prot_id = ""
        for line in fasta_file:
            if line.startswith(">"):
                compteur_seq = compteur_seq +1
                prot_id = compteur_seq
                prot_dict[prot_id] = ""
            else:
                prot_dict[prot_id] += line.strip()
                
    return prot_dict


'''
La fonction lire_blosum62 est une fonction qui permet de lire en entrée
un fichier texte contenant la matrice BLOSUM62 et va
stocker les données de la matrice sous forme de dictionnaire
'''

def lire_blosum62(fichier_blosum):
    
    '''
    Va retourner la matrice BLOSUM62 en dictionnaire
    
        Parameters :
            Fichier texte de la matrice BLOSUM

        
        Returns :
            Dictionnaire de la matrice BLOSUM qui contient tous les scores par pair

    '''
    
    with open(fichier_blosum, 'r') as f:
        lignes = f.readlines()

    new_lines = []
    for line in lignes :
        new_lines.append(line)
    lignes = new_lines

    # Récupère les acides aminés de la première ligne
    # en les splitant pour avoir une liste de string
    acides_amines = lignes[0].split()

    # Crée le dictionnaire BLOSUM62
    blosum62 = {}
    for ligne in lignes[1:]:
        valeurs = ligne.split()
        aa = valeurs[0]
        # Dans le dictionnaire blosum62, la clé contient l'acide aminé aa
        blosum62[aa] = {}
        # on itère sur chaque valeur et on ajoute la valeur correspondante
        # des acides aminés aa et acides_amines[i] dans le dictionnaire
        for i, valeur in enumerate(valeurs[1:]):
            blosum62[aa][acides_amines[i]] = int(valeur)

    return blosum62


'''
La fonction score_needleman_wunsch() est une fonction qui permet de
calculer le score d'alignement optimal entre deux séquences

Cette fonction se base sur l'algorithme de Needleman-Wunsch
Elle prend en entrée deux séquences (pour les alignées) ,
une matrice BLOSUM (ici 62) et une pénalité de gap

'''

def score_needleman_wunsch(seq1, seq2, blosum62, gap):

    '''
    Va permettre de retourner un score d'alignement optimal entre deux séquences
    
        Parameters :
            seq 1, seq 2 = string de nos séquences
            blosum62 : un dictionnaire de notre matrice BLOSUM62
            gap = un entier, ici il vaut -5
        
        Returns :
            Un score de l'alignement en int

    '''
    
    # Initialisation des variables m, n et scores
    m, n = len(seq1), len(seq2)
    score = []
    score = [[0] * (n + 1) for i in range(m + 1)]

    # Initialisation de la première ligne et première colonne
    # de la matrice des scores avec les pénalités de gap 
    for i in range(m + 1):
        score[i][0] = gap * i
    for j in range(n + 1):
        score[0][j] = gap * j
        
        
    '''
    On calcule le score pour toutes les positions i et j
    (en excluant la première ligne et la première colonne)
    On calcule pour chaque position (i,j) le score que l'on peut
    avoir suivant les cases précédentes et on prend
    le score maximal des trois scores possibles
    '''
    
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_match = score[i - 1][j - 1] + blosum62[seq1[i - 1]][seq2[j - 1]]
            score_suppr = score[i - 1][j] + gap
            score_inser = score[i][j - 1] + gap
            score[i][j] = max(score_match, score_suppr, score_inser)
            
    # On retourne la dernière valeur de la matrice (dernière case en bas à droite)
    return score[m][n]


'''
La fonction matrice_score() est une fonction qui permet de
calculer la matrice des score entre toutes les séquences
contenues dans le fichier FASTA
'''

def matrice_score(fichier_fasta) :

    '''
    Va afficher la matrice de score pour chaque paire de séquences
    
        Parameters :
            fichier_fasta = dictionnaire des fichiers fasta
            blosum62 = dictionnaire de la matrice BLOSUM62
            gap = valeur de pénalité, un entier, ici il vaut -5                       

        
        Returns :
            score = un tableau de scores de séquences deux à deux


    '''
    
    n = len(fichier_fasta)
    #On initialise la matrice carré des scores à 0
    score = np.zeros((n,n))
    
    # On itère pour chaque sur toutes les paires de séquences fasta
    # en faisant bien attention de ne pas itérer deux fois
    #le même paire (grâce à 'for j in range(i+1, n)')
    for i in range(n) :
        for j in range(i+1, n) :
            # On appelle la fonction score_needleman_wunsch afin
            # de calculer le score d'alignement de chaque paire
            score[i][j] = score_needleman_wunsch(fichier_fasta[i+1], fichier_fasta[j+1], blosum62, gap)
            # On écrit la matrice de score de maniere symétrique
            score[j][i] = score[i][j]
            
    return score

'''
La fonction matrice_distance() permet de convertir la matrice de score en
matrice de distance. Le but est de ne pas avoir de valeurs négatives
dans la matrice

Auparavant, les scores les plus élevées signifiaient deux séquences proches.

Avec la matrice des distances, une distance petite, signifie
des séquences proches.

'''

def matrice_distance(matrice_score) :

    '''
    Va convertir une matrice de score en matrice de distance
    
        Parameters :
            matrice_score = tableau de scores de séquences deux à deux
        
        Returns :
            matrice_dist = tableau de distances de séquences deux à deux

    '''
    
    n = len(matrice_score)
    #On initialise la matrice carré des scores à 0
    matrice_dist = np.zeros((n,n))
    dist_max = matrice_score.max()
    dist_min = np.min(matrice_score)
    
    for i in range(len(matrice_score)) :
        for j in range(i+1, len(matrice_score)) :
            matrice_dist[i, j] = dist_max - (matrice_score[i,j] + dist_min)
            matrice_dist[j, i] = matrice_dist[i, j]
            
    return matrice_dist
    

'''                     
La fonction ConstruireArbre permet de construire un arbre à partir d'une matrice de distances
Il utilise le principe de construction d'un arbre UPGMA
(Unweighted Pair Group Method with Arithmetic Mean) mais de manière incrémentale
'''


def ConstruireArbre(MatriceDistance):
    
    '''
    Va construire un arbre incrémentiel grâce à la méthode UPGMA
    
        Parameters :
            MatriceDistance = tableau de distances de séquences deux à deux
        
        Returns :
            tree = un arbre UPGMA de nos séquences

    '''
    
    n = len(MatriceDistance)
    
    # Initalisation d'une liste de clusters 
    clusters = [[i] for i in range(n)]
    
    distances = MatriceDistance.copy()
    
    while len(clusters) > 1:
        # On cherche les deux clusters les plus proches
        min_distance = float('inf')
        min_indices = (0, 0)
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                distance = distances[i][j]
                if distance < min_distance:
                    min_distance = distance
                    min_indices = (i, j)
        # Les coordonnées des deux clusters les plus proches sont stockés dans x et y 
        x, y = min_indices


        # Les deux clusters les plus proches sont fusionnés en un seul cluster
        new_cluster = clusters[x] + clusters[y]
        clusters.append(new_cluster)

        new_distances = [[0] * len(clusters) for i in range(len(clusters))]

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                if i < len(distances) and j < len(distances[i]):
                    # Calcul la distance entre les clusters
                    # grâce à une moyenne arithmétique
                    dist = (distances[x][i] + distances[y][j]) / 2
                    new_distances[i][j] = dist
                    new_distances[j][i] = dist
        # On met à jour la matrice des distances
        distances = new_distances
        
        # Les clusteurs les plus proches sont supprimés
        del clusters[max(x, y)]
        del clusters[min(x, y)]
        
    tree = clusters[0]
    for i in range(len(tree)) :
        tree[i] = tree[i]+1
    
    return tree

'''
La fonction alignement_needleman_wunsch() permet de réaliser
l'alignement entre deux séquences à l'aide de l'algorithme
de Needleman & Wunsch
'''

def alignement_needleman_wunsch(sequences, arbre, blosum62, gap):
    '''
    Va nous permettre de visualiser les alignements entre deux séquences
    
    Parameters:
        sequence: qui sont des séquences issues des fichiers fasta
        arbre : arbre fréquentiel des séquences
        blosum62: dictionnaire de la matrice BLOSUM62
        gap: valeur de pénalité, un entier, ici il vaut -5
        
    Returns:
        Va nous retourner les alignements de séquences
    '''
    
    # On extrait les deux séquences à alignées
    seq1 = fasta[arbre_UPGMA[0]]
    seq2 = fasta[arbre_UPGMA[1]]         

    m = len(seq1)
    n = len(seq2)
    
    score_alignement = [[0] * (n + 1) for i in range(m + 1)]

    for i in range(m + 1):
        score_alignement[i][0] = gap * i
    for j in range(n + 1):
        score_alignement[0][j] = gap * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_match = score_alignement[i - 1][j - 1] + blosum62[seq1[i - 1]][seq2[j - 1]]
            score_suppr = score_alignement[i - 1][j] + gap
            score_inser = score_alignement[i][j - 1] + gap
            score_alignement[i][j] = max(score_match, score_suppr, score_inser)
            
            
    # Tracement des alignements
    alignment_a, alignment_b = "", ""
    i, j = m, n
    
    while i > 0 and j > 0:
        #On assigne le score de score_alignement suivant la case courante (là où on se situe)
        score_case = score_alignement[i][j]
        score_diag = score_alignement[i - 1][j - 1]
        score_up = score_alignement[i][j - 1]
        score_left = score_alignement[i - 1][j]
        
        ## Ces conditions servent à verifier quel est le score maximal de ces cases ##
        
        # Cas du match de seq1 et seq2 : les deux caractères sont ajoutés dans les
        # deux alignements et on recule en diagonale
        if score_case == score_diag + blosum62[seq1[i - 1]][seq2[j - 1]]:
            alignment_a += seq1[i - 1]
            alignment_b += seq2[j - 1]
            i -= 1
            j -= 1
            
        # Cas de l'insertion/gap dans seq1 :
        # on ajoute le caractère de seq2 dans alignment_b
        # et on insert un gap à alignment_a
        # On passe à la case du haut
        elif score_case == score_up + gap:
            alignment_a += '-'
            alignment_b += seq2[j - 1]
            j -= 1

        # Cas de l'insertion/gap dans seq2 :
        # on ajoute le caractère de seq1 dans alignment_a
        # et on insert un gap à alignment_b
        # On passe à la case de gauche
        
        elif score_case == score_left + gap:
            alignment_a += seq1[i - 1]
            alignment_b += '-'
            i -= 1
            
    # On complète l'alignement jusqu'à ce que i et j atteignent
    # tous les deux 0
    
    while i > 0:
        alignment_a += seq1[i - 1]
        alignment_b += '-'
        i -= 1
    while j > 0:
        alignment_a += '-'
        alignment_b += seq2[j - 1]
        j -= 1

    # Inverse l'ordre pour avoir les bonnes séquences
    alignment_a = alignment_a[::-1]
    alignment_b = alignment_b[::-1]

    return alignment_a, alignment_b
    

if __name__ == '__main__':
    
    gap = -5
    print("La valeur du gap ici vaut : ", gap , "\n")

    fasta = lire_fasta('/LEBIB_Ines_Projet_Clustal_M2_BI/Fichier_fasta_2.txt')
    # print("Les séquences protéiques : " , fasta, "\n")
 
    
    blosum62 = lire_blosum62('/LEBIB_Ines_Projet_Clustal_M2_BI/Blosum62.txt')
    # print("Le dictionnaire de la matrice BLOSUM62 : ", blosum62,  "\n")
    
    score = score_needleman_wunsch(fasta[1], fasta[2], blosum62, gap)
    print("Score d'alignement:", score , "\n")

    matrice_des_scores = matrice_score(fasta)
    print("La matrice des scores est : \n", matrice_des_scores, "\n")
    
    matrice_des_distances = matrice_distance(matrice_des_scores)
    print("La matrice des distances est : \n", matrice_des_distances, "\n")
    
    '''

    # Pour affichier l'arbre : 

    dendro = linkage(matrice_des_distances, method='average')

    # Plot du dendrogramme :
    
    plt.figure(figsize=(20, 20))
    dendrogram(dendro)
    plt.title('Arbre UPGMA')
    plt.xlabel('Seq_ID')
    plt.ylabel('Distance')
    plt.show()
    
    '''
        
    arbre_UPGMA = ConstruireArbre(matrice_des_distances)
    print("L'arbre fréquentiel de nos séquences est : \n ", arbre_UPGMA, "\n")

    alignement_séquences = alignement_needleman_wunsch(fasta, arbre_UPGMA, blosum62, gap)
    print("Alignement des séquences : \n ", alignement_séquences, "\n")

    tps2 = time.time()
    print(f"Le temps d'exécution du code pour 10 séquences de 400 acides aminés en moyenne est : {tps2 - tps1:.3} secondes \n")

