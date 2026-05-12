import numpy as np
import math

def compter_motifs(motifs, nuc):
    """
    Construit la matrice de comptage à partir d'une liste de motifs.
    entrée motifs : liste de chaînes de même longueur k
    entrée nuc    : tuple de nucléotides, ex. ('A', 'C', 'G', 'U')
    sortie M      : tableau numpy de forme (len(nuc), k), dtype int
    """
    M = np.zeros((len(nuc), len(motifs[0])), dtype=int)  # M[a, i] = nombre de motifs avec le nucléotide a à la position i
    for i in range(len(motifs[0])):
        for j in range(len(motifs)):
            M[nuc.index(motifs[j][i])][i] += 1
    
    return M

def construire_pwm(M, pseudocount=1):
    """
    Calcule la PWM à partir d'une matrice de comptage.
    entrée M           : tableau numpy de forme (len(nuc), k) — matrice de comptage
    entrée pseudocount : entier ajouté à chaque case avant normalisation (par défaut 1)
    sortie PWM         : tableau numpy de même forme, colonnes de somme 1
    """
    PWM = np.zeros((len(M), len(M[0])))
    for i in range(len(M)):
        tot = np.sum(M[:, i]+pseudocount)
        for j in range(len(M[0])):
            PWM[i, j] = (M[i, j]+1)/tot

    return PWM

def logvraisemblance(seq, PWM, bg, nuc):
    """
    Calcule le log-rapport de vraisemblance d'un k-mer par rapport à la PWM.
    entrée seq : chaîne de longueur k
    entrée PWM : tableau numpy de forme (len(nuc), k)
    entrée bg  : dictionnaire {nucléotide: probabilité de fond}
    entrée nuc : tuple de nucléotides — sert à convertir un caractère en indice de ligne de PWM
    sortie ll  : log2-rapport de vraisemblance (float)
    """
    ll = 0.0

    p0_seq = 1
    p_seq = 1
    for i in range(len(seq)):
        p0_seq *= bg[seq[i]]
        p_seq *= PWM[nuc.index(seq[i])][i]

    ll = np.log2(p_seq/p0_seq)
    
    return ll


def rechercherPWM(sequences, k, PWM, bg, nuc):
    """
    Pour chaque séquence, trouve le k-mer qui maximise le log-rapport de vraisemblance.
    entrée sequences : liste de séquences ARN
    entrée k         : taille du motif
    entrée PWM       : matrice de poids position, forme (len(nuc), k)
    entrée bg        : distribution de fond
    entrée nuc       : tuple de nucléotides
    sortie motifsT   : liste de k-mers trouvés (un par séquence)
    """
    motifsT = []
    for seq in sequences:
        best = (None, -99999999999)
        for i in range(0,len(seq)-k+1):
            llkmer = logvraisemblance(seq[i:i+k], PWM, bg, nuc)
            if llkmer > best[1]:
                best = (seq[i:i+k],llkmer)

        motifsT.append(best[0])

    return motifsT

def performance(foundMotifs, trueMotifs):
    """
    Calcule la fraction de motifs correctement prédits.
    entrée foundMotifs : liste de motifs prédits
    entrée trueMotifs  : liste de vrais motifs
    sortie score       : fraction de motifs corrects (float entre 0 et 1)
    >>>performance(['CGAUGC', 'CGAUGC'], ['CGAUGC', 'UUAUGC'])
    0.5
    """
    score = 0.0
    correct = 0
    for i in range(len(trueMotifs)):
        if trueMotifs[i] == foundMotifs[i]:
            correct += 1

    score = correct/len(trueMotifs)
    
    return score

def hamDistancePonderee(seq1, seq2, weights):
    """
    Calcule la distance de Hamming pondérée entre deux séquences de même longueur.
    Les positions identiques contribuent 0 ; les substitutions contribuent weights[a+b].
    entrée seq1    : chaîne de caractères
    entrée seq2    : chaîne de caractères (même longueur que seq1)
    entrée weights : dictionnaire {paire de nucléotides: coût}, ex. {'AG': 0.25, 'AU': 1.0, ...}
    sortie dist    : distance pondérée (float)
    >>>hamDistancePonderee('ACGAUGC', 'GCGAUGC', weights)
    0.25
    >>>hamDistancePonderee('ACGAUGC', 'ACUAUGC', weights)
    1.0
    """
    dist = 0.0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += weights[seq1[i]+seq2[i]]

    return dist

def distanceTotale(motif, sequences, weights):
    """
    Calcule la somme des distances pondérées minimales entre un motif et chaque séquence.
    entrée motif     : chaîne de longueur k
    entrée sequences : liste de séquences
    entrée weights   : dictionnaire de coûts de substitution
    sortie td        : somme des distances minimales (float)
    >>>distanceTotale('ACGAUGC', ['GCGAUGC', 'ACUAUGC'], weights)
    1.25
    """
    td = 0.0
    for seq in sequences:
        k = len(motif)
        L = len(seq)
        minham = 99999999999999999
        for i in range(L-k+1):
            ham = hamDistancePonderee(motif, seq[i:i+k], weights)
            if ham < minham:
                minham = ham

        td += minham
    
    return td

def medianStringSearch(allkmers, sequences, weights):
    """
    Trouve le k-mer qui minimise la distance totale pondérée sur l'ensemble des séquences.
    entrée allkmers  : liste de tous les k-mers à tester
    entrée sequences : liste de séquences
    entrée weights   : dictionnaire de coûts de substitution
    sortie bestMotif : k-mer minimisant la distance totale
    """
    bestMotif = ""
    bestDist = 999999999999999
    for kmer in allkmers:
        dist = distanceTotale(kmer, sequences, weights)
        if dist < bestDist:
            bestDist = dist
            bestMotif = kmer

    return bestMotif

def extraireMotifs(motif, sequences, weights):
    """
    Pour chaque séquence, extrait le k-mer le plus proche du motif consensus.
    entrée motif     : chaîne de longueur k — motif consensus
    entrée sequences : liste de séquences
    entrée weights   : dictionnaire de coûts de substitution
    sortie listMotifs: liste de k-mers (un par séquence)
    """
    listMotifs = []
    for seq in sequences:
        k = len(motif)
        best = (None, 99999999999)
        for i in range(len(seq)-k+1):
            ham = hamDistancePonderee(motif, seq[i:i+k], weights)
            if ham < best[1]:
                best = (seq[i:i+k], ham)
                
        listMotifs.append(best[0])
    
    return listMotifs

def conclusion(perfPWM, perfMS):
    """
    Affiche lequel des deux algorithmes a obtenu la meilleure performance
    entrée perfPWM : score de performance de la PWM (float entre 0 et 1)
    entrée perfMS  : score de performance du Median String (float entre 0 et 1)
    """
    print("Conclusion :")
    print(f"Performance PWM : {perfPWM}\nPerformance MS: {perfMS}")
    algo = "PWM" if perfPWM > perfMS else "Median String"
    print(f"L'algorithme {algo} a obtenu la meilleure performance")
    return