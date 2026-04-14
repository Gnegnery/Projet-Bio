import random
import numpy as np
from tqdm import tqdm
from itertools import product
def generateRandomSequences(n:int, t:int, upper = True):
    """
    Génére plusieurs séquences nucléotidiques aléatoires 
    entrée n : nombre de séquences
    entrée t : taille des séquences
    entrée upper : bool, si True, les nucléotides seront majuscule, False, minuscule
    sortie séquences : liste de séquences nucléotidiques aléatoires ou un string si n = 1
    """
    nuc = ('A', 'C', 'G', 'T')
    sequences = []
    for _ in range(n):
        seq = ""
        for _ in range(t):
            seq += random.choice(nuc)
        if upper == False:
            seq = seq.lower()
            
        sequences.append(seq)
    
    
    return sequences[0] if n == 1 else sequences

def implantMotifs(motif:str, sequences:list, f = 0.9):
    """
    Insére un motif à des positions aléatoires des séquences
    entrée motif : motif qui va être implanté dans les séquences
    entrée séquences : liste de séquences
    entrée f : fraction des séquences qui contiendront le motif
    sortie modified_sequences: liste de séquences ayant le motif implanté
    """
    def insertMotif(motif:str, seq:str):
        tailleSeq = len(seq)
        index = random.randint(0,tailleSeq-1)
        return seq[0:index] + motif + seq[index:]

    modified_sequences = []
    
    for seq in sequences:
        if random.random()<f:
            new_seq = insertMotif(motif, seq)
        else :
            new_seq = seq
        
        modified_sequences.append(new_seq)
                
    return modified_sequences

def searchMotifs(k:int, sequences:list):
    """
    Cherche les motifs de taille k dans un ensemble de séquences
    entrée k : taille du motif
    entrée séquences : liste de séquences
    sortie motifs: dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    >>>searchMotifs(3, ['TAAGTAA', 'TATAA', 'CTATC'])
    {'TAA': 3, 'AAG': 1, 'AGT': 1, 'GTA': 1, 'TAT': 2, 'ATA': 1, 'CTA': 1, 'ATC': 1}
    """

    motifs  = {}
    for seq in sequences:
        for i in range(0, len(seq)-k+1):
            motif_actuel = seq[i:i+k]
            if motif_actuel in motifs:
                motifs[motif_actuel] += 1
            else:
                motifs[motif_actuel] = 1
    
    return motifs

def getTopMotifs(motifs:dict, top, decroissant = True):
    """
    Renvoie les top motifs le plus fréquent
    entrée motifs: dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    entrée top : les top plus fréquent, si top == "all", renvoie le dictionnaire complet trié sans filtrer
    entrée decroissant : bool, si True le dictionnaire est trié par ordre décroissant de valeur
    sortie motifsfreq: dictionnaire contenant les top motifs les plus fréquents, clé = motif, valeur = fréquence d'observation
    >>>getTopMotifs({'TAA': 3, 'AAG': 1, 'AGT': 1, 'GTA': 1, 'TAT': 2, 'ATA': 1, 'CTA': 1, 'ATC': 1}, 2)
    {'TAA': 3, 'TAT': 2}
    """
    
    motifsfreq = {}
    cpt = 0
    for k, v in sorted(motifs.items(), key=lambda item: item[1], reverse=decroissant):
        if top != 'all' and cpt == top:
            break
        motifsfreq[k] = v
        cpt += 1
            
    return motifsfreq
    
def reversecompl(seq:str):
    """Renvoie le brin complémentaire d’une séquence.
    entrée seq : sequence de nucléotides (brin sens)
    sortie     : sequence de nucléotides (brin complementaire)
    >>> reversecompl('AACGTGGCA')
    'TGCCACGTT'
    """
    compl = {'A': 'T', 'C': 'G', 'G': 'C', 'T':'A'}
    rev = ""
    for i in seq[::-1]:
        rev += compl[i]
    return rev

def removeLowComplexeHomo(motifs:list, m:int):
    """
    Enlève les motifs peu complexe ayant m fois le même nucléotide
    entrée motifs: liste de motifs
    entrée m: nombre de nucleotides répétés dans un motif peu complexe
    sortie motifsClean: liste de motifs sans les motifs peu complexe
    >>>removeLowComplexeHomo(['TAA', 'AAG', 'AGT', 'GTA', 'TAT', 'ATA', 'CTA', 'ATC'], 2)
    ['AGT', 'GTA', 'CTA', 'ATC']
    """

    motifsClean = []
    for mot in motifs:
        freqLettre = {}
        most = 0
        for c in mot:
            if c in freqLettre:
                freqLettre[c] += 1
            else :
                freqLettre[c] = 1

            if freqLettre[c] > most:
                most = freqLettre[c]

        if most < m:
            motifsClean.append(mot)
            
    return motifsClean

def removeLowComplexeHetero(motifs:list, n:int, variation = "yes"):
    """
    Enlève les motifs peu complexe ayant n fois un dinucléotide
    entrée motifs: liste de motifs, clé = motif, valeur = fréquence d'observation
    entrée n: nombre di-nucleotides répétés dans un motif peu complexe
    entrée variation : string, si "yes" permettre variation d'un nucléotide
    sortie motifsClean: liste de motifs sans les motifs peu complexe
    >>>removeLowComplexeHetero(['GGTTTGG', 'TGAGTTA', 'TGCCGTG', 'AGAGAGA', 'TCACCGA', 'TTGGTAT', 'AGGGTGG', 'TGGCTTA', 'AGAGTAG', 'GCCCCTC'], 3, variation = "yes")
    ['GGTTTGG', 'TGAGTTA', 'TGCCGTG', 'TCACCGA', 'TTGGTAT', 'TGGCTTA']
    
    """
    motifsClean = []
    for mot in motifs:
        freqNucs = {}
        most = 0
        d_nuc = 2
        t_mot = len(mot)
        lastNucs = ""
        if variation == "yes":
            for i in range(0, t_mot - d_nuc + 1):
                nucs = mot[i:i+d_nuc]
                if nucs in freqNucs:
                    freqNucs[nucs] += 1
                    
                else :
                    freqNucs[nucs] = 1
    
                if freqNucs[nucs] > most:
                    most = freqNucs[nucs]

            if most < n:
                motifsClean.append(mot)

        else:
            for shift in range(d_nuc):
                most = 0
                stop = False
                for i in range(shift, t_mot - d_nuc + 1, 2):
                    nucs = mot[i:i+d_nuc]
                    if lastNucs == nucs:
                        most += 1
                    else:
                        most = 1
                        
                    lastNucs = nucs
                    if most >= n:
                        stop = True
                        break
                        
                if stop:
                    break

            if stop == False:
                motifsClean.append(mot)
    return motifsClean

def hashTable(sequences:list, kmersV:list):
    """
    Cherche les motifs de taille k dans un ensemble de séquences avec la méthode de Hash Table
    entrée séquences : liste de séquences
    entrée kmersV: liste de Kmers valides à chercher
    sortie motifs_freq_dict: dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    >>>hashTable(['TAAGTAA', 'TATAA', 'CTATC'], ['TAA', 'AAG', 'AGT', 'GTA', 'TAT', 'ATA', 'CTA', 'ATC'])
    {'TAA': 3, 'AAG': 1, 'AGT': 1, 'GTA': 1, 'TAT': 2, 'ATA': 1, 'CTA': 1, 'ATC': 1}
    """

    if len(kmersV) == 0:
        return {}
        
    motifs_freq_dict = {k: 0 for k in kmersV}
    k = len(kmersV[0])
    for seq in sequences:
        t = len(seq)
        for i in range(0, t-k+1):
            motif = seq[i:i+k]
            if motif in motifs_freq_dict:
                motifs_freq_dict[motif] += 1
                
    return motifs_freq_dict

def gatherRevCompMotifs(motifs_dict):
    """
    Trouver les motifs reverse complementaire et rassembler leur résultats
    entrée motifs_dict : dictionaire de résultats de recherche de motifs
    sortie motifs_dict_gather : dictionaire de motifs, clé = motif, valeur = valeur motif fw, valeur motif rv
    """
    #ATCC GGAT
    d = {}

    for k, v in list(motifs_dict.items()):
        rv = reversecompl(k)
        if rv in d:
            d[rv] = ((d[rv][0]), v)
        else:
            d[k] = (v, 0)

    motifs_dict_gather = {k: v for k, v in sorted(d.items(), key=lambda mrv: mrv[1][0] + mrv[1][1], reverse=True)}
    return motifs_dict_gather

def searchGivenMotif(motifsDict, motifSpecifique, decroissant = True):
    """
    Cherche un motif specifique dans un dictionnaire de motifs trouvés
    entrée motifsDict : dictionnaire de motifs, clé = motif, valeur = fréquence d'observation
    entrée motifSpecifique: un motif specifique à chercher
    entrée decroissant : bool, si True, le dictionnaire est trié par ordre décroissant de valeur
    sortie fréquence : la fréquence du motif
    sortie ranking : dans quelle position le motif a été trouvé
    >>>searchGivenMotif(test_motifs, "TAT")
    (2, 2)
    """
    
    ranking = 0
    frequence = 0
    motifsfreq = sorted(motifsDict, key=lambda m: motifsDict[m], reverse=decroissant)
    
    for motif in motifsfreq:
        ranking += 1
        if motif == motifSpecifique:
            frequence = motifsDict[motif]
            break
    
    return ranking, frequence

def readFasta(fastaFileName, upper=False):
    """
    Lit un fichier fasta
    entrée fastaFileName: nom du fichier fasta
    sortie séquences: liste contenant toutes les séquences du fichier
    """
    
    sequence = ""
    sequences_list = []
    prev_header = ""
    header = ""
        
    for line in open(fastaFileName):
        string = line.strip()
        if string[0] != ">":
            if prev_header != header:
                prev_header = header
            sequence = sequence + string
        else:
            header = string
            if sequence != "":
                if upper:
                    sequence = sequence.upper()
                sequences_list.append(sequence)
                sequence = ""

    sequences_list.append(sequence)
    return sequences_list

def modifierMotif(motif:str, v:int, upper = True):
    """
    Modifie v positions d'un motif aléatoirement 
    entrée motif: motif à modifier
    entrée v: nombre de positions à varier, si 0 le motif n'est pas modifié
    entrée upper: bool, si True les nucléotides modifiés seront en majuscule, si False en minuscule
    sortie motif_modifie: motif modifié
    """
    if v <= 0:
        return motif
        
    nuc = ('A', 'C', 'G', 'T')
    motif_modifie = list(motif)
    t = len(motif)
    for i in range(v):
        motif_modifie[random.randint(0, t-1)] = random.choice(nuc)

    motif_modifie = ''.join(motif_modifie)
    if upper == False:
        motif_modifie = motif_modifie.lower()
    
    return motif_modifie

def implantMotifsVar(motif:str, sequences:list, v:int, f1 = 1, f2 = 0.4, upper = True):
    """
    Insère un motif dans des positions aléatoires des séquences
    entrée motif : motif qui va être implanté dans les séquences
    entrée séquences : liste de séquences
    entrée v: nombre de positions à varier, si 0 le motif sera invariable
    entrée f : fréquence d'implantation, si 1 toutes les séquences contiendront un motif
    entrée f2 : fraction de motifs variables parmi les motifs implantés, si 0 tous les motifs implantés seront identiques
    entrée upper : bool, si False le motif sera en minuscules
    sortie modified_sequences: liste de séquences ayant le motif implanté
    """
    #Fonctions interne d'insertion de motifs
    def insertMotif(motif:str, seq:str, upper:bool):
        if upper == False:
            motif = motif.lower()
        tailleSeq = len(seq)
        index = random.randint(0,tailleSeq-1)
        return seq[0:index] + motif + seq[index:]
    
    def insertMotifVar(motif:str, seq:str, v:int, upper:bool):
        tailleSeq = len(seq)
        index = random.randint(0,tailleSeq-1)
        motifVar = modifierMotif(motif, v, upper)
        return seq[0:index] + motifVar + seq[index:]

from itertools import product

def createKmers(k, m, n, p, variation):
    """
    entrée k: taille de kmers
    entrée m: taille de repetition de nucléotide
    entrée n: taille de répétition de dinucléotide
    entrée p: proportion de nucleótides T et A
    entrée variation: si True permettre variation d'un nucléotide
    sortie: liste de motifs sans les motifs peu complexe
    """
    
    motifs = ["".join(s) for s in list(product('ATGC', repeat=k))]
    motifs = removeLowComplexeHomo(motifs, m)
    motifs = removeLowComplexeHetero(motifs, n, variation)
    return removeTARich(motifs, p)

def hamDistance(str1:str, str2:str):
    """
    Calcule la distance de Hamming entre deux chaînes de caractères
    entrée str1: chaîne de caractères
    entrée str2: chaîne de caractères
    sortie distance: distance de Hamming
    >>>hamDistance("TTGGTAT", "TTGCTAA")
    2
    """
    distance = 0
    t_str2 = len(str2)
    for i in range(len(str1)):
        if i < t_str2 and str1[i] != str2[i]:
            distance += 1
    
    return distance

def totalDistance(motif:str, sequences):
    """
    Calcul la totalDistance
    entrée motif: motif à comparer, chaîne de caractères
    entrée sequences: liste de séquences
    sortie total_distance: somme de distance de hamming minimal
    """
    total_distance = 0
    motif_rv = reversecompl(motif)
    for seq in sequences:
        seq = seq.upper()
        t_seq = len(seq)
        bestHam = 999999
        bestHamRv = 999999
        for i in range(t_seq):
            bestHam = min(hamDistance(motif, seq[i:]), bestHam)
            bestHamRv = min(hamDistance(motif_rv, seq[i:]), bestHamRv)

        total_distance += min(bestHam, bestHamRv)
    return total_distance

def medianStringSearch(sequences, kmersV):
    """
    Implemente l'algorithme MedianStringSearch
    entrée séquences : liste de séquences
    entrée kmersV: Liste de Kmers à chercher
    sortie motif_dist_dict: un dictionnaire contenant les motifs et leurs distances
    """
    motif_dist_dict = {}
    size = len(kmersV)
    for i in tqdm(range(size)):
        motif = kmersV[i]
        dist = totalDistance(motif, sequences)
        motif_dist_dict[motif] = dist
        
    return motif_dist_dict

def removeTARich(motifs:list, p:int):
    """
    Enlève les motifs contenant > p de T et A en combination
    entrée motifs: liste de motifs
    entrée p: proportion de nucleótides T et A
    sortie motifsClean: liste de motifs sans les motifs TA rich
    """
    motifsClean = []
    for m in motifs:
        taille = len(m)
        count = 0
        for i in range(taille):
            if m[i] in ['A','T']:
                count += 1
        if count/taille < p:
            motifsClean.append(m)
        
    
    return motifsClean

def writeMatchesFile(path:str, matches:list):
    with open(path, 'w') as f:
        f.write("\n".join(matches))
        return 1
    return 0

def getRatio(sequences:list):
    nuc = ["A","T","C","G"]
    tot = 0
    d = {k: 0 for k in nuc}
    for seq in sequences:
        for n in nuc:
            count = seq.count(n)
            d[n] += count
            tot += count

    return d, [d[i]*100/tot for i in d], tot 
    

    
