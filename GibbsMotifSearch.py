# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 04:10:58 2016

@author: cgzheng86
"""

import sys
lines = sys.stdin.read().splitlines()

k = int(lines[0].split(" ")[0])
N = int(lines[0].split(" ")[1])
repeats = int(lines[0].split(" ")[2])
Dna = []
for i in range(1,len(lines)):
    Dna.append(lines[i])
 
def ReverseComplement(sequence):
    l = len(sequence)
    rev = ""
    rev_com = ""
    for i in range(0,l): 
        rev += sequence[l-i-1] #  making the reverse
    for i in range(0,l):
        if rev[i] == "A":
            rev_com += "T"
        elif rev[i] == "T":
            rev_com += "A"
        elif rev[i] == "C":
            rev_com += "G"
        elif rev[i] == "G":
            rev_com += "C"
        else:
            rev_com += "N"
    return (rev_com) 
    
def Score(Dna): # the lower the score the more conserved the Dna strings are
    from collections import Counter
    t = len(Dna)
    l = len(Dna[0])    
    P = 0
    for i in range(l):
        nucleotides = []
        for string in Dna:
            string = string.upper()
            nucleotides.append(string[i])
        counts = Counter(nucleotides)
        P = P + t - counts.most_common(1)[0][1] 
    return (P)

def KmerProb (pattern, matrix):
    pattern = pattern.upper()    
    probability = 1
    for i in range(len(pattern)):
        if pattern[i] == "A":
            probability = probability * matrix[0][i]
        elif pattern[i] == "C":
            probability = probability * matrix[1][i]
        elif pattern[i] == "G":
            probability = probability * matrix[2][i]
        elif pattern[i] == "T":
            probability = probability * matrix[3][i]
    return (probability)

def LaplaceProfileMatrix(Dna): # add 1 to all the counts to avoid the 0 probability
    l = len(Dna[0])
    t = len(Dna)     
    matrix = [[0 for i in range(l)] for i in range(4)]
    for i in range(l):
        nucleotides = []
        for string in Dna:
            string = string.upper()
            nucleotides.append(string[i])
        matrix[0][i] = (nucleotides.count("A") + 1)/(t + 4)
        matrix[1][i] = (nucleotides.count("C") + 1)/(t + 4)
        matrix[2][i] = (nucleotides.count("G") + 1)/(t + 4)
        matrix[3][i] = (nucleotides.count("T") + 1)/(t + 4)
    return (matrix)

def ProfileRandomKmer (string, k, matrix):
    l = len(string) 
    prob = []
    for i in range(l-k+1):
        prob.append(KmerProb(string[i:i+k], matrix))
    for i in range(l-k+1):
    	prob.append(KmerProb(ReverseComplement(string[i:i+k]), matrix))
    adjust_prob = [x/sum(prob) for x in prob]
    import numpy
    r = numpy.random.choice(range(2*(l-k+1)),p=adjust_prob)
    if r > (l-k):
    	r = r - (l-k+1)
    	result = ReverseComplement(string[r:r+k]) 
    else:
    	result = string[r:r+k]
    return (result)    

def GibbsSampler(Dna, k, N): # instead of randomly choose kmers for the entire motifs; it replaces one kmer randomly at a time
    t = len(Dna)
    l = len(Dna[0])
    motifs = []
    import random
    for string in Dna:
        i = random.randrange(0,l-k+1)
        motifs.append(string[i:i+k]) # randomly select motifs to start with
    BestMotifs = motifs
    for j in range(0, N):
        i = random.randrange(t) # choose the one to be replaced
        del motifs[i] # remove the ith kmer in the motif
        matrix = LaplaceProfileMatrix(motifs)
        motifs.insert(i, ProfileRandomKmer(Dna[i],k,matrix)) # replaced with a profile_randomly generated kmer from ith string
        if Score(motifs) < Score(BestMotifs):
            BestMotifs = motifs
    return (BestMotifs)

def GibbsMotifSearch (Dna, k, N, repeats):
    best_motifs = GibbsSampler(Dna, k, N)
    for i in range(repeats):
        new_motifs = GibbsSampler(Dna, k, N)
        if Score(new_motifs) < Score(best_motifs):
            best_motifs = new_motifs
        print (repeats - i)
    return (best_motifs)

results = GibbsMotifSearch(Dna, k, N, repeats)
for string in results:
    print (string)