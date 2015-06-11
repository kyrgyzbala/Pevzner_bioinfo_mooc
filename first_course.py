import sys
from timeit import Timer
from Bio import Seq
from Bio import SeqIO
from collections import Counter
import random

def occurences(text, word):
    count= 0
    wordlength= len(word)
    for i in range(len(text)-wordlength+1):
        if text[i:i+wordlength]==word:
            count+=1
    return count

def frequentwords(text, k):
    frequents=[]
    maxcnt=0
    for i in range(len(text)-k+1):
        word= text[i:i+k]
        cnt= occurences(text, word)

        if cnt>maxcnt:
            frequents=[word]
            maxcnt= cnt
        elif cnt==maxcnt:
            frequents.append(word)

    return set(frequents)


def fasterfrequentwords(text, k):
    frequents = []
    freqarray = computingfrequencies(text, k)
    maxcnt = max(freqarray.values())
    frequents = [number2pattern(n,k) for n,v in freqarray.items() if v==maxcnt]
    return frequents


def positions(text, word):
    poss=[]
    for i in range(len(text)-len(word)):
        if text[i:i+len(word)]==word:
            poss.append(i)
    return poss



def generate_kmer_map(k):
    alphabet=['A','T','G','C']
    kmer=[]
    key_list = []
    for i in range(k):
        if i==0:
            key_list = alphabet
        else:
            prev_list = [t for t in key_list if len(t)==i]
            new_list = [t+u for t in prev_list for u in alphabet]
            key_list = new_list

    kmers_map={key:0 for key in key_list}
    return kmers_map


def findclumps(text, k, l, t):

    kmer_map = generate_kmer_map(k)
    print len(kmer_map)


def pattern2number(pattern):

    alphabet = {'A':0,'C':1,'G':2,'T':3}
    num_str = ''.join(str(alphabet[s]) for s in pattern)
    decnum = int(num_str,4)
    return int(decnum)


def dec2base4(number):

    if number < 4:
        return str(number)
    else:
        return dec2base4(number/4)+str(number % 4)


def number2pattern(number, n):
    alphabet = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    str_num = ''.join(alphabet[int(s)] for s in str(dec2base4(number)))
    if len(str_num) < n:
        str_num = ''.join(alphabet[0] for i in range(n-len(str_num)))+str_num
    return str_num


def computingfrequencies(text, k):
    frequencyarray = {i:0 for i in range(4**k)}
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        j = pattern2number(pattern)
        frequencyarray[j]+=1
    return frequencyarray


def findclumps(genome, k, l, t):
    text = genome[:l]
    clumps = []

    freqs = computingfrequencies(text, k)
    clumps = [ number2pattern(key,k) for key,value in freqs.items() if value >=t ]

    for i in range(1, len(genome)-l+1):
        firstpattern = genome[i-1:i-1+k]
        freqs[pattern2number(firstpattern)] -= 1
        lastpattern = genome[i+l-k:i+l]
        freqs[pattern2number(lastpattern)]+=1
        if freqs[pattern2number(lastpattern)]>=t:
            clumps.append(lastpattern)

    return " ".join(s for s in set(clumps))

def clumptask():
    lines = open('data.txt').readlines()
    text = lines[0].strip()
    k,l,t = 9, 500, 3
    print findclumps(text, k, l, t)

def minimumskew(genome):
    minval,coords,gc=0,[],0
    for i in range(len(genome)):
        c=genome[i]
        if c=='G':
            gc+=1
        elif c=='C':
            gc-=1
        if gc<minval:
            minval=gc
            coords=[i+1]
        elif gc==minval:
            coords.append(i+1)
    return coords

def hammingdistance(seq1, seq2):
    if len(seq1)==len(seq2):
        return sum(ch1 != ch2 for ch1,ch2 in zip(seq1, seq2))
    else:
        raise ValueError("Two sequences are not equal")


def approximatepattern(pattern, genome, d):
    coords=[]    
    for i in range(len(genome)-len(pattern)+1):
        subs = genome[i:i+len(pattern)]
        if hammingdistance(subs ,pattern)<=d:
            coords.append(i)
    return " ".join([str(d) for d in coords])

def approximatepatterncount(text, word, d):
    count= 0
    wordlength= len(word)
    for i in range(len(text)-wordlength+1):
        if hammingdistance(text[i:i+wordlength],word)<=d:
            count+=1
    return count

def frequentwordsmismatch(text, k, d):
    frequents=[]
    maxcnt=0

    kmer_map = generate_kmer_map(k)

    for kmer in kmer_map.keys():

        kmer_map[kmer]= approximatepatterncount(text, kmer, d)

    maxval = max(kmer_map.values())
    frequents = [kmer for kmer,count in kmer_map.items() if count==maxval]

    return " ".join(set(frequents))

def frequentwordsmismatch_rc(text, k, d):
    frequents=[]

    kmer_map = {}
    for i in range(0,len(text)-k+1):
        word=text.seq[i:i+k]
        for n in neighbors(word,d):
            kmer_map[n]=0
    
    cnt=0
    for kmer in kmer_map.keys():
        tmpSeq = Seq.Seq(kmer)

        kmer_map[kmer]= approximatepatterncount(text, kmer, d) + \
                        approximatepatterncount(text, tmpSeq.reverse_complement(), d)
        cnt+=1
        if cnt%1000==0:
            print cnt
    
    maxval = max(kmer_map.values())
    frequents = [kmer for kmer,count in kmer_map.items() if count==maxval]

    return " ".join(set(frequents))

def neighbors(pattern,d):
    nctds=['A','T','C','G']
    if d==0:
        return [pattern]
    if len(pattern)==1:
        return nctds
    neighborhood=[]
    first_symbol = pattern[0]
    suffix=pattern[1:]
    suffix_neighbors = neighbors(suffix, d) 
    for neighbor in suffix_neighbors:
        if hammingdistance(neighbor, suffix)<d:
            for nctd in [n for n in nctds]:
                neighborhood.append(nctd+neighbor)
        else:
            neighborhood.append(first_symbol+neighbor)
    return neighborhood


def is_found(kmer, dna, d):
    for s in dna:
        if not approximatepattern(kmer, s, d):
            return False
    return True


def motifenumeration(dna, k, d):    
    patterns = []
    for seq in dna:
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            for kmerp in neighbors(kmer, d):
                if is_found(kmerp, dna, d):
                    patterns.append(kmerp)
    return set(patterns)


def distance_between_pattern_and_strings(pattern, dna):
    k = len(pattern)
    distance = 0

    for text in dna:
        min_ham_dist = 10e200 # very big number
        for i in range(len(text)-k+1):
            kmer = text[i:i+k]
            cur_ham_dist = hammingdistance(pattern, kmer)
            if cur_ham_dist < min_ham_dist:
                min_ham_dist = cur_ham_dist
        distance += min_ham_dist
    return distance


def median_string(dna, k):

    distance, median = 10e200, []

    for i in range(pow(4, k)):
        pattern = number2pattern(i, k)
        pattern_distance = distance_between_pattern_and_strings(pattern, dna)
        if not distance or distance > pattern_distance:
            distance = pattern_distance
            median = [pattern]
        elif distance == pattern_distance:
            median.append(pattern)
    return median, distance


def profile_probable(text, k, profile):
    prob_score, prob_kmer = 0, None
    t2n = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        tmp_score = 1
        for pos, c in enumerate(kmer):
            tmp_score *= profile[t2n[c]][pos]
        if not prob_kmer:
            prob_score = tmp_score
            prob_kmer = kmer
            continue

        if tmp_score > prob_score:
            prob_score = tmp_score
            prob_kmer = kmer

    return prob_kmer

def create_profile_matrix(motifs):
    profile = [[], [], [], []]
    letters = ['A', 'C', 'G', 'T']
    for i in range(len(motifs[0])):
        column = [l[i] for l in motifs]
<<<<<<< Updated upstream
        column_sum=0.0

        for j in range(4):
            column_sum += column.count(letters[j])+1

        for j in range(4):
            profile[j].append((column.count(letters[j])+1)/column_sum)

=======
        profile[0].append(column.count('A')+1)
        profile[1].append(column.count('C')+1)
        profile[2].append(column.count('G')+1)
        profile[3].append(column.count('T')+1)
>>>>>>> Stashed changes
    return profile


def build_consensus(motifs):
    cons = []
    for i in range(len(motifs[0])):
        column = [l[i] for l in motifs]
        cons.append(Counter(column).most_common()[0][0])
    return ''.join(cons)


def motif_score(motifs):
    consensus = build_consensus(motifs)
    score = 0
    for motif in motifs:
        score += hammingdistance(consensus, motif)

    return score


def greedy_motif_search(dna, k, t):

    best_motifs = []
    for strand in dna:
        best_motifs.append(strand[:k])

    base_strand = dna[0]
    other_strands = [l.strip() for l in dna[1:]]

    for i in range(len(base_strand)-k+1):
        kmer = base_strand[i:i+k]
        motifs = [kmer]

        for strand in other_strands:
            profile_matrix = create_profile_matrix(motifs)
            next_motif = profile_probable(strand, k, profile_matrix)
            motifs.append(next_motif)

        if motif_score(motifs)< motif_score(best_motifs):
            best_motifs = motifs
    return best_motifs


def get_motifs_using_profile(profile, dna):

    motifs = []
    k = len(profile[0])
    for strand in dna:
        motifs.append(profile_probable(strand,k,profile))
    return motifs


def randomized_motif_search(dna, k):

    best_motifs, motifs = [], []

    rand_ints = [random.randint(0, len(dna[0])-k-1) for i in range(len(dna))]
    motifs = [dna[i][p:p+k] for i,p in enumerate(rand_ints)]

    best_motifs = motifs
    while True:
        profile = create_profile_matrix(motifs)

        motifs = get_motifs_using_profile(profile, dna)

        if motif_score(motifs) < motif_score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def iterate_randomized_motif_search(dna, k, t):

    best_motif = None

    for i in range(1000):
        if not best_motif:
            best_motif = randomized_motif_search(dna, k)
            continue

        motifs = randomized_motif_search(dna, k)

        if motif_score(motifs)< motif_score(best_motif):
            best_motif = motifs

    return best_motif

if __name__=='__main__':


    # lines = open('data.txt').readlines()
    # k, t = lines[0].strip().split()
    # k, t = int(k), int(t)
    # dna = [l.strip() for l in lines[1:]]

    # for m in greedy_motif_search(dna, k, t):
    #     print m

<<<<<<< Updated upstream
    for motif in iterate_randomized_motif_search(dna, k, t):
        print motif
=======
    lines = open('data.txt').readlines()
    lines = [l.strip() for l in lines]

    for pattern in ['GTCAGCG', 'TAGTTTC', 'AACGCTG', 'ATAACGG', 'GTAGGAA', 'CGTGTAA']:
        print pattern, distance_between_pattern_and_strings(pattern, lines)

    # print median_string(lines,7)

>>>>>>> Stashed changes

