import sys
from timeit import Timer
from Bio import Seq
from Bio import SeqIO

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

    if number<4:
        return str(number)
    else:
        return dec2base4(number/4)+str(number%4)


def number2pattern(number,n):
    alphabet = {0:'A',1:'C',2:'G',3:'T'}
    str_num = ''.join(alphabet[int(s)] for s in str(dec2base4(number)))
    if len(str_num)<n:
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




if __name__=='__main__':

    lines = open('data.txt').readlines()
    k, d = lines[0].strip().split()
    k, d = int(k), int(d)
    dna = [l.strip() for l in lines[1:]]
    # dna = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
    print ' '.join(motifenumeration(dna, k, d))
    
    # print len(neighbors('CCCC',3))
    #     print n


    # gnm_seq = SeqIO.read('Salmonella_enterica.txt','fasta')
    # min_skew_coord = minimumskew(gnm_seq)[0]
    # window = gnm_seq[min_skew_coord-250:min_skew_coord+250]
    # print frequentwordsmismatch_rc(window,9,1)


    # lines = open('data.txt').readlines()
    # text=lines[0].strip()
    # k,d=lines[1].strip().split()
    # k,d = int(k), int(d)
    # print frequentwordsmismatch_rc(text, k, d)

    # lines = open('data.txt').readlines()
    # text=lines[0].strip()
    # word=lines[1].strip()
    # d=int(lines[2])
    # print approximatepatterncount('TACGCATTACAAAGCACA', 'AA', 1)

    # lines = open('data.txt').readlines()
    # pattern=lines[0].strip()
    # genome=lines[1].strip()
    # d=int(lines[2])
    # print approximatepattern(pattern, genome, d)

    # print pattern2number('ACTGCTCTAAACCTTGGT')
    # print number2pattern(8657,10)

    # text = 'GCGGTTATGCACCGTTCAAATTAGCAAACCACTAAGCGACGTAGTCTGGATTGATTTCTCCCTACCAGTGACCCAAGACGCGTTAGTGAGTTAAGTTCATATCCAGTACCTGCCGCCCTCTGTACTTGGGCGTCCGATTCGCATGCTTACTCAGGTGGAGGACACGATAATCTGATTAAACTGAGCTAAACCAGGTGGAACCAGAAACCAGGTGGGGAGTCTCGCTTCAAGCCGTTCTTGCGATCAAACCAGGTGGTCCATTATGAAACCAGGTGGCTAAACCAGGTGGTCCAGATCCTCGAATGATGTCGGTGCACATCAAAACCAGGTGGGGTGGTGGAACGTAAAACCAGGTGGCATAAACCAGGTGGGCCGGTTCGTAAACCAGGTGAAACCAGGTGGGGTGGAAACCAGGTGGGTTACAAATTACGTTGAGATGGCCCAAACCAGGTGGTGGGCTTCACCCATGTCAACAAACCACCCTATGGAACTAAACCAGGTGGAACCAGGTGGTGAAGGCTTATCCTCAGGAAAAACCAGGTGGAGGTGGTGAAATAAAACCAGGTGGACCAGGTGGATAACCCTCGCCTCGCTTCTCAACCGAGACCTGGATAAACCAGGTGGGGTGGTCCACCGATTTTTGAGACACTAGAAACCAGGTGGGCGGGGAAACCAGGTGGCAAACCAGGTGGGGTGGACGGAAACCAGGTGGATATGTCATAAAACCAAACCAGGTGGTGCACCCCCATGGTGTGTCTTATCCGTGCGTATAAACCAGGTGGTCGCACGGCTTCCACTTGCTGAGAATAGGCCCGCAGGGTCAGTGCCATGCCCTCCGTCACTCGATATGTGTTGTAAGAGTGGTTACCCCTTCATTGAAGTCGCCCACAGCCCCACCTGCATTGCTAGACTATCACCCTACAGTAGGCCTTTTCGCCTTCTTCAAGCAGCAATCTCTTATCCGCGGATGGGCGCGGCGAGCGTGGCGTCCCCGAACATTTTTACCTAACGTGTTTTGTTGGCCGCAAGCCTTCCCTCTAGTCCACCTCAGCCATTCAGCCTAGTAGCTTTCAAGCCGAGCCTTCCATATCTAATGGACCGTCCAGAATTTCACACGTTTCACAGGGCTGTGTTCGACCGCCCGTAATGCTGTTTCACAGGCGATCGCCTTGCGGTTTTTTCACAGATCGCAGCCGATGGACATGCCAACTCGATTTTCACAGAGTTTTTCACAGCGGTTTCACAGCACAGCAGTGATTGTTTCACAGCAATTTTCACTTTCACAGGGGCCCTTTTCACAGCTCAGGGCTCTTTTCACTTTCACAGTTTCACAGCGCTCCTTTCACAGAGCGGGGAAATTTAAGGGAACACTCAAGGGAACAAGGGAACACACAAAGGGAACACAACACAACACATAAGGGAACACTTTCACAGAACACAAAAGTCCGAAATCATCAGCGGCGAAGGGATTTCACAGACAGACACTTTCACAGCGCATTTCACAGATACGTACTTTCACAGGCGTACTTTCACAGACTTTCACAGAGGACAAGCTCAATTTTCACAGACAGGCTGGATAAATTTCACAGCGGTAAGGGTTTCACAGCACACATAAGGGAACACGAATTTCACAGCAGGGAACACCTCTACGAGTAATCTATTACTCTACCTACTGAAGGGAACACACCGAAGACCTACTATTACCTATTACTCTTAAAGGGAACACATTACAAGGGAACACACTCTCTCGTCATATCTCACCTCTCTATTACTCTTAAGGGAACACCTTCTCGATCAACCTATTACTCTATGGAGATAGAGATATTCCAGACATATGGAGATAACATGGAGATATGGAGATAATGGAGATGGAGATAGCTCTTATATTTATCCTATGGAGATATGATACTATTAATGGAGATAATTCTAATGGAGATATAATTACTCTAAGAGGATGGGATCTCGGGCTATTACTCTAATGGAGATAAGCACTATTACTCTAGGAAATGGAGATATGTCAATGGAGATATGTAATGGAGATAGAGGGAGATGGAGTCGCCATTTCATAATCGCCATTTCATAGTTCAGGAATCGCCATTTCCGCCATTTCTAAGATGGAGTCGCCATTTCTACGTATGGAGATAGGATCGCCATTTCATACGACCCGTTGGATATCGCCATTTCCTCGCCATTTCTGGTGACATTTCTCGCCATTTCATTTCTGGAGATAGATGGATCTCGCCATTTCATAGGAATCGCCATTTCCACGTAGGGGGGGCCACAATCCGTAGGTCGGAATTCAGACTCGCCATTTCCCATCGCCATTTCTTCACCTGTATGCCGATCCCTTCGCCATTTCTCATGGAGATAACTCTCTCTCGCCATTTCTCGCCATTTCCATTTCACTCTCATTCGCCATCGCCATTTCCATTCGCCATTTCATCGCCATTTCTTCAGGATAAGATATCGCCATTTCGACTCTCATTCGCATACTGACTCTCATTCTCATCTCGCCATTTCTCATCTGACTCTCATCCTGGGGGAAACTTGCGACTCTCATCACACTTCCGTCGACTCTCATACTGGCGGATAGCATAGGAGCCATTTAAAGACTCTCATTCTCATTCGAGACTCTCATTCAAATCCTACGAGGACTCTCATATAGACTCTCATATCATTACGAGGACTCTCATATACGAGCCATGCATGTGGCGACGACTCTCATCTACGAGCCATGCAAGCAGAATCTACGAGCGACTCTCATTACGAGCCATGTGACCGTACGAGCCATGCATGCATGCCATGCTGACTCTCATCGAGTACGAGCCATGGAAGTTCTTGTTGGTTCGTAGCCCAAGAGCTGAAGTTACGAGCCTACGAGCCATGAAGTTACTTTTACGAGCCATGAAGCTTACGATACGAGCCATGCGAGCCATGCATCCGCGCTACGAGCCATGTTCCAGTACGAGCCATGTTAGTTGCTGAAGTTAAGTTTGGCGCTGAAGTTTGTACGAGCCATGTGCCCGCTGAAGTTTGTTGTACGAGCCATGCATGCTGAAGTTAATGGCTGAAGTTAGCGTTTGCGGGCAGATCCTCATTCTACGATACGAGCCATGCCATGCAGCTGAAGTTAAGTTGGGTTACGAGCCATGCGAGCCATGTGAAGTACGAGCCATGCTGGCTGAAGTTGTTTGTGCTGCTGAAGTTGCTCTTGTCTCTAGCTGAAGTTGCCAACAGGGCTGAAGCTGAAGTTTAAGCTGAAGTTGCGAGCAGGCTGAAGTTATCGGATTGGGGCTGAAGTTCAACCTCCCGTCCCCCCACACTATATTCCCGTCCCCCCCCGCGCACGCGCCGTCTCCCGTCCCCCCTATCCCGTGCGCACGCGACGCGATCCCGTCCCCCCAGAGTGCGCGCACGCGTCCCCCTTCCCGTCCCCCTCTCCCGGGCGCACGCGTCGCTCAACATTTCCGCGCACGCGTCGCGCACGCGGGCGCACGCGGGTCCCGTCCCCCCCCCTCTTCGGCGCACGCGGAATTCCCGTCGCGCACGCGTCCCGTCCCGCGCACGCGTCGCGCACGCGACTGCCCTAACCAACAGTGCGCACGCGCCGGTAACCCGGTAACCCGGTAACCGCGCACGCGGGCGCACGCGCGTAACCCGCGCACGCGCCGCGCACGCGGCCCGGTTCCCGTCCCCCCCGGTAACCCGGTAACTCCCGTCCCCCGTAACCCGGTGCGCACGCGCCCGGCGCACGCGGAGCGCACGCGCCCCCCCCGGTAATAGCGCACGCGCCCGGGCGCACGCGCCCGGTAACCCGGTAACCCGGGCGCGCGCACGCGGCGGCGCACGCGGCGCACGCGGCGCACGCG'
    # k = 11
    # l = 566
    # t = 18
    # clumptask()

    # lines = open('data.txt').readlines()
    # text = lines[0].strip()
    # print text
    # k = int(lines[1].strip())
    # print k
    # freq = computingfrequencies(text, k)
    # print ' '.join(str(v) for v in freq.values())

    # Problem 1
    # lines = open('data.txt').readlines()
    # text= 'GACCATCAAAACTGATAAACTACTTAAAAATCAGT'
    # word= 'AAA'
    # print occurences(text, word)

    # Problem 2
    # lines = open('data.txt').readlines()
    # text= lines[0].strip()
    # k= int(lines[1].strip())
    # text= 'CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA'
    # k= 3
    # print " ".join(frequentwords(text, k))

    # text= 'CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA'
    # k= 3
    # print " ".join(fasterfrequentwords(text, k))

    # text = open('data.txt').readlines()[0].strip()
    # word = 'CTTGATCAT'
    # text = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
    # k,l,t = 9, 50, 4
    # findclumps(text,k,l,t)
    # print " ".join([str(p) for p in positions(text, word)])
    
    # text = open('data.txt').readline()
    # print " ".join([str(d) for d in minimumskew('GATACACTTCCCGAGTAGGTACTG')])

    # lines = open('data.txt').readlines()
    # seq1='ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'
    # seq2='CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG'
    # print hammingdistance(seq1, seq2)