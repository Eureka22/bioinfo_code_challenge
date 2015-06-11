from collections import Counter
import operator
import math
import sys
import numpy as np
import random
import copy
def hamming(st1,st2):
    h = 0
    for i in range(len(st1)):
        if st1[i]!=st2[i]:
            h += 1
    return h
    
def approxmatching(st, substr, d):
    pos = []
    for i in range(len(st) - len(substr) + 1):
        if hamming(st[i:i+len(substr)], substr) <= d:
            pos.append(i)
    return pos 


def freqpattern_withmismatch(st,k,d):
    freq = {}
    for item in allpatterns(k):
        freq[item] = patterndcount(st,item,d)
    maxv = max(freq.values())
    res = [item for item in freq.keys() if freq[item] == maxv]
    return res

def freqpattern_withmismatch_reverse(st,k,d):
    freq = {}
    for item in allpatterns(k):
        freq[item] = patterndcount(st,item,d) + patterndcount(st,reversecomp(item),d)
    maxv = max(freq.values())
    res = [item for item in freq.keys() if freq[item] == maxv]
    return res



def patterncount(st, substr):
    count = 0
    for i in range(len(st)-len(substr)+1):
        if st[i:i+len(substr)] == substr:
            count += 1
    return count
    
def patterndcount(st, substr,d):
    count = 0
    for i in range(len(st)-len(substr)+1):
        if hamming(st[i:i+len(substr)],substr) <= d:
            count += 1
    return count
    
def freqkmer(st,n):
    kmer = []
    for i in range(len(st)-n+1):
        kmer.append(st[i:i+n])
    cnt = Counter(kmer)
    sorted_cnt = sorted(cnt.items(), key=operator.itemgetter(1))
    sorted_cnt.reverse()
    maxfreq = sorted_cnt[0]
    solution = []
    for item in sorted_cnt:
        if item[1] == maxfreq[1]:
            solution.append(item[0])
        else:
            break
    return solution
    
    
def freqtkmer(st,n,t):
    kmer = []
    for i in range(len(st)-n+1):
        kmer.append(st[i:i+n])
    cnt = Counter(kmer)
    sorted_cnt = sorted(cnt.items(), key=operator.itemgetter(1))
    sorted_cnt.reverse()
    solution = []
    for item in sorted_cnt:
        if item[1] >= t:
            solution.append(item)
        else:
            break
    return solution
    
    
def skew(st):
    s = 0
    res = []
    for i in st:
        res.append(s)
        if i=='G':
            s += 1
        elif i=='C':
            s -= 1
    res.append(s)
    return res
    
    
def minskew(st):
    res = skew(st)
    m = min(res)
    return [i for i,value in enumerate(res) if value == m]
    
def maxskew(st):
    res = skew(st)
    m = max(res)
    return [i for i,value in enumerate(res) if value == m]
    
    
def reversecomp(st):
    rev = {'A':'T', 'C':'G','T':'A', 'G':'C'}
    st2 = []
    for i in range(len(st)):
        st2.append(rev[st[i]])
    st2.reverse()
    return ''.join(st2)
    
    
    
def findpattern(st,substr):
    result = []
    for i in range(len(st)-len(substr)+1):
        if st[i:i+len(substr)] == substr:
            result.append(i)
    return result
    
def checkclumplength(indicies, t, L):
    for i in  xrange(len(indicies)-t+1):
		if indicies[t+i-1] - indicies[i] <= L:
			return True
    return False
   

def findpatternclump(st, k, L, t):
    result = []
    patterns = [item[0] for item in freqtkmer(st,k,t)]
    print patterns
    for pattern in patterns:
        indices = findpattern(st,pattern)
        print indices
        if checkclumplength(indices, t, L):
            result.append(pattern)
    return result
    
def allpatterns(d):
    di = {0:'A',1:'T',2:'C',3:'G'}
    n = 4**d
    patternpool = []
    for i in range(n):
        pattern = ''
        q = i
        for k in range(d):
            m = q % 4
            q = q / 4
            pattern += di[m]
        patternpool.append(pattern)
    return patternpool
    
    
def motifenumeration(dnas,k,d):
    candidate = allpatterns(k)
    result = []
    #for dna in dnas:
    #    for i in range(len(dna) - k+1):
    #        candidate.append(dna[i:i+k])
    #print candidate
    
    for pattern in  set(candidate):
        select = True
        for dna in dnas:
            if patterndcount(dna,pattern,d) == 0:
                select = False
        if select:
            result.append(pattern)
    return result

def patterndistance(pattern, dna):
    distance = []
    for i in range(len(dna) - len(pattern) +1):
        distance.append(hamming(dna[i:i+len(pattern)], pattern))
    #print distance
    return min(distance)

def patterndistance_dnas(pattern, dnas):
    dis = 0
    for dna in dnas:
        dis += patterndistance(pattern,dna)
    return dis
    
    
    
def median_string(dnas, k):
    distance = sys.maxint
    for pattern in allpatterns(k):
        d = patterndistance_dnas(pattern, dnas)
        #print pattern, d
        if d < distance:
            distance = d
            median = pattern
    return median
    
def profile_fit(dna, k, profile):
    d = {'A':0,'C':1,'G':2,'T':3}
    probs = []
    patterns = []
    for i in range(len(dna)-k+1):
        pattern = dna[i:i+k]
        prob = 1
        for j,item in enumerate(pattern):
            prob *=  profile[d[item],j]
        probs.append(prob)
        patterns.append(pattern)
    
    kmer = patterns[probs.index(max(probs))]
    return kmer
    
    
def profile_prob(dna, profile):
    d = {'A':0,'C':1,'G':2,'T':3}
    probs = 1
    for i,item in enumerate(dna):
        prob *=  profile[d[item],i]
    return prob
    
def profile_generate(prof, dna):
    d = {'A':0,'C':1,'G':2,'T':3}
    probs = []
    patterns = []
    for i in range(len(dna)-k+1):
        pattern = dna[i:i+k]
        prob = 1
        for j,item in enumerate(pattern):
            prob *=  prof[d[item],j]
        probs.append(prob)
        patterns.append(pattern)
    #kmer = patterns[probs.index(max(probs))]
    
    probs /= np.sum(probs)
    #print probs
    #print patterns
    rnd = random.random()
    for i in range(len(probs)):
        #print rnd
        rnd -= probs[i]
        if rnd<0:
            break
    kmer = patterns[i] 
    return kmer
    
#def profile_generate(prof):
#    d = {'A':0,'C':1,'G':2,'T':3}
#    di = {0:'A',1:'C',2:'G',3:'T'}
#    st = ''
#    for i in range(np.shape(prof)[1]):
#        x = random.random()
#        for j in range(4):
#            x -= prof[j,i]
#            if x<0 :
#                break
#        st += di[j]
#    return st
    

def profile(dnas):
    d = {'A':0,'C':1,'G':2,'T':3}
    n = len(dnas)
    profile = np.zeros((4,len(dnas[0])))
    for i in range(len(dnas[0])):
        c = []
        for j in range(len(dnas)):
            c.append(dnas[j][i])        
        cnt = Counter(c)
        
        for k,key in enumerate(['A','C','G','T']):
            if cnt.has_key(key):
                #print cnt[key],k
                profile[k,i] = cnt[key]/float(n)
            else:
                profile[k,i] = 0
    return profile
    
    
def profile_pesudo(dnas):
    d = {'A':0,'C':1,'G':2,'T':3}
    n = len(dnas)
    profile = np.zeros((4,len(dnas[0])))
    for i in range(len(dnas[0])):
        c = ['A','C','G','T']
        for j in range(len(dnas)):
            c.append(dnas[j][i])        
        cnt = Counter(c)
        
        for k,key in enumerate(['A','C','G','T']):
            if cnt.has_key(key):
                #print cnt[key],k
                profile[k,i] = cnt[key]/float(n+4)
            else:
                profile[k,i] = 0
    return profile

def score(dnas):
    d = {'A':0,'C':1,'G':2,'T':3}
    n = len(dnas)
    score = 0
    for i in range(len(dnas[0])):
        c = []
        for j in range(len(dnas)):
            c.append(dnas[j][i])        
        cnt = Counter(c)
        #print cnt
        score += n - max(cnt.values())
    return score
    



def greedymotifsearch(dnas,k,t):
    bestmotifs = [item[0:k] for item in dnas]
    sc = score(bestmotifs)
    for i in range(len(dnas[0]) - k + 1 ):
        motif = [dnas[0][i:i+k]]
        for j in range(1,t):
            pf = profile(motif)
            newmotif = profile_fit(dnas[j],k,pf)
            motif.append(newmotif)
        newscore = score(motif)
        if newscore < sc:
            sc = newscore
            bestmotifs = motif
    return bestmotifs      
    
    
def greedymotifsearch_pesudo(dnas,k,t):
    bestmotifs = [item[0:k] for item in dnas]
    sc = score(bestmotifs)
    for i in range(len(dnas[0]) - k + 1 ):
        motif = [dnas[0][i:i+k]]
        for j in range(1,t):
            pf = profile_pesudo(motif)
            newmotif = profile_fit(dnas[j],k,pf)
            motif.append(newmotif)
        newscore = score(motif)
        if newscore < sc:
            sc = newscore
            bestmotifs = motif
    return bestmotifs                
                
                
def randomizedmotifsearch_pseudo(dnas,k,t,maxiter = 1000):
    allbest = None
    allbestscore = sys.maxint
    #candidates = []
    #for dna in dnas:
    #    for i in range(len(dna) - k + 1):
    #        candidates.append(dna[i:i+k])
    for itera in range(maxiter):
        print itera
        for dna in dnas:
            i = random.randint(0,len(dna) - k)
            candidates.append(dna[i:i+k])
        motifs = copy.deepcopy(candidates)
        
        #motifs = [candidates[i] for i in np.random.choice(len(candidates), t)]
        #print motifs
        bestmotifs = motifs
        while(1):
            pf = profile_pesudo(motifs)
            motifs = [profile_fit(item,k,pf) for item in dnas]
            #print motifs
            if score(motifs) < score(bestmotifs):
                bestmotifs = motifs
            else:
                break
        
        if allbest == None or score(bestmotifs) < score(allbest):
            allbest = bestmotifs
            allbestscore = score(bestmotifs)
    return allbest


def gibbssampler(dnas,k,t,N,maxiter = 50):
    allbest = None
    allbestscore = sys.maxint
        
    for itera in range(maxiter):
        print itera
        candidates = []
        for dna in dnas:
            i = random.randint(0,len(dna) - k)
            candidates.append(dna[i:i+k])
        motifs = copy.deepcopy(candidates)
        print motifs
        bestmotifs = copy.deepcopy(motifs)
        
        for j in range(N):
            i = random.randint(0,t-1)
            items = [item[1] for item in enumerate(motifs) if item[0] != i]
            pf = np.array(profile_pesudo(items))
            #print 1,score(bestmotifs)
            
            gpattern = profile_generate(pf,dnas[i])
            motifs[i] = gpattern
            if score(motifs) < score(bestmotifs):
                bestmotifs = copy.deepcopy(motifs)
        
        if allbest == None or score(bestmotifs) < score(allbest):
            allbest = bestmotifs
            allbestscore = score(bestmotifs)
        print allbestscore
    return allbest




if __name__=='__main__':
    
    #challenge 1
    #A = 'TATAGTAAAGTATAGTAAAGGTATAGTATTATAGTAGTTATAGTATATAGTACATATAGTATATAGTACTATAGTAGGTGTATATAGTATATAGTACTATAGTATATAGTATATATAGTAAATATATAGTACCTTGTCATATAGTAGTATAGTATATAGTACTTATAGTATATAGTAGTATAGTATATAGTATATAGTATATAGTAGTATAGTATATAGTATATAGTATTATAGTACAGCTTATAGTATGACCTATAGTATGATATATAGTAGTATAGTACATACACCGCATTTATAGTAATATAGTAAACAAATTATATAGTATAATATAGTAATATAGTAGAATAATCTGAGATATAGTATCCTATAGTAGTTATAGTATATAGTATATAGTATATAGTATATAGTATTTATAGTACTATAGTATATAGTATATAGTAGATATAGTACAGTATAGTATATAGTATCTTATAGTATCCGAGTGTTAACTATAGTAACGGCCATTCGAAAAGTATAGTAGTACGTATAGTAACGAACAAAAGTATAGTATTTATAGTATATAGTAATTATAGTAGTATAGTATATAGTAGCCGCTTTATAGTATATAGTAGGTATATAGTACGCGTGATACTATAGTATGTCTATAGTATATAGTAACTTATAGTACTTATAGTATTATAGTAGCCTATAGTATATAGTACCTATAGTATATTATAGTAGCCAATATAGTAGTATAGTACGCTCTTATAGTACGGACCTATAGTATTTTATAGTAATATAGTACTATAGTACCTTTATAGTACCTTATAGTATATAGTAGTATAGTAGATATAGTAGATCTATAGTAGTATAGTAATATAGTACTATAGTAATATAGTAGGCCTTTTATAGTAGATTAAGTATAGTATATAGTATCACAGGGGTATAGTACTAAGTCTTATAGTATTATAGTA'
    #B = 'TATAGTATA'
    #print patterncount(A,B)
    
    #challenge 2
    #A = 'CCATCCCTCCATCCCTCGGACTGCTCGGACTGCTTCCTTGAGTCCTTGAGTAGTATGGCGGACTGCTTAGTATGGCGGACTGCTTCCTTGAGCGGACTGCTTCCTTGAGCCATCCCTCCATCCCTCCCTAGAGTCCTTGAGCCATCCCTCCCTAGAGCCCTAGAGCCCTAGAGCGGACTGCTCGGACTGCTCCATCCCTTCCTTGAGTAGTATGGCCCTAGAGTCCTTGAGCGGACTGCTTAGTATGGCGGACTGCTCGGACTGCTCCATCCCTTCCTTGAGCGGACTGCTCCCTAGAGCGGACTGCTCCATCCCTCCCTAGAGTCCTTGAGTAGTATGGCCATCCCTCCATCCCTTCCTTGAGTCCTTGAGTCCTTGAGCCATCCCTCCCTAGAGCCCTAGAGCCCTAGAGTAGTATGGTAGTATGGCGGACTGCTCCCTAGAGCGGACTGCTCGGACTGCTCCCTAGAGCGGACTGCTTCCTTGAGCCATCCCTCCCTAGAGCGGACTGCTTAGTATGGCCCTAGAGTAGTATGGTAGTATGGCGGACTGCTCCCTAGAGCCATCCCTTAGTATGGTAGTATGGCCATCCCTCCATCCCTTAGTATGGTAGTATGGCCCTAGAGTAGTATGGTCCTTGAGCCATCCCTCCATCCCTCCCTAGAGCGGACTGCTCCCTAGAGCCATCCCTCCATCCCTTCCTTGAGTCCTTGAGTCCTTGAGTAGTATGGCGGACTGCTCCATCCCTTCCTTGAGCCATCCCTTCCTTGAGTCCTTGAGCCCTAGAGCCATCCCTCCCTAGAGCCATCCCTCCCTAGAGCCCTAGAGCCCTAGAG'
    #print freqkmer(A,13)
    
    #challenge 3
    #print reversecomp(A)
    
    #challenge 4
    #A = ''
    #B = 'GACGGTTGA'
    #print ' '.join([str(item) for item in findpattern(A,B)])
    
    #with open('Vibrio_cholerae.txt','r') as f:
    #    line = f.readline()
    #f.close()
    #A = line
    #B = 'CTTGATCAT'
    #print ' '.join([str(item) for item in findpattern(A,B)])
    
    #challenge 5
    #with open('E-coli.txt','r') as f:
    #    line = f.readline()
    #f.close()
    #A = line
    #print line
    #print ' '.join(findpatternclump(A,9,500,3))
    
    #challenge 6
    #print ' '.join(map(str,skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')))
    
    #challenge 7
    #print ' '.join(map(str,minskew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT')))
    
    #challenge 8
    #print hamming('GGGCCGTTGGT','GGACCGTTGAC')
    
    #challenge 9
    #print ' '.join(map(str,approxmatching('CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT', 'ATTCTGGA', 6)))
    
    #challenge 10
    #print patterndcount('AACAAGCTGATAAACATTTAAAGAG', 'AAAAA', 2)
    
    #challenge 11
    #print ' '.join(freqpattern_withmismatch(st,7,2))
    
    #challenge 12
    #print ' '.join(freqpattern_withmismatch_reverse(st,6,3))
    
    #challenge 13
    #dnas = ['ATTT','TGCC','CGGT','GAAA']
    #print ' '.join(motifenumeration(dnas,5,1))
    
    #challenge 14
    #print patterndistance_dnas('AAA',dnas)
    
    #dnas = ['TATACTGACCAATAACGGTAGGCGGTACGCGCACGGATGCGA','TAGACTGTACGGCCTATCTCACAGAGGTCGCGTAATAGTACA','ACAGAGATCGGCGCGGAGGTACGTCATTACAGAAGGGAAGTG','CACCTTATCCTCTTGCTCCTTCTAATGATCGATATTGTACGA','ATGGATTGTAACAACTCCGGGTTTACTCTTGTACGAAAACGA','AGAAACCTAGTCCAGTAGGTACGCTCAATTCTGTGTGCATCA','CAGTTGACTGATTCGGCGGAGAGAGTACGCGGTAACAGTGCC','TCATTTCGCACCCAAAGATTTTGAGTACGTACATCTGGTGAG','GTACGCCCCAAATCGAACGTGCATCGCGCTATTCGGGGTGGG','TGTGAACTTAATGTAAGTGTACGCATAACTCGTTCCCTGGTT']
    #print median_string(dnas, 6)
    
    #challenge 15
    #with open('test.txt','r') as f:
    #    dna = f.readline().strip()
    #    k = int(f.readline().strip())
    #    profile = f.readlines()
    #profile = np.array([[float(st) for st in item.strip().split()] for item in profile])
    #print profile_fit(dna,k,profile)
    
    #challenge 16
    #with open('test.txt','r') as f:
    #    k,t = map(int, f.readline().strip().split(' '))
    #    dnas = [item.strip() for item in f.readlines()]
    #    
    ##print k,t
    ##print dnas
    #print '\n'.join(greedymotifsearch_pesudo(dnas,k,t))
    
    
    #challenge 17
    #with open('test.txt','r') as f:
    #    k,t = map(int, f.readline().strip().split(' '))
    #    dnas = [item.strip() for item in f.readlines()]
    #print '\n'.join(randomizedmotifsearch_pseudo(dnas,k,t))
    
    #challenge 18
    #print profile_generate(profile_pesudo(['AATCG','AAAAA']))
    with open('test.txt','r') as f:
        k,t,N = map(int, f.readline().strip().split(' '))
        dnas = [item.strip() for item in f.readlines()]
    print k,t,N
    print '\n'.join(gibbssampler(dnas,k,t,N))
    
    #print score(['TCTCGGGG','CCAAGGTG','TACAGGCG','TTCAGGTG','TCCACGTG'])