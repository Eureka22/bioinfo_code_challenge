from collections import Counter
import operator
import math
import sys
import numpy as np
import random
import copy
import networkx as nx

def read_condon_table():
    d = {}
    with open('table.txt') as f:
        lines = f.readlines()
        for line in lines:
            a = line.strip().split(' ')
            if len(a)>1:
                d[a[0].replace('U','T')] = a[1]
    return d
            
            
def translate(dna):
    global d
    prot = ''
    for i in range(len(dna)/3):
        codon = dna[i*3:i*3+3]
        if d.has_key(codon):
            prot += d[codon]
    return prot
            
            
def translate_rna(rna):
    dna = rna.replace('U','T')
    global d
    prot = ''
    for i in range(len(dna)/3):
        codon = dna[i*3:i*3+3]
        if d.has_key(codon):
            prot += d[codon]
    return prot

def reversecomp(st):
    rev = {'A':'T', 'C':'G','T':'A', 'G':'C'}
    st2 = []
    for i in range(len(st)):
        st2.append(rev[st[i]])
    st2.reverse()
    return ''.join(st2)
    
def peptide_search(dna,pep):

    l = len(pep)*3
    candidate = [dna[i:i+l] for i in range(len(dna)-l)]
    reverse = reversecomp(dna)
    candidate2 = [reverse[i:i+l] for i in range(len(dna)-l)]
    res = []
    for item in candidate:
        if translate(item) == pep:
            res.append(item)
    for item in candidate2:
        if translate(item) == pep:
            res.append(reversecomp(item))
    return res
    
def read_mass_table():
    d = {}
    with open('masstable.txt') as f:
        lines = f.readlines()
        for line in lines:
            a = line.strip().split(' ')
            d[a[0]] = int(a[1])
    return d
    

def mass_spectrum_gen(prot):
    l = []
    mass = 0
    for i in range(len(prot)):
        l.append(mass)
        mass += m[prot[i]]
    l.append(mass)
    maxmass = mass
    res = []
    for i in range(1,len(prot)):
        for j in range(len(prot)):
            res.append( (l[(j+i)%len(prot)] - l[j] + maxmass) % maxmass )
    res.append(0)
    res.append(maxmass)
    res.sort()
    return res
    
    
def mass_spectrum_gen_num(prot):
    l = []
    mass = 0
    for i in range(len(prot)):
        l.append(mass)
        mass += prot[i]
    l.append(mass)
    maxmass = mass
    res = []
    for i in range(1,len(prot)):
        for j in range(len(prot)):
            res.append( (l[(j+i)%len(prot)] - l[j] + maxmass) % maxmass )
    res.append(0)
    res.append(maxmass)
    res.sort()
    return res
    
   
def mass_spectrum_gen_num_linear(prot):
    l = []
    mass = 0
    for i in range(len(prot)):
        l.append(mass)
        mass += prot[i]
    l.append(mass)
    maxmass = mass
    res = []
    for i in range(1,len(prot)):
        for j in range(len(prot)-i):
            res.append(l[j+i] - l[i])
    res.append(0)
    res.append(maxmass)
    res.sort()
    return res 
    
def mass_spectrum_gen_linear(prot):
    l = []
    mass = 0
    for i in range(len(prot)):
        l.append(mass)
        mass += m[prot[i]]
    l.append(mass)
    maxmass = mass
    res = []
    for i in range(1,len(prot)):
        for j in range(len(prot)-i):
            res.append(l[j+i] - l[i])
    res.append(0)
    res.append(maxmass)
    res.sort()
    return res 
    
    
def comp_comb(n):
    d2 = {}
    for key in m.keys():
        d2[m[key]] = 1
    d2[113] += 1
    d2[128] += 1
    
    print d2
    a = [0]*(n+1)
    for i in range(200):
        if i in d2.keys():
            a[i] = d2[i]
    
    for i in range(200):
        print i,a[i]
    
    
    for i in range(50,n+1):
        for key in d2.keys():
            if (i-key > 50) and (d2[key] > 0):
                a[i] += a[i-key]*d2[key]
    return a
    
    
def cyclosequencing(spectrum):
    s = set(m.values())
    sset = set(spectrum)
    print sset
    res = []   
            
    def fil(item):
        flag = True
        print item, mass_spectrum_gen_num_linear(item)
        for num in mass_spectrum_gen_num_linear(item):
            if not num in sset:
                flag = False
                break
        return flag
            
    def fil2(item):
        m = max(sset)
        return not sum(item) == m
        
    def fil3(item):
        m = max(sset)
        return  sum(item) == m
        
    for item in s:
        if item in sset:
            res.append([item])
    print res
    while len(res) >0:
        newres = []
        for j in range(len(res)):
            for i in s:
                if i in sset:
                    newres.append(res[j]+[i])
        res = newres
        # prune
        res = filter(fil,res)        
        # prune2
        candidate = []
        c = filter(fil3,res)
        for item in c:
            candidate.append(item)
        res = filter(fil2,res)
        print res
    return candidate
    
def spec_conv(spec):
    spec.sort()
    res = []
    for i in range(len(spec)):
        for j in range(i):
            res.append(spec[i] - spec[j])
    res.sort()
    return res    
    
def score(spec,seq):
    spec2 = mass_spectrum_gen(seq)
    c1 = Counter(spec)
    c2 = Counter(spec2)
    score = 0
    for i in set(c1.keys())&set(c2.keys()):
        score += min(c1[i],c2[i])
    return score
    
    
def score_linear(spec,seq):
    spec2 = mass_spectrum_gen_linear(seq)
    c1 = Counter(spec)
    c2 = Counter(spec2)
    score = 0
    for i in set(c1.keys())&set(c2.keys()):
        score += min(c1[i],c2[i])
    return score
    
    
def score_num(spec,seq):
    spec2 = mass_spectrum_gen_num(seq)
    c1 = Counter(spec)
    c2 = Counter(spec2)
    score = 0
    for i in set(c1.keys())&set(c2.keys()):
        score += min(c1[i],c2[i])
    return score
    
        
def leaderboardsequencing(spectrum,N):
    mx = max(spectrum)
    def fil(item):
        return sum(item) <= mx
        
    leaderboard = [[]]
    leaderpiptide = []
   
    s = set(m.values())
    while len(leaderboard) > 0:
        newleaderboard = []
        for j in range(len(leaderboard)):
            for i in s:
                newleaderboard.append(leaderboard[j]+[i])
        leaderboard = newleaderboard
        
        for peptide in leaderboard:
            if sum(peptide) == mx:
                if score_num(spectrum,peptide) > score_num(spectrum, leaderpiptide):
                    leaderpiptide = peptide
        leaderboard = filter(fil,leaderboard)
        pairlb = [(item,score_num(spectrum, item)) for item in leaderboard]
        pairlb.sort(key=lambda x: x[-1])
        #print pairlb
        leaderboard = [item[0] for item in pairlb[::-1][0:N]]
        for item in leaderboard:
            print item,score_num(spectrum,item)
        if len(leaderboard)>0:
            print sum(leaderboard[0]), mx
        #print leaderboard         
    return leaderpiptide
    
    
def convleaderboardsequencing(spectrum,M,N):
    conv = spec_conv(spectrum)
    def fil(item):
        return item>=57 and item <=200
    conv = filter(fil,conv)
    c = Counter(conv).items()
    c.sort(key=lambda x: x[-1])
    c = c[::-1]
    c = [item[0] for item in c]
    s = set(c[0:M])
    print s
    mx = max(spectrum)
    def fil(item):
        return sum(item) <= mx
        
    leaderboard = [[]]
    leaderpiptide = []
   
    #s = set(m.values())
    while len(leaderboard) > 0:
        newleaderboard = []
        for j in range(len(leaderboard)):
            for i in s:
                newleaderboard.append(leaderboard[j]+[i])
        leaderboard = newleaderboard
        
        for peptide in leaderboard:
            if sum(peptide) == mx:
                if score_num(spectrum,peptide) > score_num(spectrum, leaderpiptide):
                    leaderpiptide = peptide
        leaderboard = filter(fil,leaderboard)
        pairlb = [(item,score_num(spectrum, item)) for item in leaderboard]
        pairlb.sort(key=lambda x: x[-1])
        #print pairlb
        leaderboard = [item[0] for item in pairlb[::-1][0:N]]
        for item in leaderboard:
            print item,score_num(spectrum,item)
        #if len(leaderboard)>0:
        #    print sum(leaderboard[0]), mx
        #print leaderboard         
    return leaderpiptide
        
    

def spec_conv(spec):
    spec.sort()
    res = []
    for i in range(len(spec)):
        for j in range(i):
            res.append(spec[i] - spec[j])
    res.sort()
    return res

global d
d = read_condon_table()
global mass
m = read_mass_table()
print m


if __name__ == '__main__':
    
    #with open('ch1.txt') as f:
    #    line = f.readline().strip()
    #print translate(line,d)
    
    
    #with open('ch2.txt') as f:
    #    dna = f.readline().strip()
    #    pep = f.readline().strip()
    #print '\n'.join(peptide_search(dna,pep))
    
    #print ' '.join([str(item) for item in mass_spectrum_gen('SHSKHCELLWQEWW')])
    
    #with open('ch3.txt') as f:
    #    spec = [int(item) for item in f.readline().strip().split(' ')]
    #for item in cyclosequencing(spec):
    #    print '-'.join([str(i) for i in item])
        
    #with open('ch4.txt') as f:
    #    seq = f.readline().strip()
    #    spec = [int(item) for item in f.readline().strip().split(' ')]
    #print score(spec,seq)
    
    #with open('ch5.txt') as f:
    #    spec = [int(item) for item in f.readline().strip().split(' ')]
    #    
    #print ' '.join([str(item) for item in spec_conv(spec)])
    
    #with open('ch6.txt') as f:
    #    N = int(f.readline().strip())
    #    spec = [int(item) for item in f.readline().strip().split(' ')]
    #print '-'.join([str(item) for item in leaderboardsequencing(spec,N)])
    
    #with open('ch7.txt') as f:
    #    M = int(f.readline().strip())
    #    N = int(f.readline().strip())
    #    spec = [int(item) for item in f.readline().strip().split(' ')]
    #print '-'.join([str(item) for item in convleaderboardsequencing(spec,M,N)])