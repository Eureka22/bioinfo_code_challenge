from collections import Counter
import operator
import math
import sys
import numpy as np
import random
import copy
import networkx as nx

def composition(dna, k):
    return sorted(list(set([dna[i:i+k] for i in range(len(dna)-k+1)])))

def simplespell(dnas):
    l = len(dnas[0])
    res = dnas[0]
    for i in range(1,len(dnas)):
        res += dnas[i][-1]
        
    return res

def graph(dnas):
    res = []
    for i in range(len(dnas)):
        for j in range(len(dnas)):
            if dnas[i][1:-1] == dnas[j][0:-2]:
                res.append('%s %s %s\n'%(dnas[i],'->',dnas[j]))
    return res
    
def debrujin(dna, k):
    kmers = []
    res = []
    for i in range(len(dna)-(k-1)+1):
        kmers.append(dna[i:i+k-1])
    for kmer in sorted(list(set(kmers))):
        c = []
        for i in range(len(dna) - k +1):
            kkmer = dna[i:i+k]
            if kkmer[0:-1] == kmer:
                c.append(kkmer[1:k+1])
        if len(c)>0:
            res.append('%s %s %s'%(kmer,'->',','.join(sorted(c))))
    return res
    
    

    
def debrujin_patterns(patterns):
    k = len(patterns[0])
    kmers = []
    res = []
    for pattern in patterns:
        kmers.append(pattern[0:-1])
        kmers.append(pattern[1:k])
    
    for kmer in sorted(list(set(kmers))):
        c = []
        for kkmer in patterns:
            if kkmer[0:-1] == kmer:
                c.append(kkmer[1:k+1])
        if len(c)>0:
            res.append('%s %s %s'%(kmer,'->',','.join(sorted(c))))
    return res



def debrujin_patterns_pairs(patterns):
    k = len(patterns[0][0])
    kmers = []
    res = []
    for pattern in patterns:
        #print pattern
        kmers.append((pattern[0][0:-1],pattern[1][0:-1]))
        kmers.append((pattern[0][1:k],pattern[1][1:k]))
    
    #print kmers
    
    for kmer in sorted(list(set(kmers))):
        c = []
        for kkmer in patterns:
            if (kkmer[0][0:-1],kkmer[1][0:-1]) == kmer:
                c.append((kkmer[0][1:k+1], kkmer[1][1:k+1]))
        if len(c)>0:
            #res.append('%s %s %s'%(kmer,'->',','.join(sorted(c))))
            res.append([kmer, c])
    return res


def eulerian_circuit(G, source=None):
    if not nx.is_eulerian(G):
        raise nx.NetworkXError("G is not Eulerian.")

    g = G.__class__(G) # copy graph structure (not attributes)

    # set starting node
    if source is None:
        v = next(g.nodes_iter())
    else:
        v = source

    while g.size() > 0:
        print g.size()
        n = v   
        # sort nbrs here to provide stable ordering of alternate cycles
        nbrs = sorted([v for u,v in g.edges(n)])
        for v in nbrs:
            g.remove_edge(n,v)
            bridge = not nx.is_connected(g.to_undirected())
            if bridge:
                g.add_edge(n,v)  # add this edge back and try another
            else:
                break  # this edge is good, break the for loop 
        if bridge:
            g.remove_edge(n,v)            
            g.remove_node(n)
        yield (n,v)


def baseN(num,b,numerals="0123456789abcdefghijklmnopqrstuvwxyz"):
    return ((num == 0) and numerals[0]) or (baseN(num // b, b, numerals).lstrip(numerals[0]) + numerals[num % b])


if __name__ == '__main__':
    
    # challenge 1
    #with open('test.txt','r') as f:
    #    dnas = [item.strip() for item in f.readlines()]
    #print simplespell(dnas)
    
    # challenge 2
    #with open('test.txt','r') as f:
    #    dnas = [item.strip() for item in f.readlines()]
    #f.close()
    #with open('res.txt','w') as f:
    #    f.write()
    #f.close()
    
    #challenge 3
    #with open('res.txt','w') as f:
    
    #challenge 4
    #with open('test.txt','r') as f:
    #    dnas = [item.strip() for item in f.readlines()]
    #print '\n'.join(debrujin_patterns(dnas))
    #f.close()
    #with open('res.txt','w') as f:
    #    f.write('\n'.join(debrujin_patterns(dnas)))
    #f.close()
    
    
    #with open('test.txt','r') as f:
    #    lines = [item.strip() for item in f.readlines()]
    #f.close()
    #g = nx.DiGraph()
    #for line in lines:
    #    st = line.split(' -> ')
    #    for item in st[1].split(','):
    #        g.add_edge(st[0], item)
    #        print st[0],item
    #f = open('res.txt','w')       
    #edges = []
    #deg = g.degree(nx.nodes(g))
    #special = [item for item in  deg.keys() if deg[item]%2 == 1]
    #print special
    #if g.in_degree(special[0]) >= g.in_degree(special[1]):
    #    st = special[1]
    #    g.add_edge(special[0],special[1])
    #else:
    #    st = special[0]
    #    g.add_edge(special[1],special[0])
    #    
    #print st
    #
    #for edge in eulerian_circuit(g,st):
    #    edges.append(edge)
    #for edge in edges:
    #    f.write(edge[0])
    #    f.write('->')
    #f.write(edges[-1][1])
    #f.close()
    
    #challenge 4
    #with open('test.txt','r') as f:
    #    dnas = [item.strip() for item in f.readlines()]
    #print '\n'.join(debrujin_patterns(dnas))
    #f.close()
    #with open('test2.txt','w') as f:
    #    f.write('\n'.join(debrujin_patterns(dnas)))
    #f.close()
    #
    #with open('test2.txt','r') as f:
    #    lines = [item.strip() for item in f.readlines()]
    #f.close()
    #g = nx.DiGraph()
    #for line in lines:
    #    st = line.split(' -> ')
    #    for item in st[1].split(','):
    #        g.add_edge(st[0], item)
    #        print st[0],item
    #f = open('res.txt','w')       
    #edges = []
    #deg = g.degree(nx.nodes(g))
    #special = [item for item in  deg.keys() if deg[item]%2 == 1]
    #print special
    #print 'out',g.out_degree('AAG')
    #if g.in_degree(special[0]) >= g.in_degree(special[1]):
    #    st = special[1]
    #    g.add_edge(special[0],special[1])
    #else:
    #    st = special[0]
    #    g.add_edge(special[1],special[0])
    #    
    #print st
    #
    #for edge in eulerian_circuit(g,st):
    #    edges.append(edge)
    #for edge in edges[0:-1]:
    #    f.write(edge[0][0])
    #f.write(edge[1])
    #f.close()
    
    # challenge 5
    #k = 8
    #with open('test.txt','w') as f:
    #    for i in range(2**k):
    #        r  = baseN(i,2)
    #        for i in range(k-len(r)):
    #            r = '0'+r 
    #        f.write(r)
    #        f.write('\n')
    #f.close()
    #
    #
    #with open('test.txt','r') as f:
    #    dnas = [item.strip() for item in f.readlines()]
    #print '\n'.join(debrujin_patterns(dnas))
    #f.close()
    #with open('test2.txt','w') as f:
    #    f.write('\n'.join(debrujin_patterns(dnas)))
    #f.close()
    #
    #with open('test2.txt','r') as f:
    #    lines = [item.strip() for item in f.readlines()]
    #f.close()
    #g = nx.DiGraph()
    #for line in lines:
    #    print line
    #    st = line.split(' -> ')
    #    for item in st[1].split(','):
    #        g.add_edge(st[0], item)
    #f = open('res.txt','w')       
    #
    #edges = []
    #for edge in eulerian_circuit(g):
    #    edges.append(edge)
    #
    #edges2 = []
    #for edge in edges:
    #    edges2.append(edge[0]+edge[1][-1])
    #for edge in edges2:
    #    f.write(edge[0])
    #    
    
    #with open('test.txt','r') as f:
    #    k,d = map(int,f.readline().strip().split(' '))
    #    dnas = [line.strip().split('|') for line in f.readlines()]
    #f.close()
    #dnas = [(dna[0],dna[1]) for dna in dnas]
    #
    #graph = debrujin_patterns_pairs(dnas)
    #g = nx.DiGraph()
    #for edge in graph:
    #    src = edge[0]
    #    for dst in edge[1]:
    #        g.add_edge('|'.join(src), '|'.join(dst))
    #print nx.nodes(g)
    #print nx.edges(g)
    #
    #
    #deg = g.degree(nx.nodes(g))
    #special = [item for item in  deg.keys() if deg[item]%2 == 1]
    #print special
    #if g.in_degree(special[0]) >= g.in_degree(special[1]):
    #    st = special[1]
    #    g.add_edge(special[0],special[1])
    #else:
    #    st = special[0]
    #    g.add_edge(special[1],special[0])
    #    
    #print st
    #edges = []
    #for edge in eulerian_circuit(g,st):
    #    edges.append(edge)
    #print edges
    #
    #edges1 = [(item[0].split('|')[0], item[1].split('|')[0]) for item in edges]    
    #print edges1
    #
    #edges2 = [(item[0].split('|')[1], item[1].split('|')[1]) for item in edges]    
    #print edges2
    #
    #str1 = ''
    #str2 = ''
    #f = open('res.txt','w')
    #for edge in edges1[0:-1]:
    #    str1 += (edge[0][0])
    #str1 += (edge[1])
    #
    #for edge in edges2[0:-1]:
    #    str2 += (edge[0][0])
    #str2 += (edge[1])
    #
    #st = str1 + str2[len(str2)/2:len(str2)]
    #print st
    #
    #f.close()