import sys
import string
import random
import os

def EditDistance(str1, str2, indel, sub):
    p = len(str1)
    q = len(str2)
    if (p==0) | (q==0):
        return p+q
    M = [[0]*(q+1)for _ in range (p+1)]
    for i in range(p+1):
        M[i][0] = i #???
    for j in range(q+1):
        M[0][j] = j #???
    for i in range(1,p+1):
        for j in range(1, q+1):
            if str1[i] == str2[j]:
                toAdd = 0
            else:
                toAdd = sub
            M[i][j] = min(M[i-1][j]+indel, M[i][j-1]+indel, M[i-1][j-1]+toAdd)
    return M[p][q]

def neighbour(pattern, mismatch, words,alphabet,indel,sub): #mismatch is D
    if (sub>2*indel):
        v=2*indel
    else:
        v=sub
    if mismatch == 0:
        words.add(pattern)
    else:
        for i in range(len(pattern)):
            for j in range(len(alphabet)):
                new_pattern = pattern[:i] + alphabet[j] + pattern[i+1:] #substituition
                #print(pattern)
                #print(new_pattern)
                #takes the first one with D edit distance obtained for a subseq being considered- greedy
                if mismatch <= 1:
                    words.add(new_pattern)
                else:
                    #Since it is always substitution, can we consider this as deletion + insertion
                    #So can I take the largest cost as the val here
                    neighbour(new_pattern, mismatch-v, words,alphabet,indel,sub)

def DiscoverMotifs_EditDist(S, L, D, indel, sub, n, m, t, alphabet):
    patterns = []
    #For each string in S
    for n1 in range(n):
        pattern = set()
        # allows for m (str length) to vary with each string
        #For all L sized substrings
        for i in range(len(S[n1]) - L + 1):
            Lmerspattern = set()
            neighbour(S[n1][i:i + L], D, Lmerspattern, alphabet,indel,sub)
            for words in Lmerspattern:
                #pattern contains all potential Motifs in one string
                pattern.add(words)
        for j in pattern:
            #patterns contains all potential Motifs in all the strings
            patterns.append(j)
    PatternSet = []
    for element in patterns:
        #For t= n-1, we have to get n ocurrences of the pattern
        # This can be made more general by choosing RHS as t+1
        if patterns.count(element) == t+1: # n: 
            PatternSet.append(element)
    PatternSet = list(set(PatternSet))
    return PatternSet

def DiscoverMotifs_neighbor(S, L, D, indel, sub, n, m, t, alphabet):
    patterns = []
    #For each string in S
    for n1 in range(n):
        pattern = set()
        # allows for m (str length) to vary with each string
        #For all L sized substrings
        for i in range(len(S[n1]) - L + 1):
            Lmerspattern = set()
            neighbour(S[n1][i:i + L], D, Lmerspattern, alphabet,indel,sub)
            for words in Lmerspattern:
                #pattern contains all potential Motifs in one string
                pattern.add(words)
        for j in pattern:
            #patterns contains all potential Motifs in all the strings
            patterns.append(j)
    PatternSet = []
    for element in patterns:
        #For t= n-1, we have to get n ocurrences of the pattern
        # This can be made more general by choosing RHS as t+1
        if patterns.count(element) == t+1: # n: 
            PatternSet.append(element)
    PatternSet = list(set(PatternSet))
    return PatternSet

def readFile(filename):
    with open(filename, "r") as f:
        #n: number of strings in DNA database (S)
        #m: length of each string
        #L: length of substring/pattern/motif
        #D: maximum edit distance
        #t: number of strings in which pattern exists at a max distance D = n-1 for simplifies algorithm
        #indel: cost of insertion/deletion
        #sub: cost of substituition
        L, D, indel, sub, n, m, t = map(int, f.readline().strip().split())
        S = [line.strip() for line in f]        

    return [S, L, D, indel, sub, n, m, t]

def genRndFile(rndFilename, alphabet,L, D, indel, sub, n, m, t):
    #randomly generate 20 strings of size 600 through a function or read a file
    with open(rndFilename, 'w') as f:
        s= str(L)+" "+str(D)+" "+str(indel)+" "+str(sub)+" "+str(n)+" "+str(m)+" "+str(t)+"\n"
        f.write(s)
        for i in range(n):
            res = ''.join(random.choices(alphabet, k = m))+"\n"
            #print(res)
            f.write(res)

def display_menu():
    print("\nMotif Search Menu")  
    print("1. Based on DP Edit Distance")  
    print("2. Based on Neighbor")  
    opt = int(input("Enter the Choice:"))  
    return opt

def writeOutFile(outFilename,PatternSet):
    with open(outFilename, 'w') as f:
        for element in PatternSet:
            f.write(element+"\n")

if __name__ == "__main__":

    #To run enter python Motif.py
    #To run enter python Motif.py r for first generating random data

    alphabet = ['A', 'C', 'G', 'T']
    filename = "Data.txt"
    outFilename="Out.txt"
    n=20 #value given in assignment
    t=19 #value given in assignment = n-1
    m=600 #value given in assignment    
    L=15 #15 #value given in assignment    
    D=5 #5 #value given in assignment    
    sub=1 
    indel=1 

    if len(sys.argv)>1: 
        if sys.argv[1]=="r": #random file to be generated
            S=genRndFile(filename, alphabet,L, D, indel, sub, n, m, t)
    
    [S, L, D, indel, sub, n, m, t] = readFile(filename)
    
    #n,m can be computed from S also
    #For our purpose t=n-1
    #Setting these values here whether given or not in the data
    n=len(S)
    m=len(S[1]) #Assumed that each string is of same length
    t=n-1
    
    opt=display_menu()
    if opt == 1:  
        PatternSet = DiscoverMotifs_EditDist(S, L, D, indel, sub, n, m, t, alphabet)          
    elif opt == 2:  
        PatternSet = DiscoverMotifs_neighbor(S, L, D, indel, sub, n, m, t, alphabet)                                    
    else:  
        print("Incorrect choice");
        exit()
    #remove duplicates
    PatternSet = list(set(PatternSet))
    writeOutFile(outFilename,PatternSet)
    for i in PatternSet:
        print(i)
