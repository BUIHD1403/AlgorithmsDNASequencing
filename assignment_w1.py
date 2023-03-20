#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().system('wget https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa')


# In[3]:


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# In[4]:


seqs, quals = readFastq('lambda_virus.fa')


# In[5]:


print(seqs[:5])


# In[6]:


def phred33ToQ(qual):
    return ord(qual) - 33


# In[7]:


phred33ToQ('#')


# In[51]:


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


# In[52]:


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


# In[53]:


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


# In[54]:


t = 'AGGT'


# In[55]:


genome = readGenome('lambda_virus.fa')


# In[56]:


occurence = naive(t, genome)


# In[57]:


len(occurence)


# In[58]:


occurence_reverse = naive(reverseComplement(t), genome)


# In[59]:


len(occurence_reverse)


# In[60]:


count = len(occurence) + len(occurence_reverse)


# In[61]:


count


# In[62]:


count_q2 = len(naive('TTAA', genome))


# In[63]:


count_q2


# In[64]:


genome.find('AGTCGA')


# In[65]:


def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mist_match = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mist_match += 1
                if mist_match >2:
                    break
        if mist_match <= 2:
            occurrences.append(i)  # all chars matched; record
    return occurrences


# In[66]:


naive_2mm('ACTTTA', 'ACTTACTTGATAAAGT')


# In[67]:


count = len(naive_2mm('TTCAAGCC', genome))
print(count)


# In[69]:


offsets = naive_2mm('AGGAGGTT', genome)
print(offsets)


# In[70]:


get_ipython().system('wget https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq')


# In[71]:


def read_FAST_Q(filename):
    sequences = []
    qualities = []
    with open(filename, 'r') as f:
        while True:
            f.readline()
            seq = f.readline().rstrip()
            f.readline()
            qual = f.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


# In[72]:


seq, qual = read_FAST_Q('ERR037900_1.first1000.fastq')
print(seq[:5])
print(qual[:5])
print(len(seq[0]))


# In[73]:


def createHistory(qualities):
    history = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            history[q] += 1
    return history
h = createHistory(qual)
print(h)


# In[74]:


get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.bar(range(len(h)), h)
plt.show()


# In[75]:


import collections

def maxPoorQualitySequencingCycle(qualities):
    min_score = 123456789
    min_index = -1
    for i, qual in enumerate(qualities):
        score = sum(map(ord, qual))
        if min_score > score:
            min_score = score
            min_index = i
    return min_index


# In[76]:


offset = maxPoorQualitySequencingCycle(qual)
print("incorrect answer: %d" % offset)


# In[77]:


print("This my guess...\n" + qual[111])


# In[ ]:




