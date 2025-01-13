#! /usr/bin/env python3

import os
import sys
import re
import seaborn as sns 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

### filter all reads containing intact recorder sequences (edited or not)
total = []
reads_rec = []
reads_bad = []
with open ("YOUR_FILE.txt", "r") as recorder_seq:
    for line in recorder_seq:
        line = line.rstrip()
        total.append(line)
        found = re.search(r"TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TCC[A/G][A/G]TTAAATT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G](?=C)", line)
        if found:
            reads_rec.append(line)
        else:
            reads_bad.append(line)
all_reads = len(total)
good = len(reads_rec)
bad = len(reads_bad)
#print(reads_empty)

### filter all reads that contain 7N insertion
seqs_all = []
motifs = []
recorders = []
total_edits = []

for sequence in reads_rec:
    if 'N' in sequence:
        continue  # Skip this iteration if "N" is present

    found = re.search(r"^GCCAACGCTATTCTGGCTGA(.{7})ATCTATCCAACGCAATTCTGGCACG", sequence)  ### MODIFY THIS CRITERIA FOR DIFFERENT TYPE OF INSERTION. 
    if found:
        seqs_all.append(sequence)
        motif = found.group(1)  # Extract the matched motif
        motifs.append(motif)
    else:
        foundE = re.search(r"^GCCAACGCTATTCTGGCTGAATTCTATCCAACGCAATTCTGGCACG", sequence)
        if foundE:
            seqs_all.append(sequence)
            motifE = "EMPTY"  # Assign "EMPTY" as the motif
            motifs.append(motifE)

### coubt all editing events in the recorder
for sequence in seqs_all:
    found1 = re.search("TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G]TT[A/G]G[A/G](?=C)", sequence)
    if found1:
        recorder = found1.group()
        recorders.append(recorder)

        start = found1.start()
        end = found1.end()
        region = sequence[start:end]

        # Count occurrences of "TTGGA" in the extracted region
        count1 = region.count("TTGGA")
        #print(count1)
        total_edits.append(count1)
#print(total_edits)

### Count all edited reads. Create a list with either 1 or 0 values. 1 is assigned if more then 1 edit is present in a read

molecules = []
for x in total_edits:
    if x >= 1:
        molecules.append(1)
    else:
        molecules.append(0)

### Count abundance of each motif and make a dictinary 
occurrence = {}
for item in motifs:
    if item in occurrence:
        occurrence[item] +=1
    else:
        occurrence[item] =1
#occurrence_sorted = sorted(occurrence) 
#print(occurence)



### Create a dictonary of all motifs and correspond reads that have at least one edit in one molecule 

dict_0 = {}
s = iter(motifs)
p = iter(molecules)
m = list(zip(s,p))
for (s,p) in m:
    if s in dict_0:
        dict_0[s] = dict_0[s]+p
    else:
        dict_0[s] = p
#print(dict_0)

		
### Build a dictonary of all motifs and total number of editing events for each ot them 

dict_1 = {}
i = iter(motifs)
j = iter(total_edits)
k = list(zip(i,j))
for (x,y) in k:
    if x in dict_1:
        dict_1[x] = dict_1[x]+y
    else:
        dict_1[x] = y
#print(dict_1)


### Building combined dictonary:

ds = [occurrence, dict_0, dict_1]
d = {}
for k in occurrence.keys():
    d[k] = tuple(d[k] for d in ds)
#print(d)


### Convert dictonary in pands dataframe

data = pd.DataFrame.from_dict(d, orient='index')
data.columns = ['occurrence','reads_edited','total_edits']
number_of_good_reads = len(seqs_all)
data['all_reads'] = number_of_good_reads
data['expression_percent']=(data['occurrence']/data['all_reads'])*100
data['efficiency_sq'] =(data['reads_edited']/data['occurrence'])*(data['total_edits']/data['occurrence'])
data['percent_edited']=(data['total_edits']/(data['occurrence']*6))*100

### Filter df by read occurence 
data_filtered = data[data['occurrence'] >= 10]  ### MODIFY THIS PARAMETER FOR DIFFERENT CUTOFF



### Sort the DataFrame by percetn editing in descending order
final = data_filtered.sort_values(by=['percent_edited'], ascending=False)
final = final.rename_axis("Motif").reset_index()
final = final[(final['Motif'].str.len() >2) | (final['Motif'] == "EMPTY")]
print(final.to_string())

final.to_csv("YOUR_DATA.csv")

