#!/usr/bin/env python3
import pandas as pd
import re
import sys, os
import argparse
#import seqio from biopython
from Bio import SeqIO
#Funtcion to test if file is not empty
def file_is_not_empty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0
#    return os.stat(path).st_size!=0

def main():
    parser=argparse.ArgumentParser()
    #parser.add_argument("--asdbtsv", required=True, metavar='FILE')
    parser.add_argument("--plasmids", required=True, metavar='FILE')
    parser.add_argument("--clusters", required=True, metavar='FILE')
    args=parser.parse_args()

    #read plasmid fasta sequences and store in dictionary
    plasmids={}
    for record in SeqIO.parse(args.plasmids, "fasta"):
        plasmids[record.id]=record.seq

    #check whether there exist repeated bgcs contigs
    #read clusters from mcl
    repeated=set()
    representatives=set()
    if(file_is_not_empty(args.clusters)):
        clusters=[]
        with open(args.clusters) as f:
            #skip first line
            next(f)
            for line in f:                
                clusters.append(line.strip().split('\t'))
        #for each word in line add to repeated set
        for cluster in clusters:
            for word in cluster:
                repeated.add(word)
        #add keys in plasmids id not in repeated set
        for key in plasmids.keys():
            if key not in repeated:
                representatives.add(key)
        #for each cluster get sequence lenght and set as representative
        #the longest sequence
        for cluster in clusters:
            max_len=0
            max_key=''
            for word in cluster:
                if len(plasmids[word])>max_len:
                    max_len=len(plasmids[word])
                    max_key=word
            representatives.add(max_key)
     #print fasta file of all items in representatives set
        with open('representatives.fasta', 'w') as f:
            for key in representatives:
                f.write('>'+key+'\n'+str(plasmids[key])+'\n')
    

    #if cluster is empty then save all plasmids as representatives as fasta file
    else:
        with open('representative_plasmids.fna', 'w') as f:
            for key, value in plasmids.items():
                f.write('>'+key+'\n'+str(value)+'\n')
        
    #read paired distances <0.05

if __name__ == "__main__":
    main()
