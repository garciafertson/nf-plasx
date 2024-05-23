#!/usr/bin/env python3

#This script get the the contigs sequences and plamids scores
#returns the fasta sequences of the plasmids

import pandas as pd
import re
import sys, os
import argparse
#import biopython
from Bio import SeqIO

#read arguments
parser=argparse.ArgumentParser()
parser.add_argument("--plasmids", required=True, metavar='FILE')
parser.add_argument("--contigs", required=True, metavar='FILE')
parser.add_argument("--output", required=True, metavar='FILE')
args=parser.parse_args()

#read plasmid scores, the file has two columns, contig and score
df=pd.read_csv(args.plasmids, sep="\t", header=0)
#drop rows where contig column equals to 'contig'
df=df.loc[df["contig"]!="contig",:]
#read contigs
seqs=list(SeqIO.parse(args.contigs, "fasta"))
#compare sequence ids in seqs with contig column in df and keep
#the sequences that are in both
plasmids=[s for s in seqs if s.id in df["contig"].values]
#write output
SeqIO.write(plasmids, args.output, "fasta")

