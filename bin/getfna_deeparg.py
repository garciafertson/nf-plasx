#!/usr/bin/env python3
# This script reads the deepARG prediction file, and contigs file 
# and returns sequence of predicted ARGs

#import libraries
import pandas as pd
import re
import argparse
#import the fasta parser
from Bio import SeqIO

#read arguments
parser=argparse.ArgumentParser()
parser.add_argument("--deeparg", required=True, metavar='FILE')
parser.add_argument("--orfs", required=True, metavar='FILE')
parser.add_argument("--output", required=True, metavar='FILE')

#read deeparg predictions
#Column 4 has te contig id
arg_id=set()
df=pd.read_csv(args.deeparg, sep="\t", header=0)
for index, row in df.iterrows():
    arg_id.add(row[3])

#read contigs fasta file line by line and save fasta sequence where is is in arg_id
with open(args.output, "w") as out:
    for record in SeqIO.parse(args.orfs, "fasta"):
        if record.id in arg_id:
            out.write(">"+record.id+"\n"+str(record.seq)+"\n")
