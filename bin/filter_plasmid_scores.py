#!/usr/bin/env python3
# This script reads the plasx score and only keeps the plasmids with a score above a the input threshold

import pandas as pd
import re
import sys, os
import argparse

#read arguments
parser=argparse.ArgumentParser()
parser.add_argument("--plasmids", required=True, metavar='FILE')
parser.add_argument("--threshold", required=True, metavar='FLOAT')
parser.add_argument("--output", required=True, metavar='FILE')
args=parser.parse_args()

#read plasmid scores
df=pd.read_csv(args.plasmids, sep="\t", header=None)
df.columns=["plasmid", "score"]
df=df[df["score"]>float(args.threshold)]
df.to_csv( args.output, sep="\t", header=False, index=False)