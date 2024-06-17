#!/usr/bin/env python3
# Get contig file from anvio hmm fasta file headers
# and build table mge information

# import libraries
import pandas as pd
import re
import argparse
import os
# import the fasta parser
from Bio import SeqIO

def parse_header(header):
    #fasta header fields are divided by "|" character
    #remove ">" character from the beginning of the string
    header = header[1:]
    fields = header.split("|")
    #fields are divided by ":" character
    evalue= fields[2].split(":")[1]
    contig = fields[3].split(":")[1]
    gene= fields[4].split(":")[1]
    start= fields[5].split(":")[1]
    stop= fields[6].split(":")[1]
    length= fields[7].split(":")[1]
    #hmm name is in in the first field before "___"
    hmm_name = fields[0].split("___")[0]
    #return string separated by "\t" character
    return "\t".join([contig, gene, start, stop, length, evalue, hmm_name])


#def main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hmm", required=True, metavar='FILE')
    parser.add_argument("--output", required=True, metavar='FILE')
    parser.add_argument("--contigs", required=True, metavar='FILE')
    args = parser.parse_args()

    # read contigs fasta file and get fasta headers
    contigs = set()
    first = True
    with open(args.output + "_rechmm.fna", "w") as out, open(args.output + "_rechmm.tsv", "w") as table_out:
        table_out.write("\t".join(["contig", "gene_callers_id", "start", "stop", "length", "e_value", "hmm_name"]) + "\n")
        for record in SeqIO.parse(args.hmm, "fasta"):
            #divide fasta header by "|" character
            header = record.description.split("|")
            # get contig name from foutrh element
            contig = header[3].split(":")[1]
            contigs.add(contig)
            gene=header[3] + "|" + header[4]
            if first:
                #print record to output file
                SeqIO.write(record, out, "fasta")
                #print table to table output file
                table_out.write(parse_header(record.description) + "\n")
                prevgene=gene
                first = False
            else:
                print(gene, prevgene)
                if gene != prevgene:
                    #print record to output file
                    SeqIO.write(record, out, "fasta")
                    prevgene=gene
                    table_out.write(parse_header(record.description) + "\n")

    # read contigs fasta file and return conitgs in conitg set
    with open(args.output+ "_recContig.fna" , "w") as out:
        for record in SeqIO.parse(args.contigs, "fasta"):
            if record.id in contigs:
                SeqIO.write(record, out, "fasta")
if __name__ == '__main__':
    main()
