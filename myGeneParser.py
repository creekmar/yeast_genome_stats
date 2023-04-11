#!/usr/bin/env python
"""
Adapted from Gregory A Babbitt's myGeneParser.py from code-course-repo
Reads in YeastGene sequences from a folder and calculates gc frequency and gets
log2Transcription values from Yeast_RNAseq folder
Wrties this content into two files log2Transcription.txt and gcContent.txt

These two files are then analyzed through an Rscript yeast_gene_stats.r
Specifically, it makes a bar chart of the means and standard error,
a linear regression comparing gcContent means vs. log2Transcription means,
and runs ANOVA test

@authors: Gregory A Babbitt
          Ming Creekmore
"""

import fileinput
import os
import os.path
from os import path
import glob

#Reads in seqeunces as elements in list
# header holds the filename in position corresponding to position in sequence
header = list()
sequence = list()
dna_to_protein = {"UUU": "PHE", "UUC": "PHE", "UUA": "LEU", "UUG": "LEU", "UAU": "TYR", "UAC": "TYR",
                  "UCU": "SER", "UCC": "SER", "UCA": "SER", "UCG": "SER", "AGU": "SER", "AGC": "SER",
                  "UGA": "STOP", "UGG": "TRP", "CUU": "LEU", "CUC": "LEU", "CUA": "LEU", "CUG": "LEU",
                  "CCU": "PRO", "CCC": "PRO", "CCA": "PRO", "CCG": "PRO", "CAU": "HIS", "CAC": "HIS",
                  "AGA": "ARG", "AGG": "ARG", "CGU": "ARG", "CGC": "ARG", "CGA": "ARG", "CGG": "ARG",
                  "AUU": "ILE", "AUC": "ILE", "AUA": "ILE", "AUG": "MET", "AAU": "ASN", "AAC": "ASN",
                  "ACU": "THR", "ACC": "THR", "ACA": "THR", "ACG": "THR", "AAA": "LYS", "AAG": "LYS",
                  "GUU": "VAL", "GUC": "VAL", "GUA": "VAL", "GUG": "VAL", "UAA": "STOP", "UAG": "STOP",
                  "CAA": "GLN", "CAG": "GLN", "GCU": "ALA", "GCC": "ALA", "GCA": "ALA", "GCG": "ALA",
                  "GAU": "ASP", "GAC": "ASP", "GGU": "GLY", "GGC": "GLY", "GGA": "GLY", "GGG": "GLY",
                  "GAA": "GLU", "GAG": "GLU", "UGU": "CYS", "UGC": "CYS"}

# main program
def main():
    tot_file = 0
    path_foldername = os.getcwd()
    foldername = 'YeastGenes'
    files = os.listdir(foldername)
    files.sort()

    # reading yeast genes and putting them in lists
    for filename in files:
        tot_file += 1
        temp = ""
        my_path = path.join(path_foldername,foldername, filename)
        with open(my_path) as file:
            file.readline()
            temp += file.readline()
        filename = filename[:-4]
        header.append(filename)
        sequence.append(temp)
    average = 0

    # analyzing yeast genes and putting them in files
    gene = header[0][:2]
    rna_seq = ""
    with open("log2Transcription.txt", "w") as logFile, open("protein_sequences.txt", "w") as file,\
            open("gcContent.txt", "w") as gcFile:
        logFile.write(gene)
        gcFile.write(gene)
        for i in range(len(sequence)):

            # finding trans_level test
            trans_test = "NA"
            with open("Yeast_RNAseq/Nagalakshmi_2008_5UTRs_V64.gff3", "r") as rna_datafile:
                for line in rna_datafile:
                    if line.find(header[i]) != -1:
                        line = line.split()
                        line_seg = line[8]
                        trans_test = line_seg.split(";")[2][25:25+5]

            # combining gene sequences of same yeast gene and transcribing
            transcribedSeq = transcribe(sequence[i], header[i][-1])
            gc = gcContent(transcribedSeq)
            if header[i][:2] != gene or i == len(sequence)-1:
                aa_seq = full_translation(rna_seq)

                # writing data to file
                file.write(gene + "protein:\t")
                for protein in aa_seq:
                    file.write(protein + "\t")
                file.write("\n")
                if i != len(sequence)-1:
                    rna_seq = transcribedSeq
                    gene = header[i][:2]
                    logFile.write("\n" + gene + "\t" + trans_test)
                    gcFile.write("\n" + gene + "\t" + str(gc))
                else:
                    logFile.write("\t" + trans_test)
                    gcFile.write("\t" + str(gc))
            else:
                rna_seq += transcribedSeq
                logFile.write("\t" + trans_test)
                gcFile.write("\t" + str(gc))


def gcContent(seq):
    """
    Finds the frequency of GC nucleotides in a given mRNA sequence
    :param seq: the mRNA sequence
    :return: GC frequency as a float
    """
    count = 0
    tot = 0
    final = 0
    for i in seq:
        if (i == "C") or (i == "G"):
            count = count + 1
        elif (i == "A") or (i == "U"):
            tot = tot + 1
    tot = tot + count
    if tot != 0:
        final = (count/tot) * 100
    return final


def transcribe(seq, watsonCrick):
    """
    Given a DNA sequence, transcribes it into the mRNA sequence
    :param seq: DNA sequence as a string
    :param watsonCrick: the file type
    :return: mRNA sequence as a string
    """
    if watsonCrick == "W":
        return seq.replace("T", "U")
    seq = seq[::-1]
    for i in range(len(seq)):
        if seq[i] == "A":
            seq = seq[:i] + "U" + seq[i+1:]
        elif seq[i] == "T":
            seq = seq[:i] + "A" + seq[i+1:]
        elif seq[i] == "G":
            seq = seq[:i] + "C" + seq[i+1:]
        elif seq[i] == "C":
            seq = seq[:i] + "G" + seq[i+1:]
    return seq


def find_possible_start(seq):
    """
    Given an rna sequence, will return a list of all the possible start
    positions (where AUG is read)
    :param seq: the rna sequence to translate
    :return: list that represents possible start areas
    """
    start_pos = list()
    for i in range(0, len(seq)-2):
        codon = seq[i:i+3]
        if "AUG" == codon:
            start_pos.append(i)
    return start_pos


def translation(seq, index):
    """
    Given a rna sequence and an index to start at, will translate the rna to
    amino acid sequence until reaches a stopping point
    :param seq: the rna sequence to translate
    :param index: the start index to start translating
    :return: a list of amino acids
    """
    protein_seq = list()
    aa = ""
    while index < (len(seq)-2) and aa != "STOP":
        codon = seq[index:index+3]
        if codon in dna_to_protein:
            aa = dna_to_protein[codon]
            protein_seq.append(aa)
        else:
            protein_seq.append("---")
        index += 3
    return protein_seq


def get_protein_seq(rna_seq):
    """
    Given an rna_seq, will find the longest protein translation
    Starting with AUG and ending on a stop codon
    :param rna_seq: the rna_seq to be translated
    :return: list of amino acid
    """
    start_pos = find_possible_start(rna_seq)
    aa_seq = []

    # finding the longest translation
    if len(start_pos) > 0:
        for start in start_pos:
            temp = translation(rna_seq, start)
            if len(temp) > len(aa_seq):
                aa_seq = temp
    else:
        aa_seq = ["NONE"]
    return aa_seq


def full_translation(seq):
    """
    Given a rna sequence will translate the rna to amino acid sequence
    :param seq: the rna sequence to translate
    :return: a list of amino acids
    """
    index = 0
    protein_seq = list()
    aa = ""
    while index < (len(seq)-2):
        codon = seq[index:index+3]
        if codon in dna_to_protein:
            aa = dna_to_protein[codon]
            protein_seq.append(aa)
        else:
            protein_seq.append("---")
        index += 3
    return protein_seq


if __name__ == "__main__":
    main()
    os.system("Rscript yeast_gene_stats.r")
    #end script
    print ("\nend myGeneParser.py")
    exit