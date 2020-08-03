"""
Created on Feb  7, 2019
@author: blake.inderski
"""

import os
import math
import re
import random
import string

import subprocess
import numpy as np
#import statistics as stats
import sklearn.cluster
import json
from Bio.SubsMat import MatrixInfo
import shutil
import sys
import traceback

#from time import time


aa_convert_codon_di =  {
    'A':['[GRSK][CYSM].'],
    'B':['[ARWM][ARWM][CTYWKSM]', '[GRSK][ARWM][TCYWKSM]'],
    'C':['[TYWK][GRSK][TCYWKSM]'],
    'D':['[GRSK][ARWM][TCYWKSM]'],
    'E':['[GRSK][ARWM][AGRSKWM]'],
    'F':['[TYWK][TYWK][CTYWKSM]'],
    'G':['[GRSK][GRSK].'],
    'H':['[CYSM][ARWM][TCYWKSM]'],
    'I':['[ARWM][TYWK][^G]'],
    'J':['[ARWM][TYWK][^G]', '[CYSM][TYWK].', '[TYWK][TYWK][AGRSKWM]'],
    'K':['[ARWM][ARWM][AGRSKWM]'],
    'L':['[CYSM][TYWK].', '[TYWK][TYWK][AGRSKWM]'],
    'M':['[ARWM][TYWK][GRSK]'],
    'N':['[ARWM][ARWM][CTYWKSM]'],
    'O':['[TYWK][ARWM][GRSK]'],
    'P':['[CYSM][CYSM].'],
    'Q':['[CYSM][ARWM][AGRSKWM]'],
    'R':['[CYSM][GRSK].', '[ARWM][GRSK][GARSKWM]'],
    'S':['[TYWK][CYSM].', '[ARWM][GRSK][CTYWKSM]'],
    'T':['[ARWM][CYSM].'],
    'U':['[TYWK][GRSK][ARWM]'],
    'V':['[GRSK][TYWK].'],
    'W':['[TYWK][GRSK][GRSK]'],
    'X':['...'],
    'Y':['[TYWK][ARWM][CTYWKSM]'],
    'Z':['[CYSM][ARWM][AGRSKWM]','[GRSK][ARWM][AGRSKWM]'],
    '_':['[TYWK][ARWM][AGRSKWM]', '[TYWK][GRSK][ARWM]'],
    '*':['[TYWK][ARWM][AGRSKWM]', '[TYWK][GRSK][ARWM]'],
    'x':['[TYWK][ARWM][AGRSKWM]', '[TYWK][GRSK][ARWM]']}

#Stop codons are currently represented by 'X' (ambig residue) to prevent MAFFT auto deletion
dna_convert_aa_di = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'x', 'TAG':'x',
    'TGC':'C', 'TGT':'C', 'TGA':'x', 'TGG':'W'}

aa_convert_dna_di = {dna_convert_aa_di[key]:[x for x in dna_convert_aa_di if dna_convert_aa_di[x] == dna_convert_aa_di[key]] for key in dna_convert_aa_di}

ambiguous_dna_di = {
    'Y':['C', 'T'], 'R':['A', 'G'], 'W':['A', 'T'],
    'S':['G', 'C'], 'K':['T', 'G'], 'M':['C', 'A'],
    'D':['A', 'G', 'T'], 'V':['A', 'C', 'G'], 'H':['A', 'C', 'T'],
    'B':['C', 'G', 'T'], 'N':['A', 'C', 'T', 'G']}


aa_li = [x for x in "ARNDCEQGHILKMNFPSTWYV"]
ambiguous_aa_di = {"B":["N", "D"], "Z":["E", "Q"], "J":["I", "L"], "X":aa_li}
blosum_di = MatrixInfo.blosum62

dir_path = os.path.dirname(os.path.realpath(__file__))
py_file = os.path.basename(__file__)
#when debugging, set environment path to match bash
#run the following code in unix shell: echo $PATH
if len(sys.argv) == 1:
    os.environ["PATH"] = "/Users/blake.inderski/anaconda/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/sbin" #copy/paste terminal output here
valid_organism_li = [d for d in os.listdir(dir_path+"/required/") if os.path.isdir(dir_path+"/required/"+d) and d != "blastdb"]

################################################################################
################################################################################

def writeFASTA(output_filepath, read_fasta_li2):
    #create missing directories if path does not exist (requires bash)
    subprocess.call("mkdir -p "+re.sub("[^\/]*$", "", output_filepath), shell=True)
    """
    #may be useful if mkdir cannot be called
    dir_li = output_filepath.split("/")[:-1]
    dir_li = [("/").join(dir_li[:i+1]) for i, x in enumerate(dir_li)]
    miss_dir_li = [x for x in dir_li if x != "" and not os.path.isdir(x)]
    """
    write_file = open(output_filepath, "w")
    for x in range(len(read_fasta_li2[0])):
        write_file.write(">" + read_fasta_li2[0][x].replace(" ", "_") + "\n" + read_fasta_li2[1][x] + "\n")
    write_file.close()

################################################################################

def readFASTA(input_filepath):
    #supports mixed case (upper/lower); can be used to indicate residues of interest
    #example: "x" is used as (unofficial) notation for stop codons in aminio acid sequences while "X" is for ambiguous residues
    #altered notation prevents errors when using mafft, which does not support standard notation of stop codons ("*")
    read_fasta_li2 = [[], []]
    if os.path.isfile(input_filepath):
        read_file = open(input_filepath, "r")
        read_file_str = "".join([line for line in read_file])
        read_file.close()
        fasta_match = [list(x) for x in re.findall("(>.*\n)([^>]*)", read_file_str)]
        for index, match in enumerate(fasta_match):
            if match[1].replace("\n", "").islower(): #convert all characters to uppercase only if all characters are lowercase
                fasta_match[index][1] = match[1].upper()
        #convert RNA to DNA
        read_fasta_li2 = [[match[0].replace("\n", "").replace(">", "") for match in fasta_match], [match[1].replace("\n", "").replace("U", "T") for match in fasta_match]]
    return(read_fasta_li2)

################################################################################
    
def mafftAlign(input_fasta, output_fasta):
    #uses bash to call mafft to perform a residue (nucleotide or amino acid) alignment
    mafft_output = re.sub("\.[^\.]*?$", "-mafft.fasta", input_fasta) #allows output_fasta to be the same as input_fasta
    os_call = "mafft "+input_fasta+" > "+mafft_output
    subprocess.call(os_call, shell=True)
    read_msa_li2 = readFASTA(mafft_output)
    read_msa_li2[1] = [x.upper() for x in read_msa_li2[1]] #force mafft sequence output to be uppercase (expected mafft format)
    #MAFFT output is not case sensitive
    #AA translation function produces "x" as stop codon to avoid breaking scripts
    read_fasta_li2 = readFASTA(input_fasta)
    for msa_index, head in enumerate(read_msa_li2[0]):
        fasta_index = read_fasta_li2[0].index(head)
        for index1, character1 in enumerate(read_fasta_li2[1][fasta_index]):
            if character1.islower():
                for index2, charater2 in enumerate(read_msa_li2[1][msa_index][index1:]):
                    if len(read_msa_li2[1][msa_index][:index1+index2].replace("-", "")) == index1:
                        read_msa_li2[1][msa_index] = read_msa_li2[1][msa_index][:index1+index2] + \
                        read_msa_li2[1][msa_index][index1+index2].lower() + \
                        read_msa_li2[1][msa_index][index1+index2+1:]
    #MAFFT alignment can produce entirely gapped columns          
    read_msa_li2[1] = removeColumnGap(read_msa_li2[1], 1)[0]
    writeFASTA(output_fasta, read_msa_li2)
    os.remove(mafft_output)

################################################################################

def removeColumnGap(sequence_li, gap_proportion):    
    #remove columns based on allowed gap proportion; track alterations
    #mafft alignments may produce columns that contain all gaps
    #hmm profile alignments do not track excessively gapped (>50%) regions
    column_index = 0 #index conversion from input to output: index of input minus indexed value of modify_li equals the index of output
    modify_li = [0]
    #avoid overwriting input list
    read_fasta_sequences = sequence_li[:]
    while column_index < len(read_fasta_sequences[0]):
        column_li = [read_fasta_sequences[x][column_index] for x in range(len(read_fasta_sequences))]
        if column_li.count("-") >= len(column_li)*gap_proportion:
#           print("column removed", {x:column_li.count(x) for x in set(column_li)}, file=sys.stderr)
            for row_index in range(len(read_fasta_sequences)):
                read_fasta_sequences[row_index] = read_fasta_sequences[row_index][:column_index] + read_fasta_sequences[row_index][column_index+1:]
            modify_li[-1] += 1
        else:
            column_index += 1
        modify_li.append(modify_li[-1])
    return(read_fasta_sequences, modify_li[:-1])
        
################################################################################

def translateNucleotides(input_dna, reading_frame=0):
    #input_dna in string format, reading_frame enter 0, 1, or 2
    #translates nucleotides sequence into protein (stop codons = "x"); able to infer ambiguous resides
    reading_frame = reading_frame%3
    frame_input = input_dna[reading_frame:].upper()
    frame_input = frame_input.replace('U', 'T')
    frame_input = frame_input.replace('-', '')
    translate_str = ''
    if len(frame_input)%3 != 0:
        frame_input = frame_input[:-(len(frame_input)%3)]
    for i in range(0, len(frame_input), 3):
        if frame_input[i:i+3] in dna_convert_aa_di:
            translate_str += dna_convert_aa_di[frame_input[i:i+3]]
        else:
            ambiguous_combination_li = ['']
            for char in frame_input[i:i+3]:
                if char in ambiguous_dna_di:
                    ambiguous_combination_li = ambiguous_combination_li*len(ambiguous_dna_di[char])
                    z = 0
                    for y in range(len(ambiguous_dna_di[char])):
                        for x in range(len(ambiguous_combination_li)//len(ambiguous_dna_di[char])):
                            ambiguous_combination_li[z] += ambiguous_dna_di[char][y]
                            z += 1
                else:
                    for c in range(len(ambiguous_combination_li)):
                        ambiguous_combination_li[c] += char
            translate_check = dna_convert_aa_di[ambiguous_combination_li[0]]
            for codon in ambiguous_combination_li[1:]:
                if dna_convert_aa_di[codon] != translate_check:
                    translate_check = 'fail'
            if translate_check == 'fail':
                translate_str += 'X'
            else:
                translate_str += translate_check
    return(translate_str)
    
################################################################################
        
#https://stats.stackexchange.com/questions/123060/clustering-a-long-list-of-strings-words-into-similarity-groups
#https://www.programcreek.com/python/example/85778/sklearn.cluster.AffinityPropagation

def clusterExpression(sequence_li, cropped_msa_li):
    #iterative clustering of residues to produce selective (minimize incorrect matches) and specific (maximize correct matches) regular expressions
    #stepwise process of increasing selectivity whle decreasing specificity
    #input sequences are grouped by similarity; residues may be removed
    read_fasta_li = [x.replace("-", "") for x in sequence_li]    
    for crop_index, crop in enumerate(cropped_msa_li):
        crop = crop.replace("-", "") #short regex unlikely to pass duplicate check
        if len([x for x in crop if x not in ambiguous_aa_di]) >= min(10, len(cropped_msa_li[0])): #expression must be at least 10 non-ambiguous residues in length
            for read_index, read in enumerate(read_fasta_li):
                if read.count(crop) > 1: #target sequence likely contains repeat region
                    read_fasta_li[read_index] = "removed"
        else:
            cropped_msa_li[crop_index] = "removed"
    read_fasta_li = [x for x in read_fasta_li if x != "removed"]
    cropped_msa_li = [x for x in cropped_msa_li if x != "removed"]
    output = None
                
    #remove low occurrence residues prior to clustering
    #rare residues can be a result of mutation or poor alignment
    remove_threshold, remove_count = 0, 0
    word_li = list(set(cropped_msa_li))
    copy_word_li = word_li[:]
    if len(word_li) == 1:
        output = word_li[0], 0
    while len(word_li) > 1:
        if output != None and output[1] <= remove_count:
            return(output)
        distance_similarity = [[blosum_di[w1, w2] if (w1, w2) in blosum_di else \
                                [min([blosum_di[x] for x in blosum_di]) if w1 == w2 else max([blosum_di[x] for x in blosum_di])][0] \
                                for w1 in word_li] for w2 in word_li]
        
        max_n = min(10, len(word_li)-1)
        for n_cluster in range(1, max_n+1):
            #random error may occur during fit
            #grouping is random; difficult to optomize
            try:
                affprop = sklearn.cluster.SpectralClustering(n_clusters=n_cluster, affinity="precomputed").fit_predict(distance_similarity)
                segment_li2, cluster_li2 = [], []
                
                for cluster_index in list(set(affprop)):
                    cluster = [word_li[x] for x, y in enumerate(affprop) if y == cluster_index]
                    segment_li = ["".join(set([x[index] for x in cluster])) for index in range(len(cluster[0]))]
                    segment_li2.append(segment_li)
                    cluster_li2.append(cluster)
                    
                for index in range(len(segment_li2)):
                    product_li = []
                    for column in segment_li2[index]:
                        for character in ambiguous_aa_di:
                            column.replace(character, "".join(ambiguous_aa_di[character]))
                        if "-" in column:
                            product_li.append(1)
                        else:
                            product_li.append(len(column)/20)
                    if np.prod(product_li) > 10**-11:
                        segment_li2[index] = []
                    else:
                        cluster_li2[index] = []
                    #replace beginning positions containing "-" with "." to maintain start position
#                   segment_li2[index] = ["." if "-" in x and len([y for y in segment_li2[index][:i+1] if "-" in y]) == i+1 else x for i, x in enumerate(segment_li2[index])]
                    segment_li2[index] = ["["+x.replace("-", "")+"]?" if "-" in x else "["+x+"]" for x in segment_li2[index]]
                    segment_li2[index] = [re.sub("["+"".join([y for y in ambiguous_aa_di])+"]", "", x) for x in segment_li2[index]] #remove ambiguous characters
                    segment_li2[index] = [x.replace("[]?", "").replace("[]", ".") for x in segment_li2[index]] #remove "[]" (former ambiguous character) and "[]?" (former gap)
                    segment_li2[index] = [re.sub(r"\[([A-Z\.]{1})\]", r"\1", x) for x in segment_li2[index]] #remove brackets around single characters
                
                miss_count = 0
                regex = "|".join(["".join(x) for x in segment_li2 if x != []])
                for sequence in read_fasta_li:
                    #detect misalignment in protein MSA
                    if len(re.findall(regex, sequence)) > 1:
                        miss_count = len(read_fasta_li) #prevent output if multiple matches detected (reiterate)
                        break
                    elif len(re.findall(regex, sequence)) == 0:
                        miss_count += 1
                
                if miss_count == 0:
                    return(regex, miss_count)
                elif miss_count < len(read_fasta_li) and (output == None or miss_count < output[1]):
                    output = regex, miss_count
            #multiple consecutive failures unlikely
            except:
#               print(n_cluster, word_li, file=sys.stderr)
                pass
                    
        #rebuilds cropped sequences to reflect highest character (including gaps/ambiguities) homology
        while word_li == copy_word_li:
            cropped_column_li2 = [[cropped_msa_li[x][y] for x in range(len(cropped_msa_li))] for y in range(len(cropped_msa_li[0]))]
            for column_index in range(len(cropped_column_li2)):
                column_li = cropped_column_li2[column_index]
                max_residue = max(set(column_li), key = column_li.count)
                replace_li = []
                for residue in set(column_li):
                    if column_li.count(residue) == remove_threshold:
                        replace_li.append(residue)
                for row_index in range(len(cropped_column_li2[column_index])):
                    if cropped_column_li2[column_index][row_index] in replace_li:
                        split_row_li = [x for x in cropped_msa_li[row_index]]
                        split_row_li[column_index] = max_residue
                        cropped_msa_li[row_index] = "".join(split_row_li)
            remove_threshold += 1
            word_li = list(set(cropped_msa_li))
        copy_word_li = word_li[:]
        remove_count = len([1 for x in read_fasta_li if len(re.findall("|".join([y.replace("-", "") for y in word_li]), x)) == 0])
                
    return(output)
        
################################################################################  
    
#https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/shannon.htm

def calculateShannon(input_fasta):
    #calculates Shannon diversity index for each column in MSA fasta to gauge homology/conservation
    read_fasta_li = readFASTA(input_fasta)[1]
    shannon_li = []
    for column_index in range(len(read_fasta_li[0])):
        column = [read_fasta_li[x][column_index] for x in range(len(read_fasta_li))]
        gap_count = [read_fasta_li[x][column_index] for x in range(len(read_fasta_li))].count("-")
        column = [x for x in column if x != "-"] #remove gaps from target column
        proportion_li = [column.count(x)/len(column) for x in set(column)]
        if len(proportion_li) > 1:
            shannon_diversity = -sum([x*math.log(x) for x in proportion_li])
            shannon_equitability = shannon_diversity/(math.log(len(proportion_li)))
        else:
            shannon_equitability = 0
        shannon_li.append((1-shannon_equitability)*len(column)/(len(column)+gap_count))
    return(shannon_li)
    
################################################################################   
    
def regexSelect(protein_directory, cropped_regex_json):
    #from descending order of highest homology/conservation, regions are converted into regular expressions
    #regulate maximum distance/overlap between regions
    fasta_li = [f for f in os.listdir(protein_directory+"/msa") if os.path.isfile(os.path.join(protein_directory+"/msa", f)) and ".fa" in f]    
    for fasta in fasta_li:
        annotation = fasta.split(".")[0]
        msa_save_path = protein_directory.replace("protein", "")+"/hmm_profiles/msa_save/"+fasta
        if annotation not in cropped_regex_json:
            cropped_regex_json[annotation] = []
            protein_fasta = protein_directory+"/msa/"+fasta
            protein_fasta_li2 = readFASTA(protein_fasta)
            shannon_li = calculateShannon(protein_fasta)
            shannon_window_li = [sum(shannon_li[x:x+15]) for x in range(1, len(shannon_li)-16)] #allow near complete overlap with terminal ends
            shannon_window_li2 = [[x+1, y] for x, y in enumerate(shannon_window_li)] #index, value
            shannon_window_li2.sort(reverse = True, key = lambda x: x[1]) #sort by highest window sum
            #identify gaps removed from msa_save
            save_fasta_li2 = readFASTA(msa_save_path)
            modify_li = removeColumnGap(save_fasta_li2[1], 0.5)[1] #gap proportion should be the same as in checkMSA()
            
            index_li = [0, max(0, len(shannon_li)-15)]
            while len(shannon_window_li2) > 0:
                window_index = shannon_window_li2[0][0]
                min_less_int = min([x for x in index_li if x < window_index], key=lambda x:abs(x-window_index))
                min_more_int = min([x for x in index_li if x > window_index], key=lambda x:abs(x-window_index))
                overlap = max(0, min_less_int+15-window_index) + max(0, window_index+15-min_more_int)
                #prior windows are assumed to be 15 long regardless of regex length; match redundancy possible but unlikely
                if overlap <= 10: #maximum allowed overlap: 10 (2/3)
                    index_li.append(window_index)
                shannon_window_li2 = shannon_window_li2[1:]
            
            for index, index_value in enumerate(index_li):
                cropped_msa_li = [x[index_value:index_value+15] for x in protein_fasta_li2[1]]               
                cluster = clusterExpression(protein_fasta_li2[1], cropped_msa_li)
                if cluster != None: #failure when region is highly diverse, resulting in complex, non-specific regex
                    #regexSelect() requires indices of residues relative to save_msa
                    frame_mod = re.findall("#frameshift_modification=\((.*)\)", protein_fasta_li2[0][0])
                    if len(frame_mod) > 0 and len([x for x in frame_mod[0].split(",") if x.strip().isdigit()]) == 3:
                        adjust_li = [int(num.strip()) for num in frame_mod[0].split(",")]
                        adjust_li = adjust_li[:-1] #remove value used in hmmalignModified()
                        save_fasta_li2[1] = [x[:adjust_li[1]]+x[adjust_li[1]-adjust_li[0]:] for x in save_fasta_li2[1]]
                        modify_li = modify_li[:adjust_li[1]]+modify_li[adjust_li[1]-adjust_li[0]:]
                    protein_match_regex_li = [re.findall(cluster[0], x.replace("-", "")) for x in cropped_msa_li] #only consider sequences that match regex
                    combine_li = [len(x[:index_value].replace("-", ""))*3 if len(protein_match_regex_li[i]) > 0 else "" for i, x in enumerate(protein_fasta_li2[1])]
                    #assumes that protein_fasta_li2 and save_fasta_li2 have the same header order
                    #correction for gaps in msa_save nucleotide alignment
                    for combine_index, combine_value in enumerate(combine_li):
                        if combine_value != "":
                            save_value = combine_value
                            #correction for gaps that were not replaced by removeColumnGap()
                            while len(save_fasta_li2[1][combine_index][:save_value].replace("-", "")) < combine_value and \
                            len(save_fasta_li2[1][combine_index][:save_value].replace("-", "")) < len(save_fasta_li2[1][combine_index].replace("-", "")):
                                save_value += 1
                            combine_li[combine_index] = save_value
                    #correction for gaps removed by removeColumnGap()
                    combine_li = [x - modify_li[x] for x in combine_li if x != ""]              
                        
                    aa_position = round(index_value/(len(shannon_li)-15), 6)
                    if index == 0:
                        aa_position = 0
                    elif index == 1:
                        aa_position = 1                        
                    cropped_regex_json[annotation].append({
                        "regex": cluster[0],
                        "aa_position": aa_position,
                        "nucleotide_position": max(combine_li)
                    })    
            with open(protein_directory+'/protein.json', 'w') as outfile:
                for cds in cropped_regex_json:
                    cropped_regex_json[cds].sort(key=lambda x: x["aa_position"])
                json.dump(cropped_regex_json, outfile)
                    
################################################################################
    
def checkMSA(organism, dir_path=dir_path):
    #checks msa_save for consistent logic and constructs msa_build file (reference for hmm profile construction)
    #requires manual modification for translational frameshifts, limited to -1 frameshifts (not user friendly)
    input_nucleotide_msa_directory = dir_path + "/required/" + organism + "/hmm_profiles/msa_save"
    output_nucleotide_msa_directory = dir_path + "/required/" + organism + "/hmm_profiles/msa_build"
    output_protein_msa_directory = dir_path + "/required/" + organism + "/protein/msa"
    target_files = [f for f in os.listdir(input_nucleotide_msa_directory) if os.path.isfile(os.path.join(input_nucleotide_msa_directory, f))]
    for file in target_files:
        if ".fa" in file:
            read_fasta_li2 = readFASTA(input_nucleotide_msa_directory + "/" + file) #nucleotide
            read_fasta_li2[1] = [x.upper() for x in read_fasta_li2[1]]
            """
            consensus_str = ""
            index = 0
            for index in range(len(read_fasta_li2[1][0])):
                column = [read_fasta_li2[1][x][index] for x in range(len(read_fasta_li2[1]))]
                consensus_str += max(column, key=column.count)
            """
            #location of modified residue/s for frameshift correction
            frameshift_li = [] #slide distance/direction, msa_save index, msa_build index (gaps removed)
            consensus_str = "".join([max([y[x] for y in read_fasta_li2[1]], key=[y[x] for y in read_fasta_li2[1]].count) for x in range(len(read_fasta_li2[1][0]))])
            #find translation frameshift location
            if organism in ["PRRSV1", "PRRSV2", "PEDV", "PDCoV"]:
                fs_regex, fs_index = "TTTAAACTG.TAGCCGCCAGCGGCTTGACCCGCTGTGG", 10 #PRRSV
                if organism == "PEDV":
                    fs_regex, fs_index = "AGCACTGATATGGCTTATTTAAACGAGTACGGGGCTCTA", 9 #PEDV
                if organism == "PDCoV":
                    fs_regex, fs_index = "AATTCGGCTTATTTAAACG.GTAACGGGTTCTAGTGA", 16 #PDCoV
                fs_regex = "-*".join([x for x in fs_regex])
                fs_match_regex = re.findall(fs_regex, consensus_str)
                if len(fs_match_regex) > 0 and min([translateNucleotides(read_fasta_li2[1][x], 0)[:-1].count("x") for x in range(len(read_fasta_li2[1]))]) > 0:
                    add_adjust = len(re.findall("^"+"-*".join(["[A-Z]" for x in range(fs_index)]), fs_match_regex[0])[0])
                    fs_regex_adjusted_index = consensus_str.index(fs_match_regex[0]) + add_adjust
                    frameshift_li = [1, fs_regex_adjusted_index]
                    #track frequency of nucleotides adjacent to translational frameshift location
                    adjacent_nucleotide_di = {x:0 for x in "ATGC"} #"AUGC"}
                    for x in range(len(read_fasta_li2[1])):
                        for nucleotide in read_fasta_li2[1][x][fs_regex_adjusted_index-1:fs_regex_adjusted_index+1]:
                            if nucleotide in adjacent_nucleotide_di:
                                adjacent_nucleotide_di[nucleotide] += 1
                    #compensating for translational frameshift by inserting least likely nucleotide match as lower()
                    for x in range(len(read_fasta_li2[1])):
                        copy_adjacent_nucleotide_di = {x:adjacent_nucleotide_di[x] for x in adjacent_nucleotide_di}
                        for nucleotide in read_fasta_li2[1][x][fs_regex_adjusted_index-1:fs_regex_adjusted_index+1]:
                            if nucleotide in copy_adjacent_nucleotide_di:
                                del copy_adjacent_nucleotide_di[nucleotide]
                        read_fasta_li2[1][x] = read_fasta_li2[1][x][:fs_regex_adjusted_index] + min(copy_adjacent_nucleotide_di, key=copy_adjacent_nucleotide_di.get).lower() + read_fasta_li2[1][x][fs_regex_adjusted_index:]
            #remove sequences from MSA that cannot be translated (premature stop codon/s)
            protein_li = [translateNucleotides(x, 0) for x in read_fasta_li2[1]]
            start_li, end_li = [x[0] for x in protein_li], [x[-1] for x in protein_li]
            start_aa, end_aa = max(start_li, key=start_li.count), max(end_li, key=end_li.count)
            for index, protein_str in enumerate(protein_li):
                if ("x" in protein_str[:-1] or                      #premature stop codons
                    (start_aa != "M" and protein_str[0] == "M") or  #consistent start
                    (end_aa != "x" and protein_str[-1] == "x")):    #consistent end
                    protein_li[index] = ""
                    for l in range(2):
                        read_fasta_li2[l][index] = ""
            if protein_li.count("") > 0:
                print(protein_li.count(""), "sequences eliminated from", file, file=sys.stderr)
                protein_li = [x for x in protein_li if x != ""]
                for l in range(2):
                    read_fasta_li2[l] = [x for x in read_fasta_li2[l] if x != ""]
                    
            #remove columns by gap proportion (avoid problems with hmmalign)
            remove_gap_out_li = removeColumnGap(read_fasta_li2[1], 0.5) #translate indices between save_msa and build_msa
            read_fasta_li2[1] = remove_gap_out_li[0]
            if len(frameshift_li) > 0:
                frameshift_li.append(frameshift_li[1] - remove_gap_out_li[1][frameshift_li[1]])
                #add location of modified residue/s for frameshift correction
                #hmmalignModified() requires indices of modified residues relative to build_msa
                #regexSelect() requires indices of residues relative to save_msa
                for i, head in enumerate(read_fasta_li2[0]):
                    read_fasta_li2[0][i] += "#frameshift_modification=("+",".join([str(x) for x in frameshift_li])+")"
            #write protein fasta and align
            writeFASTA(output_protein_msa_directory+"/"+file, [read_fasta_li2[0], protein_li])
            mafftAlign(output_protein_msa_directory+"/"+file, output_protein_msa_directory+"/"+file)
            #write msa_build; do not re-align after removing sequences/columns (avoid disrupting frameshift index)
            writeFASTA(output_nucleotide_msa_directory + "/" + file, read_fasta_li2)
                
################################################################################

def buildHmm(organism, dir_path=dir_path):    
    #use hmmbuild to create hmm profiles (requires bash)
    target_dir = dir_path + "/required/" + organism + "/hmm_profiles/msa_build"
    build_dir = dir_path + "/required/" + organism + "/hmm_profiles/hmm_build"
    subprocess.call("mkdir -p "+build_dir, shell=True)
    target_file_li = [f for f in os.listdir(target_dir) if os.path.isfile(os.path.join(target_dir, f)) and ".fa" in f]
    for target_file in target_file_li:
        os_call = "hmmbuild "+build_dir+"/"+re.findall("^[^.]*", target_file)[0]+".hmm "+target_dir+"/"+target_file+" > /dev/null"
        subprocess.call(os_call, shell=True)
        
################################################################################

def buildGFF(organism, dir_path=dir_path):
    #builds reference template gff3 file (if one does not already exist)
    #user may make modifications to the template file or provide their own
    organism_path = dir_path+"/required/"+organism
    if not os.path.isdir(organism_path+"/gff3"):
        subprocess.call("mkdir -p "+organism_path+"/gff3", shell=True)
    if not os.path.isfile(organism_path+"/gff3/template.gff3"):
#       print("Required sample template could not be found for "+organism+". Generating..." )
        gff3_di = {}
        for protein_file in os.listdir(organism_path+"/protein/msa"):
            if ".fa" in protein_file:
                protein = re.sub("\.fa.*$", "", protein_file)
                protein_li = readFASTA(organism_path+"/protein/msa/"+protein_file)[1]
                #alignment unnecessary; variable length may impede appropriate treatment
                protein_li = [x.replace("-", "") for x in protein_li]
                #all gene and CDS are required to have start/stop codons
                if len([x for x in protein_li if x[0] == "M" and x[-1] == "x"]) > len(protein_li)/2:
                    if "gene/CDS" not in gff3_di:
                        gff3_di["gene/CDS"] = []
                    gff3_di["gene/CDS"].append(protein)
                else:
                    if "mat_peptide" not in gff3_di:
                        gff3_di["mat_peptide"] = []
                    gff3_di["mat_peptide"].append(protein)      
        #find parent of all mat_peptide
        parent_di = {}
        if "gene/CDS" in gff3_di and "mat_peptide" in gff3_di:
            #compare unaligned sequences in required reference files
            for mat_peptide in gff3_di["mat_peptide"]:
                sub_di = {re.sub("\.fa.*$", "", f):f for f in os.listdir(organism_path+"/hmm_profiles/msa_save")}
                if mat_peptide in sub_di:
                    parent_di[mat_peptide] = []
                    mat_peptide_li = list(set([x.replace("-", "") for x in readFASTA(organism_path+"/hmm_profiles/msa_save/"+sub_di[mat_peptide])[1]]))
                    for gene_cds in gff3_di["gene/CDS"]:
                        if gene_cds in sub_di:
                            match_count, cut_count_li = 0, []
                            gene_cds_li = [x.replace("-", "") for x in readFASTA(organism_path+"/hmm_profiles/msa_save/"+sub_di[gene_cds])[1]]
                            for mat_peptide_str in mat_peptide_li:
                                for gene_cds_str in gene_cds_li:
                                    find_mat_peptide = gene_cds_str.find(mat_peptide_str)
                                    if find_mat_peptide != -1:
                                        match_count += 1
                                        #number of cuts necessary to produce mat_peptide
                                        cut_count_li.append(0)
                                        if find_mat_peptide > 0:
                                            cut_count_li[-1] += 1
                                        if len(gene_cds_str) - find_mat_peptide - len(mat_peptide_str) > 3: #ignores stop codon
                                            cut_count_li[-1] += 1
                                        break
                            if match_count > len(mat_peptide_li)/2:
                                parent_di[mat_peptide].append([gene_cds, len(gene_cds_li[0]), round(sum(cut_count_li)/len(cut_count_li)), match_count/len(mat_peptide_li)])
        gff3_li2 = []
        if "mat_peptide" in gff3_di:
            for mat_peptide in gff3_di["mat_peptide"]:
                mat_peptide_value_li = ["ID="+mat_peptide, "Name="+mat_peptide, "Parent="]
                if mat_peptide in parent_di and len(parent_di[mat_peptide]) > 0:
                    #remove parents that have more than the minimum cut count (fewer cuts, more likely parent)
                    min_cut_count = min([x[2] for x in parent_di[mat_peptide]])
                    parent_di[mat_peptide] = [x for x in parent_di[mat_peptide] if x[2] == min_cut_count]
                    #add longest parent only
#                   mat_peptide_value_li[-1] += max(parent_di[mat_peptide], key = lambda x: x[1])[0]+"-gene"
                    #add all parents (sorted by length) separated by comma
                    mat_peptide_value_li[-1] += ",".join([x[0]+"-gene" for x in sorted(parent_di[mat_peptide], key = lambda x: x[1])])
                elif len(parent_di) == 0:
                    mat_peptide_value_li[-1] += "polyprotein-gene"
                else:
                    print("Missing parent for", mat_peptide, file=sys.stderr)
                    mat_peptide_value_li[-1] += "MISSING"
                gff3_li2.append(["sample", ".", "mat_peptide", "start", "end", ".", "+", ".", ";".join(mat_peptide_value_li)])
        if "gene/CDS" not in gff3_di:
            gff3_di["gene/CDS"] = ["polyprotein"]
        for x in range(2):
            for gene_cds in gff3_di["gene/CDS"]:
                if x == 0: 
                    gff3_li2.append(["sample", ".", "gene", "start", "end", ".", "+", ".", ";".join(["ID="+gene_cds+"-gene", "ID="+gene_cds, "gene="+gene_cds])])
                else: 
                    gff3_li2.append(["sample", ".", "CDS", "start", "end", ".", "+", ".", ";".join(["ID="+gene_cds+"-cds", "ID="+gene_cds+" CDS", "Parent="+gene_cds+"-gene"])])
        #write file
        write_file = open(organism_path+"/gff3/template.gff3", "w")
        write_file.write("##gff-version 3\n")
        for line in gff3_li2:
            write_file.write("\t".join(line) + "\n")
        write_file.close()
            
################################################################################  
            
def buildBlastdb(dir_path=dir_path):
    #build BLAST database from msa_save sequences
#   min_median_len = None
    read_fasta_li3 = []
    blastdb_dir = dir_path+"/required/blastdb"
    organism_dir_li = [d for d in os.listdir(dir_path+"/required") \
                       if os.path.isdir(dir_path+"/required/"+d) and d != "blastdb"]
    for organism in organism_dir_li:
        fasta_li = [dir_path+"/required/"+organism+"/hmm_profiles/msa_save/"+f for f \
                    in os.listdir(dir_path+"/required/"+organism+"/hmm_profiles/msa_save") \
                    if os.path.isfile(dir_path+"/required/"+organism+"/hmm_profiles/msa_save/"+f) and ".fa" in f]
        len_li = []
        for fasta in fasta_li:
            read_fasta_li2 = readFASTA(fasta)
            read_fasta_li2[0] = [x+"--"+organism for x in read_fasta_li2[0]]
            read_fasta_li2[1] = [x.replace("-", "") for x in read_fasta_li2[1]]
            len_li.append(len(read_fasta_li2[0]))
            read_fasta_li3.append(read_fasta_li2)
    """
        fasta_median_len = int(stats.median(len_li))
        if min_median_len == None or fasta_median_len < min_median_len:
            min_median_len = fasta_median_len
    #avoid sampling bias with equal random sampling
    for i, li2 in enumerate(read_fasta_li3):
        if len(li2[0]) > min_median_len:
            index_li = list(range(0, len(li2[0])))
            sample_index_li = random.sample(index_li, min_median_len)
            for n in range(2):
                read_fasta_li3[i][n] = [li2[n][x] for x in sample_index_li]
    """
    read_fasta_li2 = [[head for read_fasta_li2 in read_fasta_li3 for head in read_fasta_li2[0]], \
                      [seq for read_fasta_li2 in read_fasta_li3 for seq in read_fasta_li2[1]]]
    writeFASTA(blastdb_dir+"/blastdb.fasta", read_fasta_li2)
    os_call = "makeblastdb -dbtype nucl -in "+blastdb_dir+"/blastdb.fasta"
    subprocess.call(os_call, shell=True)

################################################################################   
################################################################################            

class Annotate:
    #classify all functions that call global variable/s
    def __init__(self, organism, temp_dir=""):
        #value appended to file names, avoids overwriting files when running in parallel.
        self.dir_path = dir_path
        self.required_path = self.dir_path + "/required"
        self.temporary_path = self.dir_path + "/temporary"
        if os.path.isdir(temp_dir):
            self.temporary_path = temp_dir
        if organism in valid_organism_li:
            self.organism = organism
        else:
            raise Exception("invalid input detected:", organism)
        self.total_cds = [re.findall("^[^.]*", f)[0] for f in os.listdir(self.required_path+"/"+self.organism+"/hmm_profiles/msa_save") if ".fa" in f]
        with open(self.required_path+"/"+self.organism+"/protein/protein.json") as json_file:  
            cropped_regex_json = json.load(json_file)
        for cds in cropped_regex_json:
            cropped_regex_json[cds].sort(key=lambda x: x["aa_position"])
        self.match_residues_di = cropped_regex_json
        self.msa_len_li, self.translational_frameshift_di = [], {}
        
        for cds in self.total_cds:
            self.msa_len_li.append(len(readFASTA(self.required_path+"/"+self.organism+"/protein/msa/"+cds+".fasta")[1][0]))
            
            translational_frameshift_li = []
            msa_build = self.required_path+"/"+self.organism+"/hmm_profiles/msa_build/"+cds+".fasta"
            open_file = open(msa_build, "r")
            head_line = open_file.readline()
            open_file.close()
            #hmmalignModified() requires indices of modified residues relative to build_msa
            frame_mod = re.findall("#frameshift_modification=\((.*)\)", head_line)
            if len(frame_mod) > 0 and len([x for x in frame_mod[0].split(",") if x.strip().isdigit()]) == 3:
                translational_frameshift_li = [int(num.strip()) for num in frame_mod[0].split(",")]
                translational_frameshift_li = translational_frameshift_li[:-2] + translational_frameshift_li[-1:] #remove value used in regexSelect()
                if translational_frameshift_li[0] < 0:
                    translational_frameshift_li = ["[A-Z][a-z]{"+str(translational_frameshift_li[0])+"}$"] + translational_frameshift_li
                elif translational_frameshift_li[0] > 0:
                    translational_frameshift_li = ["[A-Z]\-{"+str(translational_frameshift_li[0])+"}$"] + translational_frameshift_li
            self.translational_frameshift_di[cds] = translational_frameshift_li
            
    ################################################################################

    #https://stackoverflow.com/questions/3519565/find-the-indexes-of-all-regex-matches

    def matchRegex(self, dna_str, choose_match_li):
        #match regular expressions in protein.json to determine position of features
        #Note - regex can return more than one positive match
        total_match_di = {}
        #merge identical regex, which is possible from overlapping features
        regex_di3 = {}
        for cds in choose_match_li:
            if cds in self.match_residues_di:
                for cds_di in self.match_residues_di[cds]:
                    if cds_di["regex"] not in regex_di3:
                        regex_di3[cds_di["regex"]] = {}
                    regex_di3[cds_di["regex"]][cds] = dict(cds_di)
                    del regex_di3[cds_di["regex"]][cds]["regex"]
        #match regex against translated (three frames) input nucleotide sequence     
        protein_li = [translateNucleotides(dna_str, i) for i in range(3)]
        for regex in regex_di3:
            for reading_frame, protein in enumerate(protein_li):
                match_li2 = [(m.start(), m.group()) for m in re.finditer(regex, protein)]
                for cds in regex_di3[regex]:
                    if cds not in total_match_di:
                        total_match_di[cds] = []
                    for match in match_li2:
                        total_match_di[cds].append([reading_frame, match[0], regex_di3[regex][cds]["aa_position"], regex_di3[regex][cds]["nucleotide_position"], match[1]])
        
        #Attempt to force fit for start/stop if not already matched
        for regex in regex_di3:
            for cds in regex_di3[regex]:
                aa_position = regex_di3[regex][cds]["aa_position"]
                #if len(total_match_di[cds]) > 0 and aa_position == 0 or aa_position == 1:
                if aa_position == self.match_residues_di[cds][0]["aa_position"] or aa_position == self.match_residues_di[cds][-1]["aa_position"]:
                    if cds in total_match_di and aa_position not in [x[2] for x in total_match_di[cds]]: #minimum regex match count requirement, 1
                        for reading_frame, protein in enumerate(protein_li):
                            modify_regex_li = [re.findall("\[[A-Za-z]+\]|[A-Za-z]\??", x) for x in regex.split("|")]
                            for split_index in range(15): #maximum window size is 15; regular expression may be shorter and vary in length
                                mod_regex = "|".join(["".join(x[:split_index]+["."]+x[split_index+1:]) for x in modify_regex_li])
                                match_li2 = [(m.start(), m.group()) for m in re.finditer(mod_regex, protein)]
                                for cds in regex_di3[regex]:
                                    for match in match_li2:
                                        add = [reading_frame, match[0], regex_di3[regex][cds]["aa_position"], regex_di3[regex][cds]["nucleotide_position"], match[1]]
                                        if add not in total_match_di[cds]:
                                            total_match_di[cds].append(add)
        #reading frame, translated aa position (as index), position relative to aa conseq, position relative to hmm profile, match string
        return(total_match_di)

    ################################################################################
    
    def hmmbuild(self, segment_key, crop_msa_build_li):
        read_fasta_li2 = readFASTA(self.required_path+"/"+self.organism+"/hmm_profiles/msa_build/"+segment_key+".fasta")
        read_fasta_li2[1] = [x[crop_msa_build_li[0]:crop_msa_build_li[1]] for x in read_fasta_li2[1]]
        msa_build = self.temporary_path+"/crop.fasta"
        hmm_build = self.temporary_path+"/crop.hmm"
        writeFASTA(msa_build, read_fasta_li2)
        os_call = "hmmbuild "+hmm_build+" "+msa_build+" > /dev/null"
        subprocess.call(os_call, shell=True)
        os.remove(msa_build)
        return(hmm_build) #remove hmm_build file after hmmalign
    
    ################################################################################
    
    def hmmalign(self, dna_segment_str, hmm_build, force_match_end = False):
        #perform hmmalign on input sequence (requires bash)
        hmmalign_input = self.temporary_path+"/hmmalign_input.fasta"
        open_file = open(hmmalign_input, "w")
        open_file.write(">Sequence\n" + dna_segment_str)
        open_file.close()
        os_call = "hmmalign "+hmm_build+" "+hmmalign_input
        read_str = subprocess.check_output(os_call, shell=True).decode("utf-8")
        os.remove(hmmalign_input)
        trim_sequence_str = "".join([x[8:].strip() for x in read_str.split("\n") if "Sequence" in x and "#" not in x])
        if force_match_end: #terminal ends are known to be matching
############force align adjacent gap/lower
            
            
            #treatment for terminal ends
            #force fit condition: non-matching (by hmmalign) ends are re-defined as matching if enough nucleotides/gaps are available
            end_li = [re.findall("^[-a-z]*", trim_sequence_str)[0], re.findall("[-a-z]*$", trim_sequence_str)[0]]
            trim_sequence_str = trim_sequence_str[len(end_li[0]):len(trim_sequence_str)-len(end_li[1])]
            for i, end in enumerate(end_li):
                gap_count = end.count("-")
                lower_count = len(end.replace("-", ""))
                #add corrected ends to trim_sequence_str conditionally using multiplicative string logic
                if gap_count > 0:
#
                    #move misaligned residues to sides (easier to remove with end treatments)
                    trim_sequence_str = ("-"*(gap_count-lower_count) + end.replace("-", "")[:-gap_count] + end.replace("-", "")[-gap_count:].upper())*abs(i-1) + \
                    trim_sequence_str + (end.replace("-", "")[:gap_count].upper() + end.replace("-", "")[gap_count:] + ("-"*(gap_count-lower_count)))*i
                    #move misaligned residues to center (harder to remove with end treatments)
                    #trim_sequence_str = (end.replace("-", "")[-gap_count:].upper() + end.replace("-", "")[:-gap_count] + "-"*(gap_count-lower_count))*abs(i-1) + \
                    #trim_sequence_str + (("-"*(gap_count-lower_count)) + end.replace("-", "")[gap_count:] + end.replace("-", "")[:gap_count].upper())*i
#
        return(trim_sequence_str)

    ################################################################################
    
    def hmmalignModified(self, segment_key, dna_str, separate_match_li2, start, stop, reading_frame=0):
        #prepares and processes sequence for hmmalign, interprets result
    
        #check for reading frame changes in matched region/s
        match_li2, match_li3 = [separate_match_li2[0]], []
        for match in separate_match_li2[1:]:
            if match[0] == match_li2[0][0]:
                match_li2.append(match)
            else:
                match_li3.append(match_li2)
                match_li2 = [match]
        match_li3.append(match_li2)
        translate_str = "".join([translateNucleotides(dna_str[x[0][1]*3+x[0][0]:(x[-1][1]+len(x[-1][4]))*3+x[-1][0]], 0) for x in match_li3])
                        
        translational_frameshift_li = self.translational_frameshift_di[segment_key][:]
        hmm_build = self.required_path+"/"+self.organism+"/hmm_profiles/hmm_build/"+segment_key+".hmm"
        
        #tracks the corrected homologous sequence length relative to cropped alignments in else statement
        #necessary for finding location of translational frameshift/s; unable to reliably determine reading frame
        crop_li3 = [] #crop_dna_li, crop_msa_build_li
        if translate_str[:-1].count("x") != 0 or separate_match_li2[0][2] != 0 or separate_match_li2[-1][2] != 1:
            trim_sequence_str = self.hmmalign(dna_str[start:stop], hmm_build)
        else:
            in_sequence_li = []
            for match_index, match_li2 in enumerate(match_li3[:-1]):
                #crop input and crop/rebuild temporary hmm profile
                crop_li2 = [match_li2[-1], match_li3[match_index+1][0]]
                #crop_dna_li and crop_msa_build_li should be approximately equal lengths (treatments for when either is longer)
                crop_dna_li = [crop_li2[0][1]*3+crop_li2[0][0], (crop_li2[1][1]+len(crop_li2[1][4]))*3+crop_li2[1][0]]
                dna_segment_str = dna_str[crop_dna_li[0]:crop_dna_li[1]]
                crop_msa_build_li = [crop_li2[0][3], crop_li2[1][3]+len(crop_li2[1][-1])*3]
                crop_li3.append([crop_dna_li, crop_msa_build_li])
                
                hmm_build = self.hmmbuild(segment_key, crop_msa_build_li)
                in_sequence_str = self.hmmalign(dna_segment_str, hmm_build, True) #profile is expected to be no longer than input sequence 
                os.remove(hmm_build)
                
                #ignore/remove terminal gaps and update crop_li3 (crop_msa_build_li section)
                for i, pattern in enumerate(["^\-*", "\-*$"]):
                    	terminal_match, scale_end_li = re.findall(pattern, in_sequence_str)[0], [1, -1]
                    	if len(terminal_match) > 0: #hmm_build too long
                    		crop_li3[-1][1][i] += scale_end_li[i] * len(terminal_match)
                    		in_sequence_str = re.sub(pattern, "", in_sequence_str)
                #make regex matched residues uppercase
                for i, pattern in enumerate([["^(?:[^a-zA-Z]*[a-zA-Z]){"+str(len(crop_li2[0][-1])*3)+"}", "^(?:[\-]*[a-z])*"],
                                             ["(?:[^a-zA-Z]*[a-zA-Z]){"+str(len(crop_li2[1][-1])*3)+"}$", "(?:[\-]*[a-z])*$"]]):
                    #find characters matched by match_li3 regex
                    terminal_match_li = re.findall(pattern[0], in_sequence_str)
                    if len(terminal_match_li) > 0:
                        	terminal_match, scale_end_li = re.findall(pattern[1], terminal_match_li[0])[0], [1, -1]
                        	lower_count = len([x for x in terminal_match if x.islower()])
                        	if lower_count > 0: #dna_segment too long
                        		in_sequence_str = re.sub("^"*abs(i-1)+terminal_match+"$"*i, terminal_match.upper(), in_sequence_str)
                #compares total potential frame changes in aligned area to known reading frame of matched regex
                sub_gap = (in_sequence_str.count("-") - len([x for x in in_sequence_str if x.islower()]))%3 #remainder returns positive number
                frame_change = crop_li2[0][0] - crop_li2[1][0]
                if sub_gap != frame_change and sub_gap != 3-frame_change and sub_gap != 3+frame_change:
                    in_sequence_str = in_sequence_str.replace("-", "").upper()
                    sub_gap = frame_change
                #adds alignment correction (end treatments removed aligned information)
                if in_sequence_str.replace("-", "").upper() == in_sequence_str: #no misaligned nucleotides for frameshift
                    sub_gap = crop_li2[0][0] - crop_li2[1][0]
                    crop_location = len(crop_li2[0][4])*3
                    if sub_gap == -1 or sub_gap == 2: #add gap
                        in_sequence_str = in_sequence_str[:crop_location] + in_sequence_str[crop_location].lower() + in_sequence_str[crop_location+1:]
                    elif sub_gap == 1 or sub_gap == -2: #sub lower
                        in_sequence_str = in_sequence_str[:crop_location] + "-" + in_sequence_str[crop_location:]
                in_sequence_li.append(in_sequence_str)
            for index, in_sequence_str in enumerate(in_sequence_li[:-1]):
                #overlap correction
                overlap_int = crop_li3[index][0][1] - crop_li3[index+1][0][0]
                if overlap_int > 0:
                    if in_sequence_str[-overlap_int:].isupper(): #end uppercase
                        in_sequence_li[index] = in_sequence_str[:-overlap_int]
                        crop_li3[index][0][1] -= overlap_int
                    elif in_sequence_li[index+1][:overlap_int].isupper(): #beginning uppercase
                        in_sequence_li[index+1] = in_sequence_li[index+1][overlap_int:]
                        crop_li3[index+1][0][0] += overlap_int
                    else:
                        print("overlap correction error: evidence of frameshift found in both sequences", file=sys.stderr) #?
            #merge aligned and unaligned sections of dna_str
            flat_crop_dna_li = [start] + [y for x in crop_li3 for y in x[0]] + [stop]
            out_sequence_li = [dna_str[flat_crop_dna_li[x]:flat_crop_dna_li[x+1]] for x in range(0, len(flat_crop_dna_li), 2)]
            trim_sequence_str = out_sequence_li[0]
            for index in range(len(in_sequence_li)):
                trim_sequence_str += in_sequence_li[index]
                trim_sequence_str += out_sequence_li[index+1]  
            
        end_li = [re.findall("^[-a-z]*", trim_sequence_str)[0], re.findall("[-a-z]*$", trim_sequence_str)[0]]
        trim_sequence_str = trim_sequence_str[len(end_li[0]):len(trim_sequence_str)-len(end_li[1])]
        #fragmentation of the modified input nucleotide sequence
        #lowercase characters and gaps ("-") in the hmm alignemnt are indicative of low homology regions where reading frame shifts are most likely
        modified_nucleotide_li = re.findall("[A-Z]+[a-z\-]*", trim_sequence_str)
        #adjust start value
        start += len("".join(re.findall("[a-zA-Z]", end_li[0])))
        #modify reading frame based on end_li[0] (front)
        front_sub_gap_len = len(end_li[0].replace("-", ""))
        new_reading_frame = 0
        while (front_sub_gap_len + new_reading_frame)%3 != reading_frame:
            new_reading_frame += 1  
        
        frameshift_li2 = [] #start, stop, type (translational=0, error=1)
        read_len_correction = 0
        for modified_nucleotide_index, modified_nucleotide_str in enumerate(modified_nucleotide_li[:-1]):
            #homolog frame must be preserved; lower cannot be part of substitution with upper or gap
            homolog_len = len(re.sub("[a-z]", "", "".join(modified_nucleotide_li[:modified_nucleotide_index+1])))
            read_len = len("".join(modified_nucleotide_li[:modified_nucleotide_index+1]).replace("-", ""))
            corrected_read_len = read_len + read_len_correction 
            #correct_homolog_len references crop_msa_build_li in crop_li3; parts of modified_nucleotide_li may not have been hmmalign-ed
            pseudo_align_len = corrected_read_len
            for crop_index, crop_li2 in enumerate(crop_li3[::-1]): #[crop_dna_li, crop_msa_build_li]
                if crop_li2[0][0]-start < pseudo_align_len:
                    pseudo_align_len -= crop_li2[0][0]-start - crop_li2[1][0]
                    break
            #unique condition for translational frameshits
            if len(translational_frameshift_li) == 3 \
            and len(re.findall(translational_frameshift_li[0], modified_nucleotide_str)) > 0 \
            and abs(translational_frameshift_li[2] - pseudo_align_len) <= 1:  #aligned length, target location, indel length correction
                location_difference = translational_frameshift_li[2] - pseudo_align_len
                slip_len = translational_frameshift_li[1]
                if slip_len > 0: #check slip direction (alternative direction not currently needed)
                    modified_nucleotide_li[modified_nucleotide_index] = modified_nucleotide_str[:-slip_len] #crop out gap/s
                    #corrective measure if hmmalignment slightly misplaces site of translational frameshift
                    if location_difference > 0: #move frameshift upstream
                        modified_nucleotide_li[modified_nucleotide_index] += modified_nucleotide_li[modified_nucleotide_index+1][:location_difference]
                        modified_nucleotide_li[modified_nucleotide_index+1] = modified_nucleotide_li[modified_nucleotide_index+1][location_difference:]
                    elif location_difference < 0: #move frameshift downstream
                        modified_nucleotide_li[modified_nucleotide_index+1] = modified_nucleotide_li[modified_nucleotide_index][location_difference:] + modified_nucleotide_li[modified_nucleotide_index+1]
                        modified_nucleotide_li[modified_nucleotide_index] = modified_nucleotide_li[modified_nucleotide_index][:location_difference]
                    corrected_read_len += location_difference
                    modified_nucleotide_li[modified_nucleotide_index] += modified_nucleotide_li[modified_nucleotide_index][-slip_len:]
                    frameshift_li2.append([corrected_read_len, corrected_read_len-slip_len, 0])
                    read_len_correction -= slip_len
            else: #must read next segment for presence of stop codons
                #reading frame determination
                sub_gap = (homolog_len-read_len)%3 #remainder returns positive number
                if sub_gap != 0:
                    end_gap_str = re.findall("[a-z\-]*$", modified_nucleotide_str)[0]
                    if len(end_gap_str) > 0:
                        modified_nucleotide_str = modified_nucleotide_str[:-len(end_gap_str)]
                    #add gap/s to preserve homolog_len; track changes in read_len with slip_len variable
                    if sub_gap == 2: # or sub_gap == -1: #skip nucleotide
                        modified_nucleotide_str = modified_nucleotide_str[:-1]+"-" + end_gap_str
                        frameshift_li = [corrected_read_len-1, corrected_read_len, 1]
                        slip_len = 1
                    elif sub_gap == 1: # or sub_gap == -2: #re-read nucleotide
                        modified_nucleotide_str += modified_nucleotide_str[-1].lower() + end_gap_str
                        frameshift_li = [corrected_read_len, corrected_read_len-1, 1]
                        slip_len = -1
                    read_str = "".join(modified_nucleotide_li[:modified_nucleotide_index+2])
                    sub_str = "".join(modified_nucleotide_li[:modified_nucleotide_index])+modified_nucleotide_str+modified_nucleotide_li[modified_nucleotide_index+1]
                    translate_read_str = translateNucleotides(read_str.replace("-", ""), new_reading_frame)
                    translate_sub_str = translateNucleotides(sub_str.replace("-", ""), new_reading_frame)
                    if translate_sub_str[:-1].count("x") < translate_read_str[:-1].count("x"):
                        modified_nucleotide_li[modified_nucleotide_index] = modified_nucleotide_str
                        frameshift_li2.append(frameshift_li)
                        read_len_correction += slip_len
        modified_nucleotide_str = "".join(modified_nucleotide_li).replace("-", "").upper()
        #correct frameshift locations based on adjusted start
        frameshift_li2 = [[x[0]+start, x[1]+start, x[2]] for x in frameshift_li2]
        nucleotide_out_li = end_li
        return(modified_nucleotide_str, frameshift_li2, nucleotide_out_li, new_reading_frame)

    ################################################################################
    
    def annotateSequence(self, dna_str, choose_match_li):
        #processes input sequence based on provided list of possible feature/s
        fix_dna_str = dna_str.upper().replace("U", "T").replace("-", "")
        total_annotate_li2 = []
        total_match_di = self.matchRegex(fix_dna_str, choose_match_li)
        for key in total_match_di:
            total_match_di[key].sort(key=lambda x: x[1])
            if len(total_match_di[key]) > 1: #check if enough conserved matches, minimum 2
                #separate values that fall out of order (possible duplicate sequences or false matches)
                separate_match_li3 = [[total_match_di[key][0]]]
                for value in total_match_di[key]:
                    if value not in separate_match_li3[-1]:
                        if value[2] > separate_match_li3[-1][-1][2]:
                            separate_match_li3[-1].append(value)            
                        else:
                            separate_match_li3.append([value])
                
                location_li3 = []
                msa_len = self.msa_len_li[self.total_cds.index(key)]
                for separate_match_index, separate_match_li2 in enumerate(separate_match_li3):
                    if len(separate_match_li2) > 1: #recheck count requirement post separation, minimum 2
                        match_start = separate_match_li2[0][0]+(separate_match_li2[0][1]*3) #reading frame + amino acid index
                        match_stop = separate_match_li2[-1][0]+(separate_match_li2[-1][1]+len(separate_match_li2[-1][-1]))*3 #reading frame + amino acid index + adjustment to end of regex pattern
                        #estimated start assuming feature is full length; helps match discontinuity (frame change/s)
                        estimate_start = match_start - msa_len*3*(0.10 + separate_match_li2[0][2])
                        estimate_stop = match_stop + msa_len*3*(0.10 + 1-separate_match_li2[-1][2])
                        location_li3.append([[match_start, match_stop], [estimate_start, estimate_stop]])
                    else:
                        separate_match_li3[separate_match_index] = []
                location_li3 = [[[int(z) for z in y] for y in x] for x in location_li3] #estimated values may be decimals
                separate_match_li3 = [x for x in separate_match_li3 if x != []]

                #merge overlapping locations
                len_location_li3 = None
                while len_location_li3 == None or len_location_li3 != len(location_li3):
                    for loc1_index, loc1 in enumerate(location_li3[:-1]):
                        loc2 = location_li3[loc1_index+1]
                        if loc1 != [] and loc2 != []:
                            if loc1[0][0] >= loc2[1][0] and loc1[0][1] <= loc2[1][1] or loc2[0][0] >= loc1[1][0] and loc2[0][1] <= loc1[1][1]:
                                location_li3[loc1_index] = [[min(loc1[0][0], loc2[0][0]), max(loc1[0][1], loc2[0][1])], \
                                                            [min(loc1[1][0], loc2[1][0]), max(loc1[1][1], loc2[1][1])]]
                                separate_match_li3[loc1_index] = separate_match_li3[loc1_index]+separate_match_li3[loc1_index+1]
                                location_li3[loc1_index+1], separate_match_li3[loc1_index+1] = [], [] 
                                break
                    len_location_li3 = len(location_li3)
                    location_li3 = [x for x in location_li3 if x != []]
                    separate_match_li3 = [x for x in separate_match_li3 if x != []]
                for loc_index, separate_match_li2 in enumerate(separate_match_li3):
                    total_annotate_li = ["" for x in range(7)]
                    total_annotate_li[0] = key+" full"
                    reading_frame = 0
                    start, stop = location_li3[loc_index][0][0], location_li3[loc_index][0][1]
                    if separate_match_li2[0][2] != 0:
                        #estimated start must be in same reading frame as separate_match_li2[0] value
                        start = location_li3[loc_index][1][0]
                        while abs(start%3) != separate_match_li2[0][0]:
                            start -= 1
                        if start < 0:
                            start = 0
                            reading_frame = separate_match_li2[0][0]
                    if separate_match_li2[-1][2] != 1:
                        stop = location_li3[loc_index][1][1]
                        if stop > len(fix_dna_str):
                            stop = len(fix_dna_str)
                    translate_str = translateNucleotides(fix_dna_str[start:stop], reading_frame)
                    frameshift_li2 = []
                    if translate_str[:-1].count("x") != 0 or separate_match_li2[0][2] != 0 or separate_match_li2[-1][2] != 1:
                        modified_nucleotide_str, frameshift_li2, nucleotide_out_li, reading_frame = self.hmmalignModified(key, fix_dna_str, separate_match_li2, start, stop, reading_frame)
                        translate_str = translateNucleotides(modified_nucleotide_str, reading_frame)
                        align_len, crop_len = len(re.sub("[a-z]", "", modified_nucleotide_str)), len(re.sub("[a-z]", "", "".join(nucleotide_out_li)))
                        if crop_len > 0:
                            if crop_len/(align_len+crop_len) < 0.05:
                                total_annotate_li[0] = key+" near complete"
                            else:
                                total_annotate_li[0] = key+" partial"                           
                        for i in range(2):
                            total_annotate_li[i+3] = "; ".join([",".join([str(x[0]), str(x[1]+1)]) for x in frameshift_li2 if x[-1] == i])
                    
                    #scale output values if input sequence contains gaps
                    if dna_str.count("-") > 0:
                        location_li = [start] + [y for x in frameshift_li2 for y in x[:2]] + [stop]
                        for index, location in enumerate(location_li):
                            while len(dna_str[:location].replace("-", "")) != location_li[index]:
                                location += 1
                            location_li[index] = location
                        start, stop = location_li[0], location_li[-1]
                        frameshift_li2 = [location_li[1:-1][x:x+2]+[frameshift_li2[i][-1]] for i, x in enumerate([y for y in range(0, len(location_li[1:-1]), 2)])]
                    
                    total_annotate_li[1] = str(start+1)+","+str(stop) #add one to start to account for genbank format
                    total_annotate_li[2] = str(reading_frame)
                    for i in range(2):
                        total_annotate_li[i+3] = "; ".join([",".join([str(x[0]), str(x[1]+1)]) for x in frameshift_li2 if x[-1] == i])
                    if translate_str[-1] == "x":
                        translate_str = translate_str[:-1]
                    total_annotate_li[-2:] = str(len(translate_str)), translate_str
                    total_annotate_li2.append(total_annotate_li)

        #0:translation type, 1:start/stop range, 2:reading frame, 3:location of translational frameshift, 4:homology adjustment (error correction), 5:sequence length, 6:AA sequence
        return(total_annotate_li2)

################################################################################
################################################################################
        
def continueMessage(input_value="Additional setup required. Continue? (y/n) ", \
                    response_y="You have selected to continue. Please wait.", \
                    response_n="Process terminated."):
    while True:
        output = input(input_value)
        if output.lower() == "y" or output.lower() == "yes":
            if len(response_y) > 0: 
                print(response_y, file=sys.stderr)
            return(True)
        elif output.lower() == "n" or output.lower() == "no":
            if len(response_n) > 0: 
                print(response_n, file=sys.stderr)
            return()
        else:
            print("Invalid response.", file=sys.stderr)  

################################################################################

def organismDependencies(organism):
    #checks that all necessary dependencies exist; prompts user to install missing programs or build missing reference files
    global valid_organism_li
    if organism in valid_organism_li:
        msa_save_li, msa_build_li, hmm_build_li, protein_msa_li = [], [], [], []
        if os.path.isdir(dir_path+"/required/"+organism+"/hmm_profiles/msa_save"):
            msa_save_li =    [re.sub("\.fa.*", "", f) for f in os.listdir(dir_path+"/required/"+organism+"/hmm_profiles/msa_save") \
                              if os.path.isfile(dir_path+"/required/"+organism+"/hmm_profiles/msa_save/"+f) and f[0] != "."]
        if os.path.isdir(dir_path+"/required/"+organism+"/hmm_profiles/msa_build"):
            msa_build_li =   [re.sub("\.fa.*", "", f) for f in os.listdir(dir_path+"/required/"+organism+"/hmm_profiles/msa_build") \
                              if os.path.isfile(dir_path+"/required/"+organism+"/hmm_profiles/msa_build/"+f)]
        if os.path.isdir(dir_path+"/required/"+organism+"/hmm_profiles/hmm_build"):
            hmm_build_li =   [re.sub("\.hmm.*", "", f) for f in os.listdir(dir_path+"/required/"+organism+"/hmm_profiles/hmm_build") \
                              if os.path.isfile(dir_path+"/required/"+organism+"/hmm_profiles/hmm_build/"+f)]
        if os.path.isdir(dir_path+"/required/"+organism+"/protein/msa"):
            protein_msa_li = [re.sub("\.fa.*", "", f) for f in os.listdir(dir_path+"/required/"+organism+"/protein/msa") \
                              if os.path.isfile(dir_path+"/required/"+organism+"/protein/msa/"+f)]            
        try:
            with open(dir_path+"/required/"+organism+"/protein/protein.json") as json_file:  
                cropped_regex_json = json.load(json_file)
        except:
            cropped_regex_json = {}        
        if (len([x for x in msa_save_li if x not in msa_build_li]) != 0 or \
            len([x for x in msa_save_li if x not in hmm_build_li]) != 0 or \
            len([x for x in msa_save_li if x not in protein_msa_li]) != 0 or \
            len([x for x in msa_save_li if x not in cropped_regex_json]) != 0) or \
            not os.path.isfile(dir_path+"/required/"+organism+"/gff3/template.gff3"):
            if continueMessage(input_value="Additional setup required for "+organism+". Continue? (y/n) ", \
                               response_y="You have selected to continue. Please wait. Notice: this is a slow process.", \
                               response_n="You have chosen to skip setup for "+organism+"."):
                if len([x for x in msa_save_li if x not in msa_build_li]) != 0:
                    checkMSA(organism, dir_path)
                if len([x for x in msa_save_li if x not in hmm_build_li]) != 0:
                    buildHmm(organism, dir_path)
                """
                if len([x for x in msa_save_li if x not in protein_msa_li]) != 0:
                    #manually check protein MSAs, removing suspect recombinant/mutant sequences
                    #nucleotide hmm profile alignments more resilient to poor alignment and do not require similar treatment
                    exclude_di2 = { "PRRSV1": {}, #"ORF1ab":["A26843"], "nsp11":["A26843"]},
                                    "PRRSV2": {},
                                    "PEDV": {"ORF4":["KX580953"]},
                                    "SVA": {"2C":["KX173339"]},
                                    "CSFV": {"E0":["HI516623"],
                                             "E1":["HI516623"],
                                             "NS5B":["LC016722"]},
                                    "FMDV": {"2A":["MG372731", "AY304994"],
                                             "3B3":["MH426555"],
                                             "L":["DQ409183", "DQ409186", "DQ409189"]},
                                    "PDCoV": {} }
                    for key in exclude_di2[organism]:
                        protein_msa_path = dir_path+"/required/"+organism+"/protein/msa/"+key+".fasta"
                        read_fasta_li2 = readFASTA(protein_msa_path)
                        for accession in exclude_di2[organism][key]:
                            for head_index, head in enumerate(read_fasta_li2[0]):
                                if accession in head:
                                    read_fasta_li2[0] = read_fasta_li2[0][:head_index] + read_fasta_li2[0][head_index+1:]
                                    read_fasta_li2[1] = read_fasta_li2[1][:head_index] + read_fasta_li2[1][head_index+1:]
                        #may need to re-align instead of simply removing excess gaps
                        read_fasta_li2[1] = removeColumnGap(read_fasta_li2[1], 1)[0]
                        writeFASTA(protein_msa_path, read_fasta_li2)
                """
                if len([x for x in msa_save_li if x not in cropped_regex_json]) != 0:
                    regexSelect(dir_path+"/required/"+organism+"/protein", cropped_regex_json)
                if not os.path.isfile(dir_path+"/required/"+organism+"/gff3/template.gff3"):
                    buildGFF(organism, dir_path)
                print("Setup for "+organism+" finished successfully.", file=sys.stderr)
                return(True)
        else:
            return(True)
            
################################################################################

def annotateFASTA(input_organism_li, read_fasta_li2, choose_match_li=[]):
    #processes input fasta based on provided list of possible organism/s and feature/s
    global valid_organism_li
    annotate_li4 = [[] for x in read_fasta_li2[0]]
    input_organism_li = [x for x in input_organism_li if x in valid_organism_li]
    if len(input_organism_li) == 0:
        input_organism_li = valid_organism_li[:]
    #check required bash utilities
    required_program_li = ["hmmalign", "hmmbuild", "mafft", "awk"]
    if "prrsv_orf5_RFLP" in sys.argv[0] or len(input_organism_li) != 1:
        required_program_li += ["blastn", "makeblastdb", "blastdbcmd"]
    missing_program_li = [x for x in required_program_li if shutil.which(x) == None]
    if len(missing_program_li) > 0:
        print("The following dependencies could not be found:", ", ".join(missing_program_li)+".", file=sys.stderr)
        print("Please install and rerun.", file=sys.stderr)
    else:
        #check blast database
        blastdb = dir_path+"/required/blastdb/blastdb.fasta"
        os_call = "blastdbcmd -info -db "+blastdb+" > /dev/null"
        check_blastdb_li = [subprocess.call(os_call.replace(".fasta", ""), shell=True), subprocess.call(os_call, shell=True)]
        if "blastn" in required_program_li and 0 not in check_blastdb_li and continueMessage(input_value="BLAST database must be constructed. Continue? (y/n) ", \
                                                                                             response_n="Not an optional step. Process terminated."):
            buildBlastdb(dir_path)
            print("BLAST database successfully constructed.", file=sys.stderr)
            check_blastdb_li = [subprocess.call(os_call.replace(".fasta", ""), shell=True), subprocess.call(os_call, shell=True)]
        if check_blastdb_li[0] == 0:
            blastdb = blastdb.replace(".fasta", "")
        if "blastn" not in required_program_li or 0 in check_blastdb_li:
            #check selected organism dependencies
            for organism in input_organism_li:
                if not organismDependencies(organism):
                    print("Removing "+organism+" from searchable organisms.", file=sys.stderr)
                    input_organism_li = [x for x in input_organism_li if x != organism]
            #build random, non-existent directory name
            while True:
                temp_id = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for r in range(10))
                if not os.path.isdir(dir_path+"/temporary/temp"+temp_id):
                    break
            temp_dir = dir_path+"/temporary/temp"+temp_id
            subprocess.call("mkdir -p "+temp_dir, shell=True)
            try:
                #identify organism with BLAST
                fasta_di = {} #track original FASTA order
                if len(input_organism_li) > 1:
                    if "blastn" in required_program_li:
                        for index, sequence in enumerate(read_fasta_li2[1]):
                            writeFASTA(temp_dir+"/blast_in.fasta", [[str(index)], [sequence.replace("-", "")]])
############################blastn call/interpretation is largely reliant upon default behaviors; may require customized options
                            os_call = "blastn -db "+blastdb +" -query "+temp_dir+"/blast_in.fasta -outfmt 6 | awk -F'\t' '{print $2}' | awk -F'--' '{print $NF}' > "+temp_dir+"/blast_out.txt"
                            subprocess.call(os_call, shell=True)
                            with open(temp_dir+"/blast_out.txt") as blast_out: # Use file to refer to the file object
                                read_out = blast_out.read().split("\n")[:-1] #ignore last newline character
                                if len(read_out) > 0:
                                    organism_count_di = {x:read_out.count(x) for x in set(read_out)}
                                    max_key_count = max(organism_count_di.items(), key=lambda x : x[1])
                                    total_key_count = sum(organism_count_di.values())
                                    if max_key_count[0] in input_organism_li and max_key_count[1]/total_key_count > 0.5:
                                        if max_key_count[0] not in fasta_di:
                                            fasta_di[max_key_count[0]] = []
                                        fasta_di[max_key_count[0]].append([index, sequence])
                elif len(input_organism_li) == 1:
                    fasta_di[input_organism_li[0]] = [[index, sequence] for index, sequence in enumerate(read_fasta_li2[1])]
                #annotate by organism
                for organism in fasta_di:
                    annotate = Annotate(organism, temp_dir)
                    choose_match_li = [x for x in choose_match_li if x in annotate.total_cds]
                    if choose_match_li == []:
                        choose_match_li = annotate.total_cds
                    for di_li in fasta_di[organism]: #use original FASTA order
                        total_annotate_li2 = annotate.annotateSequence(di_li[1], choose_match_li)
                        if len(total_annotate_li2) > 0:
                            annotate_li4[di_li[0]] = [organism, total_annotate_li2]
            except:
                traceback.print_exc(file=sys.stderr)
                print(read_fasta_li2[0], file=sys.stderr) #?
            subprocess.call("rm -rf "+temp_dir, shell=True)
    return(annotate_li4)

################################################################################

#https://m.ensembl.org/info/website/upload/gff3.html
#http://gmod.org/wiki/GFF3#GFF3_Format
#https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
#https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
    
def outputGFF(read_fasta_li2, annotate_li4):
    #accepts annotated output (non-standard format) and uses reference template to produce gff3 file (standard format)
    organism_li = list(set([x[0] for x in annotate_li4 if len(x) == 2]))
    gff_di = {}
    for organism in organism_li:
        gff_path = dir_path+"/required/"+organism+"/gff3/template.gff3"
        if os.path.isfile(gff_path):
            read_gff = open(gff_path, "r")
            read_gff_li2 = [line.strip().split("\t") for line in read_gff]
            read_gff.close()
            gff_di[organism] = read_gff_li2
        else:
            print("GFF3 template file could not be found for "+organism, file=sys.stderr)
    
    gff_output_str = "##gff-version 3\n"
    for fasta_index, annotate_li3 in enumerate(annotate_li4):
        if len(annotate_li3) == 2 and annotate_li3[0] in gff_di:
            annotate_li2 = annotate_li3[1]
            feature_li = [re.sub(" full| near complete| partial", "", x[0]) for x in annotate_li2]
            gff_li3 = [[x] for x in gff_di[annotate_li3[0]]]
            for gff_index, gff_li2 in enumerate(gff_li3):
                gff_template_li = gff_li2[0][:]  #first list always template (maintain order)
                #0:translation type, 1:start/stop range, 2:reading frame, 3:location of translational frameshift, 4:homology adjustment (error correction), 5:sequence length, 6:AA sequence
                for feature_index, annotate_li in enumerate(annotate_li2):
                    feature = feature_li[feature_index]
                    note_str = re.sub(feature, "", annotate_li[0]).strip() #completeness
                    find_parent = re.findall("[Pp][Aa][Rr][Ee][Nn][Tt]\=([^\;]*)|$", gff_template_li[-1])[0]
                    if len(re.findall("[^a-zA-Z0-9]"+feature+"[^a-zA-Z0-9]", gff_template_li[-1].replace(find_parent, ""))) > 0: #regex should be specific enough for single character feature names
                        #apply feature annotations to gff file line
                        gff_template_li[0] = read_fasta_li2[0][fasta_index]
                        gff_template_li[1] = "United States Swine Pathogen Database"
                        #edit/modify note
                        if len(annotate_li[3]) > 0:
                            translational_frameshift_li = [int(x.strip()) for x in annotate_li[3].split(",")] #assumes maximum of one translational frameshift per genomic feature
                            note_str += ", "+"{:+}".format(translational_frameshift_li[0]-translational_frameshift_li[1]-1)+" translational frameshift at "+str(translational_frameshift_li[0])
                        if len(annotate_li[4]) > 0:
                            note_str += ", adjustments made to account for premature stop codon/s"
                        note_index = gff_template_li[-1].upper().find("NOTE=")
                        if note_index >= 0:
                            gff_template_li[-1] = gff_template_li[-1][:note_index]+note_str+", "+gff_template_li[-1][note_index:]
                        else:
                            gff_template_li[-1] += ";Note="+note_str
                            
                        #translational frameshift and homology adjustment treatment
                        frameshift_str = annotate_li[3] + "; " + annotate_li[4]
                        frameshift_li2 = [[y.strip() for y in x.split(",")] for x in frameshift_str.split(";") if x.strip() != ""]
                        frameshift_li2 = sorted(frameshift_li2, key=lambda x: int(x[0]))
                        
                        start, stop = tuple([x.strip() for x in annotate_li[1].split(",")])
                        frameshift_li = [start]+[x for y in frameshift_li2 for x in y]+[stop]
                        frameshift_li2 = [[frameshift_li[x], frameshift_li[x+1]] for x in range(0, len(frameshift_li), 2)]
                        index_gff_li2 = [gff_template_li[:3]+x+gff_template_li[5:] for x in frameshift_li2]
                        if annotate_li[2] != 0:
                            index_gff_li2[0][-2] = annotate_li[2]
                        gff_li3[gff_index] += index_gff_li2
                       
                        #missing parent treatment (inference)
                        if len(find_parent) > 0:
                            parent_li = find_parent.split(",")
                            for parent in parent_li:
                                #check that parent does not exist in annotate_li2
                                if len([x for x in feature_li if len(re.findall("(?:^|[^a-zA-Z0-9])"+x+"(?:[^a-zA-Z0-9]|$)", parent)) > 0]) == 0:
                                    #find appropriate line in gff3 for gene
                                    for index, li2 in enumerate(gff_li3):
                                        template_li = li2[0][:]
                                        if parent in re.sub("[Pp][Aa][Rr][Ee][Nn][Tt]\=([^\;]*)", template_li[-1]):
                                            if template_li[0] != read_fasta_li2[0][fasta_index]:
                                                template_li[0] = read_fasta_li2[0][fasta_index]
#                                               template_li[1] = "United States Swine Pathogen Database" #not supported
                                                #edit/modify note
                                                note_index = template_li[-1].upper().find("NOTE=")
                                                if note_index >= 0:
                                                    template_li[-1] = template_li[-1][:note_index]+"Infered, "+template_li[-1][note_index:]
                                                else:
                                                    template_li[-1] += ";Note=Infered"
                                            
                                            index_gff_li2 = [template_li[:3]+x+template_li[5:] for x in frameshift_li2]
                                            if annotate_li[2] != 0:
                                                index_gff_li2[0][-2] = annotate_li[2]
                                            gff_li3[index] += index_gff_li2
                                            #order and merge values
                                            gff_li3[index] = gff_li3[index][0] + sorted(gff_li3[index][1:], key=lambda x: int(x[3]))
                                            i = 0
                                            while i < len(gff_li3[index])-1:
                                                if int(gff_li3[index][i+1][3]) - int(gff_li3[index][i][4]) == 1:
                                                    gff_li3[index][i][4] = gff_li3[index][i+1][4]
                                                    gff_li3[index] = gff_li3[index][:i+1]+gff_li3[index][i+2:]
                                                else:
                                                    i += 1
            for gff_li2 in gff_li3:
                if len(gff_li2) > 1:
                    gff_output_str += "\n".join(["\t".join(x) for x in gff_li2[1:]])+"\n"
    return(gff_output_str)

################################################################################
################################################################################

def readable_list(_s):
    #https://stackoverflow.com/questions/53981845/grammatically-correct-human-readable-string-from-list-with-oxford-comma
    if len(_s) < 3:
        return ' and '.join(map(str, _s))
    *a, b = _s
    return f"{', '.join(map(str, a))}, and {b}"

#gather system arguments
valid_output_format_li = ["GFF3"]
valid_output_format_str = ", ".join(valid_output_format_li)
if "sequence_annotate" in sys.argv[0]:
    if "-h" in sys.argv[1:] or "-help" in sys.argv[1:]:
        print("""USAGE
  python3 sequence_annotate.py in > out

OPTIONAL ARGUMENTS
 -h or -help
   Print USAGE, DESCRIPTION and ARGUMENTS
 -organism=<Organisms>
   List (comma delimited) possible organism identity of sequences in input fasta
 -in=<Input_File>
   Input FASTA file name
 -out=<Output_File>
   Output file name
 -out_format=<Output_Format>
   Output file format ("""+", ".join(valid_output_format_li)+""")

DESCRIPTION
  Virus genome annotation package, last updated July 30th 2020""", file=sys.stderr)
        
    input_organism_li2 = [x[10:].split(",") for x in sys.argv[1:] if x[:10] == "-organism="]
    input_organism_li = list(set([y for x in input_organism_li2 for y in x]))  
    input_fasta_li = [f for f in [x[4:] for x in sys.argv[1:] if x[:4] == "-in="] if os.path.isfile(f) and ".fa" in f]
    output_fasta_li = [f for f in [x[5:] for x in sys.argv[1:] if x[:5] == "-out="] if os.path.isfile(f) and ".fa" in f]
    output_format = [x[12:] for x in sys.argv[1:] if x[:12] == "-out_format="]
    if len(input_fasta_li) > 0:
        if len(output_fasta_li) > 1: 
            print("Only one output file is supported. Please make revisions and try again.", file=sys.stderr)
        elif len(output_format) > 1:
            print("Only one output file format is supported. Please make revisions and try again.", file=sys.stderr)
        elif len([f for f in output_format if f not in valid_output_format_li]) == 1:
            print("Invalid output file format option provided. Supported options are "+\
                  readable_list(valid_output_format_li)+". Please make revisions and try again.", file=sys.stderr)
        else:
            read_fasta_li2 = [[], []]
            for fasta in input_fasta_li:
                temp_read_fasta_li2 = readFASTA(fasta)
                read_fasta_li2[0] += temp_read_fasta_li2[0]
                read_fasta_li2[1] += temp_read_fasta_li2[1]
            annotate_li4 = annotateFASTA(input_organism_li, read_fasta_li2)
            output_gff3 = outputGFF(read_fasta_li2, annotate_li4)
            if len(output_fasta_li) == 1:
                write_file = open(output_fasta_li[0], "w")
                write_file.write(output_gff3)
                write_file.close()
            else:
                print(output_gff3)

