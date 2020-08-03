"""
Created on Mar 18, 2019
@author: blake.inderski
"""

import pandas as pd
import re
import os
import subprocess
from sequence_annotate import *
import sys
import traceback

################################################################################
################################################################################

dir_path = os.path.dirname(os.path.realpath(__file__))
py_file = os.path.basename(__file__)

df = pd.read_excel(dir_path + "/prrsv_orf5_cut_pattern_determination.xlsx")
crop_head_li = [i[i.index("= ")+2:].replace("^", "").replace("(", "[").replace(")", "]") for i in list(df) if "Unnamed" not in i]
cut_site_li = [1, 3, 4]
rflp_di = {i:[[]] for i in crop_head_li}
df = df.rename(columns=df.iloc[0]).drop(df.index[0])
cut_sites_li = df["Cut Sites"].values.tolist()
for x in range(len(cut_sites_li)):
    for y in range(len(cut_sites_li[x])):
        if str(cut_sites_li[x][y]).lower() != "nan" and str(cut_sites_li[x][y]).lower() != "-":
            rflp_di[crop_head_li[y]].append([int(i) for i in str(cut_sites_li[x][y]).split(",")])
        
################################################################################

def findRFLP(nucleotide_sequence, annotate_li2, require_full=True):
    #assigns restriction fragment length polymorphism (RFLP) pattern values to type 2 PRRSV ORF5 genes
    def matchRFLP(dna_sequence, rflp_value):
        for r in range(len(crop_head_li)):
            regex = crop_head_li[r]
            #find restriction sites
            finditer_li = re.finditer(r"(?=("+regex+"))", dna_sequence.replace("U", "T"))
            result = [i.start(1) for i in finditer_li]
            fix_result = [i+cut_site_li[r]+1 for i in result]
            if fix_result in rflp_di[regex]:
                rflp_value[r] = str(rflp_di[regex].index(fix_result)+1)
        return(rflp_value)
    rflp_value = ["null"] #does not contain full length ORF5
    orf5_li = [x for x in annotate_li2 if "ORF5" in x[0]]
    if orf5_li != [] and (require_full != True or "full" in orf5_li[0][0]):
        start, stop = [int(x) for x in orf5_li[0][1].split(",")]
        dna_sequence = nucleotide_sequence[start-1:stop]
        rflp_value = matchRFLP(dna_sequence, ["X" for i in range(3)])
        if "X" in "".join(rflp_value) or len(nucleotide_sequence) != 603:
            annotate = Annotate("PRRSV2")
            hmm_build = dir_path+"/required/PRRSV2/hmm_profiles/hmm_build/ORF5.hmm"    
            hmmalign_simple = annotate.hmmalign(dna_sequence, hmm_build)
            hmmalign_simple = re.sub("[a-z]+", "", hmmalign_simple)
            rflp_value = matchRFLP(hmmalign_simple, ["X" for i in range(3)])
    return("-".join(rflp_value))

################################################################################
################################################################################

#gather system arguments
if "prrsv_orf5_RFLP" in sys.argv[0]:
    if "-h" in sys.argv[1:] or "-help" in sys.argv[1:]:
        print("""USAGE
  python3 prrsv_orf5_RFLP.py input > output

OPTIONAL ARGUMENTS
 -h or -help
   Print USAGE, DESCRIPTION and ARGUMENTS
 -in=<Input_File>
   Input FASTA file name
 -out=<Output_File>
   Output FASTA file name
 -require_full=<Bool>
   Only assign RFLP pattern to input sequences that are complete (default=True)
 
DESCRIPTION
  Virus genome annotation package, last updated July 30th 2020""", file=sys.stderr)
        
    input_fasta_li = [f for f in [x[4:] for x in sys.argv[1:] if x[:4] == "-in="] if os.path.isfile(f) and ".fa" in f]
    output_fasta_li = [f for f in [x[5:] for x in sys.argv[1:] if x[:5] == "-out="] if os.path.isfile(f) and ".fa" in f]
    require_full = [y[0] for y in [re.findall("^\-require_full\=", "", x) for x in sys.argv[1:]] if y != []]
    if len(input_fasta_li) > 0:
        if len(output_fasta_li) > 1: 
            print("Only one output file is supported. Please make revisions and try again.", file=sys.stderr)
        else:
            read_fasta_li2 = [[], []]
            for fasta in input_fasta_li:
                temp_read_fasta_li2 = readFASTA(fasta)
                read_fasta_li2[0] += temp_read_fasta_li2[0]
                read_fasta_li2[1] += temp_read_fasta_li2[1]
            annotate_li4 = annotateFASTA(["PRRSV1", "PRRSV2"], read_fasta_li2, ["ORF5"]) 
            require_full = True
            if len(require_full) == 1 and require_full[0] == "f":
                require_full = False
            for fasta_index, annotate_li3 in enumerate(annotate_li4):
                if len(annotate_li3) > 0 and annotate_li3[0] == "PRRSV2":
                    read_fasta_li2[0][fasta_index] += "/"+findRFLP(read_fasta_li2[1][fasta_index], annotate_li3[1], require_full)
                else:
                    read_fasta_li2[0][fasta_index] += "/na" #invalid organism
                    
            if len(output_fasta_li) == 1:
                writeFASTA(output_fasta_li[0], read_fasta_li2)
            else:
                for fasta_index, head in enumerate(read_fasta_li2[0]):
                    print(">"+head)
                    print(read_fasta_li2[1][fasta_index])
                
            