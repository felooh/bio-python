
##SPECIAL NOTES: DO NOT INCLUDE ANY FUNCTION CALLS IN YOUR FINAL CODE
## ALSO, USE THE FUNCTION NAMES EXACTLY AS I HAVE WRITTEN THEM
##  AND DO NOT CHANGE THE NUMBER OR TYPES OF PARAMETERS

# NOTE: All the functions below do not return values. Instead, they write files.
#       You can return 1 but it is not necessary.

#Global Variables used in functions 1 and 2 below:
standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

wm={"A":[1.5,0,-1.5,0,-1.5],
    "G":[-1.5,1,1,-1.5,0.5],
    "C":[-1.5,-1.5,-1.5,0,0],
    "T":[-1.5,-1.5,0,0.5,0]}

#Function 1: TRANSLATE Reads in a DNA sequence fasta file,
#            and outputs a Fasta file of protein translations.
#            This function uses the global variable standard_code above.
#            The output fasta will have the same titles as the input file.


def dna_2_prot(f1, f2="translated_fasta.txt"):
    """f1=input file name, f2=output file name"""
    return

#Example function call:

#dna_2_prot("testDNA.txt")









#Function 2: Reads in a fasta file of DNA sequences,and uses the weight matrix
#            dictionary (global variable wm - see above) to find
#            all the scores in each sequence that have a score greater than zero.

#The function needs to calculate the weight matrix scores for all overlapping
#sequences of a given length. In the wm case, the length is 5.
#For example: "AGTGTCA" has 3 overlapping sequences of 5 bases in length:
"""
AGTGTCA
AGTGT
 GTGTC
  TGTCA
"""

#The function will write the sequences and the scores to a file.
#The output file should look something like this. Note the spaces and commmas:
"""
Protein1: AAGGC,1.0 AGGCA,2.0 CGGAT,0.5 ...
Protein2: AAGGC,1.0 AGGCA,2.0 AAGGC,1.0 ...
...
"""

def weight_matrix_scores(f1, w=5, f2="wm_scores.txt"):
    """f1=input file name, w=weight matrix length, f2=output file name"""
    return

#Example function call:

#weight_matrix_scores("testDNA.txt")

#Function 3: Reads in a fasta file and searches for a given motif. Outputs a
#    file of the number of times that the motif was found in each sequence. The
#    output file will be tab-delimited.

#Here is a little code snippet that might be useful. 

"""
import re
Seq1="MNGMNNMNRVFRLVQRM"
y=re.findall("V[A-Z]R[ML]",Seq1)
print(y)
['VFRL', 'VQRM']
"""
#And the output should look like this
"""
SeqName     Motif       Hits            
Seq1        MN[A-Z]        2
Seq2        MN[A-Z]        0
Seq3        MN[A-Z]        1
""" 

def motif_finder(f1, motif, f2="motifs.txt"):
    """f1=input file name, motif=the motif pattern, f2=output file name"""
    return

#Example function call:

#motif_finder("testDNA.txt","MN[A-Z]")
