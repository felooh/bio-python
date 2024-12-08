# Kindly note that the testDNA.txt file should be in the same folder/directory as this file for it to work successfully.

import os
import re

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


# Function 1: Translate DNA to protein
def dna_2_prot(f1, f2="translated_fasta.txt"):
    """
    Translate a DNA file (FASTA format) to a protein file.

    Args:
        f1 (str): The name of the input DNA FASTA file.
        f2 (str): The name of the output protein FASTA file. Defaults to "translated_fasta.txt".

    This function reads a DNA sequence file in FASTA format, translates each sequence
    into its corresponding protein sequence, and writes the result to a new FASTA file.
    """

    # The input file name
    filename = f1

    # Get the directory of the current script
    this_dir = os.path.dirname(__file__)

    # Generate the full file path for the input file
    file_path = os.path.join(this_dir, filename)

    try:
        # Open the input DNA file for reading and the output file for writing
        with open(file_path, 'r') as infile, open(f2, 'w') as outfile:
            current_header = None  # Variable to store the current FASTA header
            current_sequence = ""  # Variable to accumulate DNA sequences under the current header

            for line in infile:
                line = line.strip()  # Remove any leading or trailing whitespace
                if line.startswith(">"):
                    # If the line starts with ">", it's a FASTA header
                    if current_header and current_sequence:
                        # Translate the accumulated DNA sequence to protein if there's a header and sequence
                        protein = translate_dna(current_sequence)
                        # Write the translated protein sequence to the output file in FASTA format
                        outfile.write(f"{current_header}\n{protein}\n")
                    
                    # Update the current header and reset the sequence accumulator
                    current_header = line
                    current_sequence = ""
                else:
                    # Otherwise, it's a DNA sequence line
                    # Append the uppercased sequence to the accumulator
                    current_sequence += line.upper()

            # After the loop, check if there's any remaining sequence to process
            if current_header and current_sequence:
                protein = translate_dna(current_sequence)
                outfile.write(f"{current_header}\n{protein}\n")

        # Success message
        print(f"✅ Translated protein file saved to: {f2}")

    except FileNotFoundError:
        # Handle the case where the input file does not exist
        print(f"❌ Error: The file '{file_path}' was not found. Please check the file name.")
    except Exception as e:
        # Catch any other exceptions and print the error message
        print(f"❌ An unexpected error occurred: {e}")

def translate_dna(sequence):
    """Translates a DNA sequence to a protein sequence using the standard genetic code."""
    protein = []
    for i in range(0, len(sequence) - 2, 3):
        # Replace thymine with uracil
        codon = sequence[i:i + 3].replace("T", "U")
        # Use '?' for unknown codons  
        amino_acid = standard_code.get(codon, "?")  
        # Stop translation at stop codons
        if amino_acid == "*":  
            break
        protein.append(amino_acid)
    return "".join(protein)



# Function 2: Calculate weight matrix scores
def weight_matrix_scores(f1, w=5, f2="wm_scores.txt"):
    """
    Calculate weight matrix scores for DNA sequences in a FASTA file.

    Args:
        f1 (str): The name of the input DNA FASTA file.
        w (int): The length of the window for scoring subsequences. Default is 5.
        f2 (str): The name of the output file to save scores. Default is "wm_scores.txt".

    This function computes scores for subsequences of length `w` in each DNA sequence 
    using a predefined weight matrix (`wm`) and writes the results to a file.
    """

    filename = f1

    # Get the directory of the current script
    this_dir = os.path.dirname(__file__)

    # Generate the full file path for the input file
    file_path = os.path.join(this_dir, filename)

    try:
        # Open the input DNA file for reading and the output file for writing
        with open(file_path, 'r') as infile, open(f2, 'w') as outfile:
            current_header = None  # Variable to store the current FASTA header
            current_sequence = ""  # Variable to accumulate DNA sequences under the current header

            for line in infile:
                line = line.strip()  # Remove any leading or trailing whitespace
                if line.startswith(">"):
                    # If the line starts with ">", it's a FASTA header
                    if current_header and current_sequence:
                        # Compute scores for all subsequences of length `w`
                        scores = []
                        for i in range(len(current_sequence) - w + 1):
                            subseq = current_sequence[i:i + w]
                            # Calculate the score for the current subsequence using the weight matrix
                            score = sum(
                                wm[base][j] if base in wm and j < len(wm[base]) else 0
                                for j, base in enumerate(subseq)
                            )
                            if score > 0:
                                scores.append(f"{subseq},{score:.1f}")
                        # Write scores to the output file if there are any
                        if scores:
                            outfile.write(f"{current_header} {' '.join(scores)}\n")

                    # Update the current header and reset the sequence accumulator
                    current_header = line
                    current_sequence = ""
                else:
                    # Otherwise, it's a DNA sequence line
                    # Append the uppercased sequence to the accumulator
                    current_sequence += line.upper()

            # After the loop, process any remaining sequence
            if current_header and current_sequence:
                scores = []
                for i in range(len(current_sequence) - w + 1):
                    subseq = current_sequence[i:i + w]
                    score = sum(
                        wm[base][j] if base in wm and j < len(wm[base]) else 0
                        for j, base in enumerate(subseq)
                    )
                    if score > 0:
                        scores.append(f"{subseq},{score:.1f}")
                if scores:
                    outfile.write(f"{current_header} {' '.join(scores)}\n")

        # Success message
        print(f"✅ Weight matrix scores written to: {f2}")

    except FileNotFoundError:
        # Handle the case where the input file does not exist
        print(f"❌ Error: The file '{file_path}' was not found. Please check the file name.")
    except Exception as e:
        # Catch any other exceptions and print the error message
        print(f"❌ An unexpected error occurred: {e}")
      

# Function 3: Find motifs
def motif_finder(f1, motif="MN[A-Z]", f2="motifs.txt"):
    """
    Find motifs in a DNA file using regular expressions.

    Args:
        f1 (str): The name of the input DNA FASTA file.
        motif (str): The regular expression pattern to search for. Default is "MN[A-Z]".
        f2 (str): The name of the output file to save motif search results. Default is "motifs.txt".

    This function identifies occurrences of the specified motif in each DNA sequence
    in the input FASTA file and writes the results to an output file.
    """
    filename = f1

    # Get the directory of the current script
    this_dir = os.path.dirname(__file__)

    # Generate the full file path for the input file
    file_path = os.path.join(this_dir, filename)

    try:
        # Open the input DNA file for reading and the output file for writing
        with open(file_path, 'r') as infile, open(f2, 'w') as outfile:
            # Write a header line to the output file
            outfile.write("SeqName\t\tMotif\t\tHits\n")

            current_header = None  # Variable to store the current FASTA header
            current_sequence = ""  # Variable to accumulate DNA sequences under the current header

            for line in infile:
                line = line.strip()  # Remove any leading or trailing whitespace
                if line.startswith(">"):
                    # If the line starts with ">", it's a FASTA header
                    if current_header and current_sequence:
                        # Search for motifs in the accumulated DNA sequence
                        hits = re.findall(motif, current_sequence)
                        # Write the results to the output file
                        outfile.write(f"{current_header[1:]}\t{motif}\t\t{len(hits)}\n")

                    # Update the current header and reset the sequence accumulator
                    current_header = line
                    current_sequence = ""
                else:
                    # Otherwise, it's a DNA sequence line
                    # Append the uppercased sequence to the accumulator
                    current_sequence += line.upper()

            # After the loop, process any remaining sequence
            if current_header and current_sequence:
                hits = re.findall(motif, current_sequence)
                outfile.write(f"{current_header[1:]}\t{motif}\t\t{len(hits)}\n")

        # Success message
        print(f"✅ Motif search results written to: {f2}")

    except FileNotFoundError:
        # Handle the case where the input file does not exist
        print(f"❌ Error: The file '{file_path}' was not found. Please check the file name.")
    except Exception as e:
        # Catch any other exceptions and print the error message
        print(f"❌ An unexpected error occurred: {e}")


