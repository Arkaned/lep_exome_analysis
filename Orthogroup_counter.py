from Bio import SeqIO
from collections import Counter
import os

def parse_fasta_file(fasta_file):
    """ Parse a FASTA file and return a dictionary of sequences """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def calculate_amino_acid_content(sequence):
    """Calculate the amino acid content of a given protein sequence"""
    length = len(sequence)
    aa_counts = Counter(sequence)
    aa_content = {aa: count / length for aa, count in aa_counts.items()}
    return aa_content

def calculate_total_no_amino_acids(sequence):
    length = len(sequence)
    aa_counts = Counter(sequence)
    total_amino_acids = sum(aa_counts.values())
    return total_amino_acids
    
def process_directory(directory):
    """Process all Fasta files in a directory and calculate amino acid content."""
    aa_content_dict = {}
    for filename in os.listdir(directory):
        if filename.endswith(".faa") or filename.endswith("fa"):
            file_path = os.path.join(directory, filename)
            sequences = parse_fasta_file(file_path)
            for seq_id, seq in sequences.items():
                #aa_content_dict[seq_id] = calculate_amino_acid_content(seq) # use this if you want to calculate % for each amino acid
                
                aa_content_dict[seq_id] = calculate_total_no_amino_acids(seq) 
    return aa_content_dict

# example usage
directory_path = "/home/arkaned/tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/OrthoFinder/Results_May15/Orthogroup_Sequences/"
amino_acid_content = process_directory(directory_path)

sum_of_all_genes = 0
for seq_id, aa_content in amino_acid_content.items():
    print(f"Sequence ID: {seq_id}, {aa_content}") 
    sum_of_all_genes += aa_content
    ##for aa, content in aa_content: these two lines are useful if each amino acid's composition is wanted
    ##    print(f"  {aa}: {content:.4f}")     # if so, remove aa_content from above

print(" sum of all genes is: ", sum_of_all_genes)