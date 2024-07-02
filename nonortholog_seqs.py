import pandas as pd
import os
import numpy as np
from Bio import SeqIO
from collections import Counter

filename = "/home/arkaned/tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/OrthoFinder/Results_May15/Orthogroups/Orthogroups_UnassignedGenes.tsv"

df = pd.read_csv(filename, sep='\t')

new_df = []

for i in range(len(df)):
    for j in df.iloc[i, 1:36]:
        if not pd.isna(j):
            new_df.append(j)
            

new_df = pd.DataFrame(new_df)

###########################################################
# now that we have our dataset of non-orthologues, let's get the genetic sequences from the primary transcripts...
# Since these are non-orthologous, each of the gene sequences will be unique to their species. 


def parse_fasta_file(fasta_file):
    """ Parse a FASTA file and return a dictionary of sequences """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

seqs_to_append = []
species_names = []
def find_genes(dir, file):
    filepath = os.path.join(dir, file)
    sequences = parse_fasta_file(filepath)
    for i in df.columns[1:36]:
        if i in file:
            for gene_id, gen_seq in sequences.items(): 
                for j in df[str(i)]:
                    if j == gene_id:
                        seqs_to_append.append(gen_seq)
                        species_names.append(file)
                        break 
                            
                            
directory = "/home/arkaned/tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/"
#df = df.to_numpy()

def collate_sequences(dir):
    counter = 0
    for file in os.listdir(directory):
        if file.endswith(".faa") or file.endswith("fa"):
            find_genes(dir, file)
            counter += 1
            print(counter, file)

collate_sequences(directory)

new_df['non_ortho_seqs'] = pd.DataFrame(seqs_to_append)
new_df['species_file'] = pd.DataFrame(species_names)
#new_df.to_csv("non_orthologues_with_seq.csv", index=False)

# cool, so all of this works so far... 
# now I guess we can calculate amino acid content for each species? 
# oh, but do we need to do this per species?
# that might require a bit of restructuring of the previous code, e.g. labelling each column.

# A tricky thing about the counter function is that it doesn't read things in a consistent order... 
# amino acid symbols: A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V
aa_symbs = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
aa_counts_percentage = {aa: [] for aa in aa_symbs}

for seq in new_df['non_ortho_seqs']:
    total_length = len(seq)
    for aa in aa_symbs:
        count = seq.count(aa)
        percentage = (count / total_length)
        
        aa_counts_percentage[aa].append((count, percentage))

aa_counts_df = pd.DataFrame(aa_counts_percentage)

result_df = pd.concat([new_df, aa_counts_df], axis=1)

result_df.to_csv("tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/nonorth_seq_and_aa.csv", index=False)

# Great, so a list of things that need to be done from here for non-orthologs (and orthologs as well):
# get perc of each aa for each gene (done!)
# get sum and avg of such values (both percentage and total count) for each species
# 
#
# finally sum of AA for all species orthologues vs non-orthologues

# When done, don't forget to confirm that the sum of orthologs and non-orthologs adds up to the entire exome.



