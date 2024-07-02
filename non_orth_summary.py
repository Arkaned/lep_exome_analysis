import pandas as pd
import re
from collections import defaultdict
# the main objective of this programming file is to get the sum for species ideally...

df = pd.read_csv("/home/arkaned/tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/nonorth_seq_and_aa.csv")
df_n = df.values

aa_symbs = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def extract_count(value):
    match = re.search(r"\((\d+),", value)
    return int(match.group(1)) if match else 0

# Initialize a dictionary to store the total counts for each species
species_aa_counts = defaultdict(lambda: defaultdict(int))
#aa_counts = defaultdict(int)

# Iterate through each row in the DataFrame
for index, row in df.iterrows():
    species = row['species_file']
    
    # Sum the counts for the columns of interest
    for col in aa_symbs:  # Adjust this list if you have more columns
        count = extract_count(row[col])
        aa = col
        #species_counts[species] += count
        #aa_counts[aa] += count
        species_aa_counts[species][aa] += count
        
        
for species, aa_counts in species_aa_counts.items():
    print(f"Species: {species}")
    for aa, total_count in aa_counts.items():
        print(f"  Amino Acid: {aa}, Total Count: {total_count}")
#print("Total counts per species:")
#for species, total_count in species_counts.items():
#    print(f"Species: {species}, Total Amino Acid Count: {total_count}")

# Print the total counts for each amino acid
#print("\nTotal counts per amino acid:")
#for aa, total_count in aa_counts.items():
#    print(f"Amino Acid: {aa}, Total Count: {total_count}")
cheese = pd.DataFrame(species_aa_counts)
cheese['amino_acid'] = aa_symbs
cheese.to_csv("/home/arkaned/tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/nonorth_aa_counts_sp_summary.csv", index=False)
