import pandas as pd
import math

filename = "/home/arkaned/tyota__ju_lle/orthofinder_tutorial/Exomes/primary_transcripts/OrthoFinder/Results_May15/Orthogroups/Orthogroups_UnassignedGenes.tsv"

df = pd.read_csv(filename, sep='\t')

new_df = []

for i in range(len(df)):
    for j in df.iloc[i, 1:36]:
        if not pd.isna(j):
            new_df.append(j)
            

print(new_df)  
print("length of non-orthologues: ", len(new_df))          
            