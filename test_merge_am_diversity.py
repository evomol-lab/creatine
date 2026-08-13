import glob
import os
import pandas as pd
import numpy as np

slc6a8_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /SLC6A8"
counts_file = os.path.join(slc6a8_dir, "AM_classification_counts (1).tsv")

df_am = pd.read_csv(counts_file, sep="\t")
print("AlphaMissense Counts DF shape:", df_am.shape)
print(df_am.head())

# Read FASTA diversity maps
def parse_fasta_and_map_p48029(fasta_file):
    sequences = {}
    current_header = ""
    current_seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = "".join(current_seq)
            
    p48029_header = None
    for h in sequences:
        if "P48029" in h:
            p48029_header = h
            break
            
    p48029_seq = sequences[p48029_header]
    seq_list = list(sequences.values())
    align_len = len(p48029_seq)
    
    pos_map = {}
    p48029_pos = 0
    
    for col_idx in range(align_len):
        char_p48029 = p48029_seq[col_idx]
        if char_p48029 != '-':
            p48029_pos += 1
            col_chars = [s[col_idx] for s in seq_list if col_idx < len(s) and s[col_idx] != '-']
            div = len(set(col_chars))
            pos_map[p48029_pos] = div
            
    return pos_map

fasta_files = sorted(glob.glob(os.path.join(slc6a8_dir, "*.fasta")))

for f in fasta_files:
    fname = os.path.basename(f)
    pos_map = parse_fasta_and_map_p48029(f)
    df_div = pd.DataFrame(list(pos_map.items()), columns=["Position", "Residue_Diversity"])
    merged = pd.merge(df_div, df_am, on="Position", how="inner")
    
    print(f"\n=================== {fname} ===================")
    print(f"Total merged positions: {len(merged)}")
    
    # Group by Diversity == 1 vs > 1
    low_div = merged[merged["Residue_Diversity"] == 1]
    high_div = merged[merged["Residue_Diversity"] > 1]
    
    print("Diversity == 1 (Lowest diversity / Conserved):")
    print("  Sites count:", len(low_div))
    print("  Pathogenic sum:", low_div["pathogenic"].sum(), "avg:", low_div["pathogenic"].mean())
    print("  Benign sum:    ", low_div["benign"].sum(), "avg:", low_div["benign"].mean())
    print("  Ambiguous sum: ", low_div["ambiguous"].sum(), "avg:", low_div["ambiguous"].mean())
    
    print("Diversity > 1:")
    print("  Sites count:", len(high_div))
    print("  Pathogenic sum:", high_div["pathogenic"].sum(), "avg:", high_div["pathogenic"].mean())
    print("  Benign sum:    ", high_div["benign"].sum(), "avg:", high_div["benign"].mean())
    print("  Ambiguous sum: ", high_div["ambiguous"].sum(), "avg:", high_div["ambiguous"].mean())
