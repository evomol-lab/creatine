import glob
import os
import collections

slc6a8_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /SLC6A8"
fasta_files = sorted(glob.glob(os.path.join(slc6a8_dir, "*.fasta")))

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
            
    # Find P48029 header
    p48029_header = None
    for h in sequences:
        if "P48029" in h:
            p48029_header = h
            break
            
    if not p48029_header:
        print(f"P48029 not found in {fasta_file}")
        return None
        
    p48029_seq = sequences[p48029_header]
    seq_list = list(sequences.values())
    num_seqs = len(seq_list)
    align_len = len(p48029_seq)
    
    pos_map = {} # p48029_pos (1-indexed) -> diversity
    p48029_pos = 0
    
    for col_idx in range(align_len):
        char_p48029 = p48029_seq[col_idx]
        if char_p48029 != '-':
            p48029_pos += 1
            col_chars = [s[col_idx] for s in seq_list if col_idx < len(s) and s[col_idx] != '-']
            div = len(set(col_chars))
            pos_map[p48029_pos] = div
            
    return pos_map

for f in fasta_files:
    pos_map = parse_fasta_and_map_p48029(f)
    if pos_map:
        div_counts = collections.Counter(pos_map.values())
        print(f"\n{os.path.basename(f)} (Mapped P48029 positions: {len(pos_map)}):")
        for d in sorted(div_counts.keys()):
            print(f"  Diversity {d}: {div_counts[d]} sites")
