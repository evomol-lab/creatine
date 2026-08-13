import numpy as np
import matplotlib.pyplot as plt
import collections
import math
import glob
import os

def read_fasta(filename):
    sequences = []
    current_seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append("".join(current_seq))
    return sequences

def calculate_metrics(sequences):
    num_seqs = len(sequences)
    seq_len = len(sequences[0])
    entropies = []
    diversities = []
    
    for i in range(seq_len):
        col = [seq[i] for seq in sequences if seq[i] != '-']
        
        # Diversity
        unique_residues = set(col)
        diversities.append(len(unique_residues))
        
        # Entropy
        counts = collections.Counter(col)
        total = len(col)
        entropy = 0
        if total > 0:
            for count in counts.values():
                p = count / total
                entropy -= p * math.log2(p)
        entropies.append(entropy)
        
    return entropies, diversities

def plot_metrics(fasta_file):
    print(f"Processing {fasta_file}...")
    sequences = read_fasta(fasta_file)
    if not sequences:
        print("No sequences found.")
        return
        
    seq_len = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_len:
            print(f"Error: Sequences are not of the same length in {fasta_file}.")
            return
            
    entropies, diversities = calculate_metrics(sequences)
    positions = np.arange(1, seq_len + 1)
    
    # Plot Entropy
    plt.figure(figsize=(15, 5))
    plt.plot(positions, entropies, color='blue', linewidth=1)
    plt.title(f"Shannon Entropy - {os.path.basename(fasta_file)}")
    plt.xlabel("Alignment Position")
    plt.ylabel("Entropy (bits)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(fasta_file.replace('.fasta', '_entropy.png'), dpi=300)
    plt.close()
    
    # Plot Diversity
    plt.figure(figsize=(15, 5))
    plt.plot(positions, diversities, color='red', linewidth=1)
    plt.title(f"Residue Diversity - {os.path.basename(fasta_file)}")
    plt.xlabel("Alignment Position")
    plt.ylabel("Number of Unique Residues")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(fasta_file.replace('.fasta', '_diversity.png'), dpi=300)
    plt.close()
    
    print(f"Saved plots for {fasta_file}")

def main():
    base_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 "
    fasta_files = glob.glob(os.path.join(base_dir, "**/*.fasta"), recursive=True)
    
    if not fasta_files:
        print(f"No fasta files found in {base_dir}")
        return
        
    for fasta in fasta_files:
        plot_metrics(fasta)

if __name__ == '__main__':
    main()
