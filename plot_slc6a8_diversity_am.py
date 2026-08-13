import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# Set stylish plot theme
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

slc6a8_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /SLC6A8"
counts_file = os.path.join(slc6a8_dir, "AM_classification_counts (1).tsv")

df_am = pd.read_csv(counts_file, sep="\t")

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
            
    return pos_map, len(seq_list)

fasta_files = sorted(glob.glob(os.path.join(slc6a8_dir, "*.fasta")))

# Color palette
colors = {
    'pathogenic': '#d62828', # Crimson Red
    'ambiguous': '#f4a261',  # Soft Amber/Orange
    'benign': '#2a9d8f'      # Teal Green
}

# 1. Generate Individual Plots for each FASTA file
for f in fasta_files:
    fname = os.path.basename(f)
    base_name = fname.replace(".fasta", "")
    pos_map, num_seqs = parse_fasta_and_map_p48029(f)
    
    df_div = pd.DataFrame(list(pos_map.items()), columns=["Position", "Residue_Diversity"])
    merged = pd.merge(df_div, df_am, on="Position", how="inner")
    
    # Bin diversity if max diversity is high
    max_div = merged["Residue_Diversity"].max()
    if max_div > 8:
        # Group > 8 as '8+' or bins
        merged["Div_Group"] = merged["Residue_Diversity"].apply(lambda x: str(x) if x <= 7 else "8+")
        div_order = [str(i) for i in range(1, 8)] + ["8+"]
    else:
        merged["Div_Group"] = merged["Residue_Diversity"].astype(str)
        div_order = [str(i) for i in sorted(merged["Residue_Diversity"].unique())]
        
    grouped_counts = merged.groupby("Div_Group")[["pathogenic", "ambiguous", "benign"]].sum().reindex(div_order)
    grouped_means = merged.groupby("Div_Group")[["pathogenic", "ambiguous", "benign"]].mean().reindex(div_order)
    grouped_sites = merged.groupby("Div_Group")["Position"].count().reindex(div_order)
    
    grouped_props = grouped_means.div(grouped_means.sum(axis=1), axis=0) * 100

    fig, axes = plt.subplots(3, 1, figsize=(14, 15), gridspec_kw={'height_ratios': [1.2, 1, 1.2]})
    
    # --- Subplot 1: Mean Predictions per Site by Diversity ---
    ax1 = axes[0]
    x_pos = np.arange(len(div_order))
    width = 0.65
    
    bottom_val = np.zeros(len(div_order))
    for cls in ['pathogenic', 'ambiguous', 'benign']:
        vals = grouped_means[cls].values
        ax1.bar(x_pos, vals, width, bottom=bottom_val, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
        bottom_val += vals
        
    # Annotate site counts on top of bars
    for idx, (p, total_mean) in enumerate(zip(x_pos, bottom_val)):
        n_sites = grouped_sites.iloc[idx]
        if not np.isnan(total_mean) and total_mean > 0:
            ax1.text(p, total_mean + 0.3, f"n={n_sites}", ha='center', va='bottom', fontsize=10, fontweight='bold', color='#333333')

    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(div_order, fontsize=11, fontweight='bold')
    ax1.set_xlabel("Residue Diversity Level (Unique Amino Acids in Alignment)", fontsize=12, labelpad=8)
    ax1.set_ylabel("Mean Predictions per Site", fontsize=12)
    ax1.set_title(f"AlphaMissense Predictions vs. Residue Diversity\nAlignment: {fname} ({num_seqs} sequences)", fontsize=14, fontweight='bold', pad=12)
    ax1.legend(title="AlphaMissense Class", frameon=True, facecolor='white', framealpha=0.9, fontsize=11, title_fontsize=11)
    ax1.grid(axis='y', linestyle='--', alpha=0.5)
    
    # --- Subplot 2: Percentage (%) Distribution of AlphaMissense Classes ---
    ax2 = axes[1]
    bottom_prop = np.zeros(len(div_order))
    for cls in ['pathogenic', 'ambiguous', 'benign']:
        props = grouped_props[cls].values
        ax2.bar(x_pos, props, width, bottom=bottom_prop, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
        # Add percentage labels inside bars
        for idx, val in enumerate(props):
            if val > 8: # Only show label if bar segment is wide enough
                ax2.text(x_pos[idx], bottom_prop[idx] + val/2, f"{val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=10)
        bottom_prop += props
        
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(div_order, fontsize=11, fontweight='bold')
    ax2.set_xlabel("Residue Diversity Level", fontsize=12, labelpad=8)
    ax2.set_ylabel("Percentage of Predictions (%)", fontsize=12)
    ax2.set_ylim(0, 105)
    ax2.set_title("Proportional Distribution of Pathogenic, Ambiguous & Benign Predictions", fontsize=13, fontweight='bold')
    ax2.legend(title="AlphaMissense Class", frameon=True, facecolor='white', framealpha=0.9, fontsize=11)
    ax2.grid(axis='y', linestyle='--', alpha=0.5)

    # --- Subplot 3: Sequence Track along P48029 (Position vs Diversity & Pathogenic Count) ---
    ax3 = axes[2]
    ax3_twin = ax3.twinx()
    
    positions = merged["Position"].values
    diversity = merged["Residue_Diversity"].values
    patho_counts = merged["pathogenic"].values
    
    l1 = ax3.plot(positions, diversity, color='#457b9d', linewidth=1.2, label='Residue Diversity', alpha=0.85)
    l2 = ax3_twin.plot(positions, patho_counts, color='#d62828', linewidth=1.2, label='Pathogenic Mutations Count', alpha=0.75)
    
    # Highlight lowest diversity sites (Diversity == 1)
    low_div_mask = (diversity == 1)
    ax3.scatter(positions[low_div_mask], diversity[low_div_mask], color='#1d3557', s=12, zorder=5, label='Conserved Sites (Diversity = 1)')

    ax3.set_xlabel("P48029 Sequence Position", fontsize=12)
    ax3.set_ylabel("Residue Diversity", fontsize=12, color='#457b9d', fontweight='bold')
    ax3_twin.set_ylabel("Pathogenic Variants Count", fontsize=12, color='#d62828', fontweight='bold')
    ax3.tick_params(axis='y', labelcolor='#457b9d')
    ax3_twin.tick_params(axis='y', labelcolor='#d62828')
    ax3.set_title("Sequence Position Track: Residue Diversity vs. Pathogenic Prediction Density", fontsize=13, fontweight='bold')
    
    # Combined legend for Subplot 3
    lines = l1 + l2 + [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1d3557', markersize=6, label='Conserved Sites (D=1)')]
    labels = [l.get_label() for l in lines]
    ax3.legend(lines, labels, loc='upper right', frameon=True, facecolor='white', framealpha=0.9)
    ax3.grid(True, linestyle='--', alpha=0.4)

    plt.tight_layout()
    output_png = os.path.join(slc6a8_dir, f"{base_name}_diversity_vs_alphamissense.png")
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved individual plot: {output_png}")


# 2. Generate a 2x3 Grid Plot comparing ALL 6 FASTA files side-by-side
fig, axes = plt.subplots(2, 3, figsize=(20, 12), sharey=True)
axes = axes.flatten()

for idx, f in enumerate(fasta_files):
    fname = os.path.basename(f)
    base_name = fname.replace(".fasta", "")
    pos_map, num_seqs = parse_fasta_and_map_p48029(f)
    
    df_div = pd.DataFrame(list(pos_map.items()), columns=["Position", "Residue_Diversity"])
    merged = pd.merge(df_div, df_am, on="Position", how="inner")
    
    # Diversity grouping for comparison: 1, 2, 3, 4, 5+
    merged["Div_Group"] = merged["Residue_Diversity"].apply(lambda x: str(x) if x <= 4 else "5+")
    div_order = ["1", "2", "3", "4", "5+"]
    
    grouped_means = merged.groupby("Div_Group")[["pathogenic", "ambiguous", "benign"]].mean().reindex(div_order)
    grouped_props = grouped_means.div(grouped_means.sum(axis=1), axis=0) * 100
    grouped_sites = merged.groupby("Div_Group")["Position"].count().reindex(div_order)
    
    ax = axes[idx]
    x_pos = np.arange(len(div_order))
    width = 0.65
    
    bottom_prop = np.zeros(len(div_order))
    for cls in ['pathogenic', 'ambiguous', 'benign']:
        props = grouped_props[cls].values
        ax.bar(x_pos, props, width, bottom=bottom_prop, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
        bottom_prop += props
        
    # Annotate site counts on top
    for i_pos, n_sites in enumerate(grouped_sites):
        if not np.isnan(n_sites) and n_sites > 0:
            ax.text(x_pos[i_pos], 102, f"n={n_sites}", ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(div_order, fontsize=11, fontweight='bold')
    ax.set_xlabel("Residue Diversity Level", fontsize=11)
    if idx % 3 == 0:
        ax.set_ylabel("Proportion of Mutations (%)", fontsize=12, fontweight='bold')
    ax.set_ylim(0, 115)
    ax.set_title(f"{fname}\n({num_seqs} aligned sequences)", fontsize=12, fontweight='bold', pad=8)
    ax.grid(axis='y', linestyle='--', alpha=0.4)

# Create a single unified legend for the grid figure
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=3, frameon=True, facecolor='white', fontsize=13, title="AlphaMissense Classification", title_fontsize=13)
fig.suptitle("AlphaMissense Pathogenicity vs. Alignment Residue Diversity across all SLC6A8 Alignments (P48029 Reference)", fontsize=16, fontweight='bold', y=1.06)

plt.tight_layout()
grid_output = os.path.join(slc6a8_dir, "SLC6A8_all_alignments_diversity_vs_alphamissense.png")
plt.savefig(grid_output, dpi=300, bbox_inches='tight')
plt.close()
print(f"\nSaved comprehensive multi-alignment grid plot: {grid_output}")
