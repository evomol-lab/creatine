import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set plot styling
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#d0d0d0'
plt.rcParams['axes.linewidth'] = 1.1

base_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 "
gatm_dir = os.path.join(base_dir, "GATM")
slc6a8_dir = os.path.join(base_dir, "SLC6A8")

# Colors
colors = {
    'pathogenic': '#d62828', # Crimson Red
    'ambiguous': '#f4a261',  # Soft Amber/Orange
    'benign': '#2a9d8f'      # Teal Green
}

def parse_fasta_and_map_ref(fasta_file, target_id):
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
            
    ref_header = None
    for h in sequences:
        if target_id in h:
            ref_header = h
            break
            
    if not ref_header:
        raise ValueError(f"Target {target_id} not found in {fasta_file}")
        
    ref_seq = sequences[ref_header]
    seq_list = list(sequences.values())
    align_len = len(ref_seq)
    
    pos_map = {}
    ref_pos = 0
    
    for col_idx in range(align_len):
        char_ref = ref_seq[col_idx]
        if char_ref != '-':
            ref_pos += 1
            col_chars = [s[col_idx] for s in seq_list if col_idx < len(s) and s[col_idx] != '-']
            div = len(set(col_chars))
            pos_map[ref_pos] = div
            
    return pos_map, len(seq_list)

# Data processing for GATM (4668)
gatm_fasta = os.path.join(gatm_dir, "guest_4668_seq_aligned_trimmed.fasta")
gatm_am_file = os.path.join(gatm_dir, "AM_classification_counts (2).tsv")
pos_map_gatm, num_seqs_gatm = parse_fasta_and_map_ref(gatm_fasta, "P50440")
df_am_gatm = pd.read_csv(gatm_am_file, sep="\t")
df_div_gatm = pd.DataFrame(list(pos_map_gatm.items()), columns=["Position", "Residue_Diversity"])
merged_gatm = pd.merge(df_div_gatm, df_am_gatm, on="Position", how="inner")

max_div_gatm = merged_gatm["Residue_Diversity"].max()
if max_div_gatm > 8:
    merged_gatm["Div_Group"] = merged_gatm["Residue_Diversity"].apply(lambda x: str(x) if x <= 7 else "8+")
    div_order_gatm = [str(i) for i in range(1, 8)] + ["8+"]
else:
    merged_gatm["Div_Group"] = merged_gatm["Residue_Diversity"].astype(str)
    div_order_gatm = [str(i) for i in sorted(merged_gatm["Residue_Diversity"].unique())]

grouped_means_gatm = merged_gatm.groupby("Div_Group")[["pathogenic", "ambiguous", "benign"]].mean().reindex(div_order_gatm)
grouped_props_gatm = grouped_means_gatm.div(grouped_means_gatm.sum(axis=1), axis=0) * 100

# Data processing for SLC6A8 (4682)
slc6a8_fasta = os.path.join(slc6a8_dir, "guest_4682_seq_aligned_trimmed.fasta")
slc6a8_am_file = os.path.join(slc6a8_dir, "AM_classification_counts (1).tsv")
pos_map_slc, num_seqs_slc = parse_fasta_and_map_ref(slc6a8_fasta, "P48029")
df_am_slc = pd.read_csv(slc6a8_am_file, sep="\t")
df_div_slc = pd.DataFrame(list(pos_map_slc.items()), columns=["Position", "Residue_Diversity"])
merged_slc = pd.merge(df_div_slc, df_am_slc, on="Position", how="inner")

max_div_slc = merged_slc["Residue_Diversity"].max()
if max_div_slc > 8:
    merged_slc["Div_Group"] = merged_slc["Residue_Diversity"].apply(lambda x: str(x) if x <= 7 else "8+")
    div_order_slc = [str(i) for i in range(1, 8)] + ["8+"]
else:
    merged_slc["Div_Group"] = merged_slc["Residue_Diversity"].astype(str)
    div_order_slc = [str(i) for i in sorted(merged_slc["Residue_Diversity"].unique())]

grouped_means_slc = merged_slc.groupby("Div_Group")[["pathogenic", "ambiguous", "benign"]].mean().reindex(div_order_slc)
grouped_props_slc = grouped_means_slc.div(grouped_means_slc.sum(axis=1), axis=0) * 100


# ---------------------------------------------------------
# Build Combined 2x2 Panel Figure
# ---------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(18, 12), dpi=300)
fig.suptitle("Residue Diversity & AlphaMissense Pathogenicity Analysis\nGATM (Alignment 4668) vs. SLC6A8 (Alignment 4682)", 
             fontsize=16, fontweight='bold', y=0.98)

width = 0.65

# --- Panel A: GATM (4668) Percentage Distribution ---
ax_a = axes[0, 0]
x_pos_gatm = np.arange(len(div_order_gatm))
bottom_prop_gatm = np.zeros(len(div_order_gatm))

for cls in ['pathogenic', 'ambiguous', 'benign']:
    props = grouped_props_gatm[cls].values
    ax_a.bar(x_pos_gatm, props, width, bottom=bottom_prop_gatm, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
    for idx, val in enumerate(props):
        if val > 7:
            ax_a.text(x_pos_gatm[idx], bottom_prop_gatm[idx] + val/2, f"{val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    bottom_prop_gatm += props

ax_a.set_xticks(x_pos_gatm)
ax_a.set_xticklabels(div_order_gatm, fontsize=10, fontweight='bold')
ax_a.set_xlabel("Residue Diversity Level (Alignment 4668 - GATM)", fontsize=11, labelpad=6)
ax_a.set_ylabel("Percentage of Predictions (%)", fontsize=11, fontweight='bold')
ax_a.set_ylim(0, 105)
ax_a.set_title("A. Proportional Distribution of AlphaMissense Classes (GATM - 415 seqs)", fontsize=12, fontweight='bold', pad=10)
ax_a.legend(title="AlphaMissense Class", frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
ax_a.grid(axis='y', linestyle='--', alpha=0.5)

# --- Panel B: SLC6A8 (4682) Percentage Distribution ---
ax_b = axes[0, 1]
x_pos_slc = np.arange(len(div_order_slc))
bottom_prop_slc = np.zeros(len(div_order_slc))

for cls in ['pathogenic', 'ambiguous', 'benign']:
    props = grouped_props_slc[cls].values
    ax_b.bar(x_pos_slc, props, width, bottom=bottom_prop_slc, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
    for idx, val in enumerate(props):
        if val > 7:
            ax_b.text(x_pos_slc[idx], bottom_prop_slc[idx] + val/2, f"{val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    bottom_prop_slc += props

ax_b.set_xticks(x_pos_slc)
ax_b.set_xticklabels(div_order_slc, fontsize=10, fontweight='bold')
ax_b.set_xlabel("Residue Diversity Level (Alignment 4682 - SLC6A8)", fontsize=11, labelpad=6)
ax_b.set_ylabel("Percentage of Predictions (%)", fontsize=11, fontweight='bold')
ax_b.set_ylim(0, 105)
ax_b.set_title("B. Proportional Distribution of AlphaMissense Classes (SLC6A8 - 200 seqs)", fontsize=12, fontweight='bold', pad=10)
ax_b.legend(title="AlphaMissense Class", frameon=True, facecolor='white', framealpha=0.9, fontsize=9.5)
ax_b.grid(axis='y', linestyle='--', alpha=0.5)

# --- Panel C: GATM (4668) Sequence Track along P50440 ---
ax_c = axes[1, 0]
ax_c_twin = ax_c.twinx()

pos_g = merged_gatm["Position"].values
div_g = merged_gatm["Residue_Diversity"].values
patho_g = merged_gatm["pathogenic"].values

l1_g = ax_c.plot(pos_g, div_g, color='#457b9d', linewidth=1.2, label='Residue Diversity', alpha=0.85)
l2_g = ax_c_twin.plot(pos_g, patho_g, color='#d62828', linewidth=1.2, label='Pathogenic Mutations Count', alpha=0.75)

low_div_mask_g = (div_g == 1)
ax_c.scatter(pos_g[low_div_mask_g], div_g[low_div_mask_g], color='#1d3557', s=10, zorder=5, label='Conserved Sites (D=1)')

ax_c.set_xlabel("P50440 (GATM) Sequence Position", fontsize=11, fontweight='bold')
ax_c.set_ylabel("Residue Diversity", fontsize=11, color='#457b9d', fontweight='bold')
ax_c_twin.set_ylabel("Pathogenic Variants Count", fontsize=11, color='#d62828', fontweight='bold')
ax_c.tick_params(axis='y', labelcolor='#457b9d')
ax_c_twin.tick_params(axis='y', labelcolor='#d62828')
ax_c.set_title("C. Sequence Position Track: GATM (P50440)", fontsize=12, fontweight='bold', pad=10)

lines_g = l1_g + l2_g + [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1d3557', markersize=5, label='Conserved Sites (D=1)')]
labels_g = [l.get_label() for l in lines_g]
ax_c.legend(lines_g, labels_g, loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=8.5)
ax_c.grid(True, linestyle='--', alpha=0.4)

# --- Panel D: SLC6A8 (4682) Sequence Track along P48029 ---
ax_d = axes[1, 1]
ax_d_twin = ax_d.twinx()

pos_s = merged_slc["Position"].values
div_s = merged_slc["Residue_Diversity"].values
patho_s = merged_slc["pathogenic"].values

l1_s = ax_d.plot(pos_s, div_s, color='#457b9d', linewidth=1.2, label='Residue Diversity', alpha=0.85)
l2_s = ax_d_twin.plot(pos_s, patho_s, color='#d62828', linewidth=1.2, label='Pathogenic Mutations Count', alpha=0.75)

low_div_mask_s = (div_s == 1)
ax_d.scatter(pos_s[low_div_mask_s], div_s[low_div_mask_s], color='#1d3557', s=10, zorder=5, label='Conserved Sites (D=1)')

ax_d.set_xlabel("P48029 (SLC6A8) Sequence Position", fontsize=11, fontweight='bold')
ax_d.set_ylabel("Residue Diversity", fontsize=11, color='#457b9d', fontweight='bold')
ax_d_twin.set_ylabel("Pathogenic Variants Count", fontsize=11, color='#d62828', fontweight='bold')
ax_d.tick_params(axis='y', labelcolor='#457b9d')
ax_d_twin.tick_params(axis='y', labelcolor='#d62828')
ax_d.set_title("D. Sequence Position Track: SLC6A8 (P48029)", fontsize=12, fontweight='bold', pad=10)

lines_s = l1_s + l2_s + [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1d3557', markersize=5, label='Conserved Sites (D=1)')]
labels_s = [l.get_label() for l in lines_s]
ax_d.legend(lines_s, labels_s, loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=8.5)
ax_d.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout(rect=[0, 0, 1, 0.95])

output_panel_png = "/media/jpmslima/home2/coding/creatine-1/Evol/panel_plot_gatm4668_slc6a84682.png"
plt.savefig(output_panel_png, dpi=300, bbox_inches='tight')
plt.close()
print(f"Combined 2x2 panel plot saved successfully: {output_panel_png}")


# ---------------------------------------------------------
# Also generate 2-subplot individual figures (without 1st subplot) for each gene
# ---------------------------------------------------------

# 1. GATM 4668 2-panel figure
fig_g, axes_g = plt.subplots(2, 1, figsize=(12, 10), dpi=300, gridspec_kw={'height_ratios': [1, 1.2]})
fig_g.suptitle("AlphaMissense Pathogenicity & Diversity Analysis: GATM (Alignment 4668 - P50440)", fontsize=14, fontweight='bold')

# Top: Proportional Distribution %
ax_g1 = axes_g[0]
bottom_prop_g = np.zeros(len(div_order_gatm))
for cls in ['pathogenic', 'ambiguous', 'benign']:
    props = grouped_props_gatm[cls].values
    ax_g1.bar(x_pos_gatm, props, width, bottom=bottom_prop_g, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
    for idx, val in enumerate(props):
        if val > 8:
            ax_g1.text(x_pos_gatm[idx], bottom_prop_g[idx] + val/2, f"{val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=10)
    bottom_prop_g += props

ax_g1.set_xticks(x_pos_gatm)
ax_g1.set_xticklabels(div_order_gatm, fontsize=11, fontweight='bold')
ax_g1.set_xlabel("Residue Diversity Level", fontsize=11)
ax_g1.set_ylabel("Percentage of Predictions (%)", fontsize=11, fontweight='bold')
ax_g1.set_ylim(0, 105)
ax_g1.set_title("Proportional Distribution of Pathogenic, Ambiguous & Benign Predictions", fontsize=12, fontweight='bold')
ax_g1.legend(title="AlphaMissense Class", frameon=True, facecolor='white', framealpha=0.9, fontsize=10)
ax_g1.grid(axis='y', linestyle='--', alpha=0.5)

# Bottom: Sequence track along P50440
ax_g2 = axes_g[1]
ax_g2_twin = ax_g2.twinx()
l1 = ax_g2.plot(pos_g, div_g, color='#457b9d', linewidth=1.2, label='Residue Diversity', alpha=0.85)
l2 = ax_g2_twin.plot(pos_g, patho_g, color='#d62828', linewidth=1.2, label='Pathogenic Mutations Count', alpha=0.75)
ax_g2.scatter(pos_g[low_div_mask_g], div_g[low_div_mask_g], color='#1d3557', s=12, zorder=5, label='Conserved Sites (D=1)')

ax_g2.set_xlabel("P50440 (GATM) Sequence Position", fontsize=11, fontweight='bold')
ax_g2.set_ylabel("Residue Diversity", fontsize=11, color='#457b9d', fontweight='bold')
ax_g2_twin.set_ylabel("Pathogenic Variants Count", fontsize=11, color='#d62828', fontweight='bold')
ax_g2.tick_params(axis='y', labelcolor='#457b9d')
ax_g2_twin.tick_params(axis='y', labelcolor='#d62828')
ax_g2.set_title("Sequence Position Track: Residue Diversity vs. Pathogenic Prediction Density", fontsize=12, fontweight='bold')
lines = l1 + l2 + [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1d3557', markersize=6, label='Conserved Sites (D=1)')]
labels = [l.get_label() for l in lines]
ax_g2.legend(lines, labels, loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
ax_g2.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
gatm_2panel_out = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /GATM/guest_4668_seq_aligned_trimmed_diversity_vs_alphamissense_2panel.png"
plt.savefig(gatm_2panel_out, dpi=300, bbox_inches='tight')
plt.close()
print(f"GATM 2-panel plot saved: {gatm_2panel_out}")


# 2. SLC6A8 4682 2-panel figure
fig_s, axes_s = plt.subplots(2, 1, figsize=(12, 10), dpi=300, gridspec_kw={'height_ratios': [1, 1.2]})
fig_s.suptitle("AlphaMissense Pathogenicity & Diversity Analysis: SLC6A8 (Alignment 4682 - P48029)", fontsize=14, fontweight='bold')

# Top: Proportional Distribution %
ax_s1 = axes_s[0]
bottom_prop_s = np.zeros(len(div_order_slc))
for cls in ['pathogenic', 'ambiguous', 'benign']:
    props = grouped_props_slc[cls].values
    ax_s1.bar(x_pos_slc, props, width, bottom=bottom_prop_s, label=cls.capitalize(), color=colors[cls], edgecolor='white', alpha=0.9)
    for idx, val in enumerate(props):
        if val > 8:
            ax_s1.text(x_pos_slc[idx], bottom_prop_s[idx] + val/2, f"{val:.1f}%", ha='center', va='center', color='white', fontweight='bold', fontsize=10)
    bottom_prop_s += props

ax_s1.set_xticks(x_pos_slc)
ax_s1.set_xticklabels(div_order_slc, fontsize=11, fontweight='bold')
ax_s1.set_xlabel("Residue Diversity Level", fontsize=11)
ax_s1.set_ylabel("Percentage of Predictions (%)", fontsize=11, fontweight='bold')
ax_s1.set_ylim(0, 105)
ax_s1.set_title("Proportional Distribution of Pathogenic, Ambiguous & Benign Predictions", fontsize=12, fontweight='bold')
ax_s1.legend(title="AlphaMissense Class", frameon=True, facecolor='white', framealpha=0.9, fontsize=10)
ax_s1.grid(axis='y', linestyle='--', alpha=0.5)

# Bottom: Sequence track along P48029
ax_s2 = axes_s[1]
ax_s2_twin = ax_s2.twinx()
l1 = ax_s2.plot(pos_s, div_s, color='#457b9d', linewidth=1.2, label='Residue Diversity', alpha=0.85)
l2 = ax_s2_twin.plot(pos_s, patho_s, color='#d62828', linewidth=1.2, label='Pathogenic Mutations Count', alpha=0.75)
ax_s2.scatter(pos_s[low_div_mask_s], div_s[low_div_mask_s], color='#1d3557', s=12, zorder=5, label='Conserved Sites (D=1)')

ax_s2.set_xlabel("P48029 (SLC6A8) Sequence Position", fontsize=11, fontweight='bold')
ax_s2.set_ylabel("Residue Diversity", fontsize=11, color='#457b9d', fontweight='bold')
ax_s2_twin.set_ylabel("Pathogenic Variants Count", fontsize=11, color='#d62828', fontweight='bold')
ax_s2.tick_params(axis='y', labelcolor='#457b9d')
ax_s2_twin.tick_params(axis='y', labelcolor='#d62828')
ax_s2.set_title("Sequence Position Track: Residue Diversity vs. Pathogenic Prediction Density", fontsize=12, fontweight='bold')
lines = l1 + l2 + [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1d3557', markersize=6, label='Conserved Sites (D=1)')]
labels = [l.get_label() for l in lines]
ax_s2.legend(lines, labels, loc='upper right', frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
ax_s2.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout()
slc_2panel_out = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /SLC6A8/guest_4682_seq_aligned_trimmed_diversity_vs_alphamissense_2panel.png"
plt.savefig(slc_2panel_out, dpi=300, bbox_inches='tight')
plt.close()
print(f"SLC6A8 2-panel plot saved: {slc_2panel_out}")
