import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# ---------------------------------------------------------
# Visual Styling Configuration
# ---------------------------------------------------------
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#d0d0d0'
plt.rcParams['axes.linewidth'] = 1.1

COLOR_PALETTE = {
    'Benign': '#2ecc71',           # Emerald Green
    'Likely benign': '#3498db',    # Soft Blue
    'VUS / Uncertain': '#95a5a6',  # Neutral Gray
    'Likely pathogenic': '#e67e22',# Warm Orange
    'Pathogenic': '#e74c3c'        # Crimson Red
}

CLIN_ORDER = ['Benign', 'Likely benign', 'VUS / Uncertain', 'Likely pathogenic', 'Pathogenic']

def simplify_clin(val):
    val_str = str(val)
    if 'Pathogenic' in val_str and 'Likely' not in val_str:
        return 'Pathogenic'
    elif 'Likely pathogenic' in val_str:
        return 'Likely pathogenic'
    elif 'Benign' in val_str and 'Likely' not in val_str:
        return 'Benign'
    elif 'Likely benign' in val_str:
        return 'Likely benign'
    else:
        return 'VUS / Uncertain'

base_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 "
gatm_dir = os.path.join(base_dir, "GATM")
slc6a8_dir = os.path.join(base_dir, "SLC6A8")

# 1. Map UniProt to FASTA Alignment
def map_uniprot_to_alignment(alignment_file, target_id):
    records = list(SeqIO.parse(alignment_file, 'fasta'))
    ref_rec = None
    for rec in records:
        if target_id in rec.id:
            ref_rec = rec
            break
    if not ref_rec:
        raise ValueError(f'Target {target_id} not found in {alignment_file}')
        
    pos_map = {}
    u_pos = 0
    for col_idx, char in enumerate(str(ref_rec.seq)):
        if char != '-' and char != '.':
            u_pos += 1
            pos_map[u_pos] = col_idx
            
    return records, ref_rec, pos_map

def analyze_dbsnp_alignment(csv_file, alignment_file, target_id, aln_label):
    df_dbsnp = pd.read_csv(csv_file)
    df_dbsnp['significancia_clinica'] = df_dbsnp['significancia_clinica'].fillna('Uncertain/Not provided')
    df_dbsnp['categoria_clinica'] = df_dbsnp['significancia_clinica'].apply(simplify_clin)
    
    records, ref_rec, pos_map = map_uniprot_to_alignment(alignment_file, target_id)
    other_records = [r for r in records if target_id not in r.id and 'HUMAN' not in r.id]
    num_other_seqs = len(other_records)
    
    results = []
    for idx, row in df_dbsnp.iterrows():
        pos = int(row['posicao'])
        wild = str(row['aa_wild'])
        mut = str(row['aa_mut'])
        var_code = str(row['AAposAA'])
        clin_sig = str(row['significancia_clinica'])
        clin_cat = str(row['categoria_clinica'])
        rs_id = str(row['dbSNP_id'])
        
        if pos not in pos_map:
            continue
            
        col_idx = pos_map[pos]
        other_aas = [str(r.seq)[col_idx] for r in other_records]
        
        mut_count = other_aas.count(mut)
        mut_freq = mut_count / num_other_seqs if num_other_seqs > 0 else 0.0
        
        results.append({
            'alignment': aln_label,
            'AAposAA': var_code,
            'position': pos,
            'wild': wild,
            'mut': mut,
            'dbSNP_id': rs_id,
            'clinical_significance': clin_sig,
            'clinical_category': clin_cat,
            'num_organisms': mut_count,
            'num_total_other_seqs': num_other_seqs,
            'freq_organisms': mut_freq,
            'observed_in_other_orgs': mut_count > 0
        })
    return pd.DataFrame(results)

# Analyze GATM 4668
gatm_csv = os.path.join(gatm_dir, "P50440_missense_dbSNP_AAposAA.csv")
gatm_fasta = os.path.join(gatm_dir, "guest_4668_seq_aligned_trimmed.fasta")
df_gatm = analyze_dbsnp_alignment(gatm_csv, gatm_fasta, "P50440", "GATM 4668 (415 seqs)")

# Analyze SLC6A8 4682
slc6a8_csv = os.path.join(slc6a8_dir, "P48029_missense_dbSNP_AAposAA.csv")
slc6a8_fasta = os.path.join(slc6a8_dir, "guest_4682_seq_aligned_trimmed.fasta")
df_slc = analyze_dbsnp_alignment(slc6a8_csv, slc6a8_fasta, "P48029", "SLC6A8 4682 (200 seqs)")


# ---------------------------------------------------------
# Build Combined 2x2 Panel Figure (English)
# ---------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)
fig.suptitle('Natural Occurrence of Human Missense Variants (dbSNP) in Non-Human Orthologs\nComparing GATM (Alignment 4668) and SLC6A8 (Alignment 4682)', 
             fontsize=16, fontweight='bold', y=0.98)

# --- PANEL A: GATM 4668 Bar Chart (% Observed by Clinical Category) ---
ax_a = axes[0, 0]
cat_counts_gatm = []
for cat in CLIN_ORDER:
    sub = df_gatm[df_gatm['clinical_category'] == cat]
    total = len(sub)
    obs = sub['observed_in_other_orgs'].sum()
    pct_obs = (obs / total * 100) if total > 0 else 0
    cat_counts_gatm.append({'category': cat, 'total': total, 'obs': obs, 'pct_obs': pct_obs})

df_cat_gatm = pd.DataFrame(cat_counts_gatm)
bars_a = ax_a.bar(df_cat_gatm['category'], df_cat_gatm['pct_obs'], 
                  color=[COLOR_PALETTE[c] for c in CLIN_ORDER], alpha=0.85, edgecolor='black', width=0.55)

for bar, row in zip(bars_a, df_cat_gatm.itertuples()):
    yval = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{row.obs}/{row.total}\n({row.pct_obs:.1f}%)", 
              ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_a.set_title('A. Human Variants Observed in Non-Human Species (GATM 4668 — 415 seqs)', fontsize=11.5, fontweight='bold', pad=10)
ax_a.set_ylabel('% Human Missense Variants Observed', fontsize=10.5, fontweight='bold')
ax_a.set_ylim(0, 105)
ax_a.tick_params(axis='x', rotation=15, labelsize=9)
ax_a.grid(axis='y', linestyle='--', alpha=0.5)

# --- PANEL B: SLC6A8 4682 Bar Chart (% Observed by Clinical Category) ---
ax_b = axes[0, 1]
cat_counts_slc = []
for cat in CLIN_ORDER:
    sub = df_slc[df_slc['clinical_category'] == cat]
    total = len(sub)
    obs = sub['observed_in_other_orgs'].sum()
    pct_obs = (obs / total * 100) if total > 0 else 0
    cat_counts_slc.append({'category': cat, 'total': total, 'obs': obs, 'pct_obs': pct_obs})

df_cat_slc = pd.DataFrame(cat_counts_slc)
bars_b = ax_b.bar(df_cat_slc['category'], df_cat_slc['pct_obs'], 
                  color=[COLOR_PALETTE[c] for c in CLIN_ORDER], alpha=0.85, edgecolor='black', width=0.55)

for bar, row in zip(bars_b, df_cat_slc.itertuples()):
    yval = bar.get_height()
    ax_b.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{row.obs}/{row.total}\n({row.pct_obs:.1f}%)", 
              ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_b.set_title('B. Human Variants Observed in Non-Human Species (SLC6A8 4682 — 200 seqs)', fontsize=11.5, fontweight='bold', pad=10)
ax_b.set_ylabel('% Human Missense Variants Observed', fontsize=10.5, fontweight='bold')
ax_b.set_ylim(0, 105)
ax_b.tick_params(axis='x', rotation=15, labelsize=9)
ax_b.grid(axis='y', linestyle='--', alpha=0.5)

# --- PANEL C: GATM 4668 Scatter Track along P50440 ---
ax_c = axes[1, 0]
df_obs_gatm = df_gatm[df_gatm['observed_in_other_orgs']].copy()

for cat in CLIN_ORDER:
    sub = df_obs_gatm[df_obs_gatm['clinical_category'] == cat]
    if not sub.empty:
        ax_c.scatter(sub['position'], sub['freq_organisms'] * 100, 
                     color=COLOR_PALETTE[cat], label=cat, s=sub['num_organisms']*0.5 + 25, 
                     alpha=0.75, edgecolors='black', linewidth=0.6)

patog_obs_gatm = df_obs_gatm[df_obs_gatm['clinical_category'].isin(['Pathogenic', 'Likely pathogenic'])]
for _, r in patog_obs_gatm.iterrows():
    ax_c.annotate(f"{r['AAposAA']} ({r['num_organisms']} orgs)", 
                  (r['position'], r['freq_organisms'] * 100),
                  textcoords="offset points", xytext=(0, 8), ha='center',
                  fontsize=8, fontweight='bold', color='#c0392b',
                  arrowprops=dict(arrowstyle="->", color='#c0392b', lw=0.8))

ax_c.set_title('C. Occurrence Frequency across Non-Human Species along GATM (P50440)', fontsize=11.5, fontweight='bold', pad=10)
ax_c.set_xlabel('Amino Acid Position on P50440 (Human GATM)', fontsize=10, fontweight='bold')
ax_c.set_ylabel('% Non-Human Species Carrying Variant', fontsize=10, fontweight='bold')
ax_c.set_xlim(0, 430)
ax_c.set_ylim(-2, 100)
ax_c.legend(title='dbSNP Clinical Category', loc='upper right', frameon=True, fontsize=8, title_fontsize=8.5)
ax_c.grid(True, linestyle='--', alpha=0.4)

# --- PANEL D: SLC6A8 4682 Scatter Track along P48029 ---
ax_d = axes[1, 1]
df_obs_slc = df_slc[df_slc['observed_in_other_orgs']].copy()

for cat in CLIN_ORDER:
    sub = df_obs_slc[df_obs_slc['clinical_category'] == cat]
    if not sub.empty:
        ax_d.scatter(sub['position'], sub['freq_organisms'] * 100, 
                     color=COLOR_PALETTE[cat], label=cat, s=sub['num_organisms']*0.8 + 25, 
                     alpha=0.75, edgecolors='black', linewidth=0.6)

patog_obs_slc = df_obs_slc[df_obs_slc['clinical_category'].isin(['Pathogenic', 'Likely pathogenic'])]
for _, r in patog_obs_slc.iterrows():
    ax_d.annotate(f"{r['AAposAA']} ({r['num_organisms']} orgs)", 
                  (r['position'], r['freq_organisms'] * 100),
                  textcoords="offset points", xytext=(0, 8), ha='center',
                  fontsize=8, fontweight='bold', color='#c0392b',
                  arrowprops=dict(arrowstyle="->", color='#c0392b', lw=0.8))

ax_d.set_title('D. Occurrence Frequency across Non-Human Species along SLC6A8 (P48029)', fontsize=11.5, fontweight='bold', pad=10)
ax_d.set_xlabel('Amino Acid Position on P48029 (Human SLC6A8)', fontsize=10, fontweight='bold')
ax_d.set_ylabel('% Non-Human Species Carrying Variant', fontsize=10, fontweight='bold')
ax_d.set_xlim(0, 640)
ax_d.set_ylim(-2, 85)
ax_d.legend(title='dbSNP Clinical Category', loc='upper right', frameon=True, fontsize=8, title_fontsize=8.5)
ax_d.grid(True, linestyle='--', alpha=0.4)

plt.tight_layout(rect=[0, 0, 1, 0.95])

output_panel_2x2 = "/media/jpmslima/home2/coding/creatine-1/Evol/panel_dbsnp_cross_species_4668_4682_english.png"
plt.savefig(output_panel_2x2, dpi=300, bbox_inches='tight')
plt.close()
print(f"Combined 2x2 dbSNP panel plot saved: {output_panel_2x2}")


# ---------------------------------------------------------
# Build Comprehensive 2x3 Panel Figure (English) including Top Abundant Variants
# ---------------------------------------------------------
fig_2x3, axes_3 = plt.subplots(3, 2, figsize=(16, 17), dpi=300)
fig_2x3.suptitle('Cross-Species Conservation of Human Missense Variants (dbSNP)\nComparative Analysis: GATM (4668 Alignment) vs. SLC6A8 (4682 Alignment)', 
                 fontsize=16, fontweight='bold', y=0.99)

# Row 1: Panels A & B (Percentage Observed)
ax1_a = axes_3[0, 0]
bars_1a = ax1_a.bar(df_cat_gatm['category'], df_cat_gatm['pct_obs'], color=[COLOR_PALETTE[c] for c in CLIN_ORDER], alpha=0.85, edgecolor='black', width=0.55)
for bar, row in zip(bars_1a, df_cat_gatm.itertuples()):
    yval = bar.get_height()
    ax1_a.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{row.obs}/{row.total}\n({row.pct_obs:.1f}%)", ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax1_a.set_title('A. Human Variants Observed in Non-Human Species (GATM 4668 — 415 seqs)', fontsize=11, fontweight='bold', pad=8)
ax1_a.set_ylabel('% Human Variants Observed', fontsize=10, fontweight='bold')
ax1_a.set_ylim(0, 105)
ax1_a.tick_params(axis='x', rotation=15, labelsize=8.5)

ax1_b = axes_3[0, 1]
bars_1b = ax1_b.bar(df_cat_slc['category'], df_cat_slc['pct_obs'], color=[COLOR_PALETTE[c] for c in CLIN_ORDER], alpha=0.85, edgecolor='black', width=0.55)
for bar, row in zip(bars_1b, df_cat_slc.itertuples()):
    yval = bar.get_height()
    ax1_b.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{row.obs}/{row.total}\n({row.pct_obs:.1f}%)", ha='center', va='bottom', fontsize=8.5, fontweight='bold')
ax1_b.set_title('B. Human Variants Observed in Non-Human Species (SLC6A8 4682 — 200 seqs)', fontsize=11, fontweight='bold', pad=8)
ax1_b.set_ylabel('% Human Variants Observed', fontsize=10, fontweight='bold')
ax1_b.set_ylim(0, 105)
ax1_b.tick_params(axis='x', rotation=15, labelsize=8.5)

# Row 2: Panels C & D (Scatter Track along Protein)
ax2_c = axes_3[1, 0]
for cat in CLIN_ORDER:
    sub = df_obs_gatm[df_obs_gatm['clinical_category'] == cat]
    if not sub.empty:
        ax2_c.scatter(sub['position'], sub['freq_organisms'] * 100, color=COLOR_PALETTE[cat], label=cat, s=sub['num_organisms']*0.5 + 25, alpha=0.75, edgecolors='black', linewidth=0.6)
for _, r in patog_obs_gatm.iterrows():
    ax2_c.annotate(f"{r['AAposAA']} ({r['num_organisms']} orgs)", (r['position'], r['freq_organisms'] * 100), textcoords="offset points", xytext=(0, 8), ha='center', fontsize=7.5, fontweight='bold', color='#c0392b', arrowprops=dict(arrowstyle="->", color='#c0392b', lw=0.7))
ax2_c.set_title('C. Frequency across Non-Human Species along GATM (P50440)', fontsize=11, fontweight='bold', pad=8)
ax2_c.set_xlabel('Amino Acid Position on P50440 (GATM)', fontsize=9.5, fontweight='bold')
ax2_c.set_ylabel('% Non-Human Species Carrying Variant', fontsize=9.5, fontweight='bold')
ax2_c.set_xlim(0, 430)
ax2_c.set_ylim(-2, 100)
ax2_c.legend(title='dbSNP Category', loc='upper right', frameon=True, fontsize=7.5, title_fontsize=8)

ax2_d = axes_3[1, 1]
for cat in CLIN_ORDER:
    sub = df_obs_slc[df_obs_slc['clinical_category'] == cat]
    if not sub.empty:
        ax2_d.scatter(sub['position'], sub['freq_organisms'] * 100, color=COLOR_PALETTE[cat], label=cat, s=sub['num_organisms']*0.8 + 25, alpha=0.75, edgecolors='black', linewidth=0.6)
for _, r in patog_obs_slc.iterrows():
    ax2_d.annotate(f"{r['AAposAA']} ({r['num_organisms']} orgs)", (r['position'], r['freq_organisms'] * 100), textcoords="offset points", xytext=(0, 8), ha='center', fontsize=7.5, fontweight='bold', color='#c0392b', arrowprops=dict(arrowstyle="->", color='#c0392b', lw=0.7))
ax2_d.set_title('D. Frequency across Non-Human Species along SLC6A8 (P48029)', fontsize=11, fontweight='bold', pad=8)
ax2_d.set_xlabel('Amino Acid Position on P48029 (SLC6A8)', fontsize=9.5, fontweight='bold')
ax2_d.set_ylabel('% Non-Human Species Carrying Variant', fontsize=9.5, fontweight='bold')
ax2_d.set_xlim(0, 640)
ax2_d.set_ylim(-2, 85)
ax2_d.legend(title='dbSNP Category', loc='upper right', frameon=True, fontsize=7.5, title_fontsize=8)

# Row 3: Panels E & F (Top 10 Most Abundant Variants in Non-Human Species)
ax3_e = axes_3[2, 0]
top_gatm = df_gatm.sort_values(by='num_organisms', ascending=False).head(10).iloc[::-1]
colors_top_g = [COLOR_PALETTE[cat] for cat in top_gatm['clinical_category']]
bars_e = ax3_e.barh(top_gatm['AAposAA'], top_gatm['num_organisms'], color=colors_top_g, alpha=0.85, edgecolor='black', height=0.6)
for bar, row in zip(bars_e, top_gatm.itertuples()):
    xval = bar.get_width()
    ax3_e.text(xval + 3, bar.get_y() + bar.get_height()/2.0, f"{row.num_organisms} species ({row.freq_organisms*100:.1f}%) | {row.clinical_category}", ha='left', va='center', fontsize=7.5, fontweight='bold')
ax3_e.set_title('E. Top 10 Most Abundant Human Variants in Non-Human Species (GATM)', fontsize=11, fontweight='bold', pad=8)
ax3_e.set_xlabel('Number of Non-Human Species Carrying the Variant', fontsize=9.5, fontweight='bold')
ax3_e.set_xlim(0, 430)

ax3_f = axes_3[2, 1]
top_slc = df_slc.sort_values(by='num_organisms', ascending=False).head(10).iloc[::-1]
colors_top_s = [COLOR_PALETTE[cat] for cat in top_slc['clinical_category']]
bars_f = ax3_f.barh(top_slc['AAposAA'], top_slc['num_organisms'], color=colors_top_s, alpha=0.85, edgecolor='black', height=0.6)
for bar, row in zip(bars_f, top_slc.itertuples()):
    xval = bar.get_width()
    ax3_f.text(xval + 2, bar.get_y() + bar.get_height()/2.0, f"{row.num_organisms} species ({row.freq_organisms*100:.1f}%) | {row.clinical_category}", ha='left', va='center', fontsize=7.5, fontweight='bold')
ax3_f.set_title('F. Top 10 Most Abundant Human Variants in Non-Human Species (SLC6A8)', fontsize=11, fontweight='bold', pad=8)
ax3_f.set_xlabel('Number of Non-Human Species Carrying the Variant', fontsize=9.5, fontweight='bold')
ax3_f.set_xlim(0, 180)

plt.tight_layout(rect=[0, 0, 1, 0.97])

output_panel_2x3 = "/media/jpmslima/home2/coding/creatine-1/Evol/panel_dbsnp_cross_species_4668_4682_2x3_english.png"
plt.savefig(output_panel_2x3, dpi=300, bbox_inches='tight')
plt.close()
print(f"Comprehensive 2x3 dbSNP panel plot saved: {output_panel_2x3}")
