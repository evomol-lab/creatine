import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

# ---------------------------------------------------------
# Configuração de Estilo Visual
# ---------------------------------------------------------
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.edgecolor'] = '#d0d0d0'
plt.rcParams['axes.linewidth'] = 1.1

COLOR_PALETTE = {
    'Benign': '#2ecc71',           # Verde
    'Likely benign': '#3498db',    # Azul
    'VUS / Uncertain': '#95a5a6',  # Cinza
    'Likely pathogenic': '#e67e22',# Laranja
    'Pathogenic': '#e74c3c'        # Vermelho
}

CLIN_ORDER = ['Benign', 'Likely benign', 'VUS / Uncertain', 'Likely pathogenic', 'Pathogenic']

# 1. Carregar variantes dbSNP humanas para P50440 (GATM)
df_dbsnp = pd.read_csv('P50440_missense_dbSNP_AAposAA.csv')
df_dbsnp['significancia_clinica'] = df_dbsnp['significancia_clinica'].fillna('Uncertain/Not provided')

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

df_dbsnp['categoria_clinica'] = df_dbsnp['significancia_clinica'].apply(simplify_clin)

# 2. Mapeamento UniProt P50440 -> Alinhamento FASTA
def map_uniprot_to_alignment(alignment_file, target_id='P50440'):
    records = list(SeqIO.parse(alignment_file, 'fasta'))
    ref_rec = None
    for rec in records:
        if rec.id == target_id:
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

def analyze_alignment(alignment_file, aln_label):
    records, ref_rec, pos_map = map_uniprot_to_alignment(alignment_file)
    other_records = [r for r in records if r.id != 'P50440' and 'HUMAN' not in r.id]
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
            'alinhamento': aln_label,
            'AAposAA': var_code,
            'posicao': pos,
            'wild': wild,
            'mut': mut,
            'dbSNP_id': rs_id,
            'significancia_clinica': clin_sig,
            'categoria_clinica': clin_cat,
            'num_organismos': mut_count,
            'num_total_outras_seqs': num_other_seqs,
            'freq_organismos': mut_freq,
            'observada_outros_org': mut_count > 0
        })
    return pd.DataFrame(results)

df_4668 = analyze_alignment('guest_4668_seq_aligned_trimmed.fasta', 'Alinhamento 4668 (415 seqs)')
df_4725 = analyze_alignment('guest_4725_seq_aligned_trimmed.fasta', 'Alinhamento 4725 (371 seqs)')

# Salvar tabelas processadas
df_4668.to_csv('variantes_humanas_gatm_em_outros_organismos_4668.csv', index=False)
df_4725.to_csv('variantes_humanas_gatm_em_outros_organismos_4725.csv', index=False)

# ---------------------------------------------------------
# CONSTRUÇÃO DA FIGURA DE PAINEL MULTIPLO (2x2)
# ---------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)
fig.suptitle('Presença de Variantes Missense Humanas (dbSNP P50440 - GATM) em Sequências de Outros Organismos', 
             fontsize=16, fontweight='bold', y=0.98)

# ---------------------------------------------------------
# PAINEL A: Proporção de Variantes Observadas por Categoria Clínica (Alinhamento 4668)
# ---------------------------------------------------------
ax_a = axes[0, 0]

cat_counts_4668 = []
for cat in CLIN_ORDER:
    sub = df_4668[df_4668['categoria_clinica'] == cat]
    total = len(sub)
    obs = sub['observada_outros_org'].sum()
    pct_obs = (obs / total * 100) if total > 0 else 0
    cat_counts_4668.append({'categoria': cat, 'total': total, 'obs': obs, 'pct_obs': pct_obs})

df_cat_4668 = pd.DataFrame(cat_counts_4668)

bars_obs_a = ax_a.bar(df_cat_4668['categoria'], df_cat_4668['pct_obs'], color=[COLOR_PALETTE[c] for c in CLIN_ORDER], alpha=0.85, edgecolor='black', width=0.55)

for bar, row in zip(bars_obs_a, df_cat_4668.itertuples()):
    yval = bar.get_height()
    ax_a.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{row.obs}/{row.total}\n({row.pct_obs:.1f}%)", 
              ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_a.set_title('A. Variantes Observadas em Outros Organismos (Alinhamento 4668 - 415 seqs)', fontsize=12, fontweight='bold', pad=10)
ax_a.set_ylabel('% de Variantes Humanas Observadas em Outros Organismos', fontsize=10, fontweight='bold')
ax_a.set_ylim(0, 100)
ax_a.tick_params(axis='x', rotation=15, labelsize=9)

# ---------------------------------------------------------
# PAINEL B: Proporção de Variantes Observadas por Categoria Clínica (Alinhamento 4725)
# ---------------------------------------------------------
ax_b = axes[0, 1]

cat_counts_4725 = []
for cat in CLIN_ORDER:
    sub = df_4725[df_4725['categoria_clinica'] == cat]
    total = len(sub)
    obs = sub['observada_outros_org'].sum()
    pct_obs = (obs / total * 100) if total > 0 else 0
    cat_counts_4725.append({'categoria': cat, 'total': total, 'obs': obs, 'pct_obs': pct_obs})

df_cat_4725 = pd.DataFrame(cat_counts_4725)

bars_obs_b = ax_b.bar(df_cat_4725['categoria'], df_cat_4725['pct_obs'], color=[COLOR_PALETTE[c] for c in CLIN_ORDER], alpha=0.85, edgecolor='black', width=0.55)

for bar, row in zip(bars_obs_b, df_cat_4725.itertuples()):
    yval = bar.get_height()
    ax_b.text(bar.get_x() + bar.get_width()/2.0, yval + 1.5, f"{row.obs}/{row.total}\n({row.pct_obs:.1f}%)", 
              ha='center', va='bottom', fontsize=9, fontweight='bold')

ax_b.set_title('B. Variantes Observadas em Outros Organismos (Alinhamento 4725 - 371 seqs)', fontsize=12, fontweight='bold', pad=10)
ax_b.set_ylabel('% de Variantes Humanas Observadas em Outros Organismos', fontsize=10, fontweight='bold')
ax_b.set_ylim(0, 100)
ax_b.tick_params(axis='x', rotation=15, labelsize=9)

# ---------------------------------------------------------
# PAINEL C: Distribuição ao Longo da Proteína P50440 (Posição vs Frequência)
# ---------------------------------------------------------
ax_c = axes[1, 0]

df_obs_4668 = df_4668[df_4668['observada_outros_org']].copy()

for cat in CLIN_ORDER:
    sub = df_obs_4668[df_obs_4668['categoria_clinica'] == cat]
    if not sub.empty:
        ax_c.scatter(sub['posicao'], sub['freq_organismos'] * 100, 
                     color=COLOR_PALETTE[cat], label=cat, s=sub['num_organismos']*0.5 + 25, 
                     alpha=0.75, edgecolors='black', linewidth=0.6)

# Destacar mutações patogênicas humanas presentes em outros organismos no GATM
patog_obs = df_obs_4668[df_obs_4668['categoria_clinica'].isin(['Pathogenic', 'Likely pathogenic'])]
for _, r in patog_obs.iterrows():
    ax_c.annotate(f"{r['AAposAA']} ({r['num_organismos']} orgs)", 
                  (r['posicao'], r['freq_organismos'] * 100),
                  textcoords="offset points", xytext=(0, 8), ha='center',
                  fontsize=8, fontweight='bold', color='#c0392b',
                  arrowprops=dict(arrowstyle="->", color='#c0392b', lw=0.8))

ax_c.set_title('C. Frequência de Ocorrência em Outras Espécies por Posição na Proteína P50440 (Aln 4668)', fontsize=12, fontweight='bold', pad=10)
ax_c.set_xlabel('Posição do Aminoácido na Proteína P50440 (Human GATM)', fontsize=10, fontweight='bold')
ax_c.set_ylabel('% de Espécies Não-Humanas com a Variante', fontsize=10, fontweight='bold')
ax_c.set_xlim(0, 430)
ax_c.set_ylim(-2, 100)
ax_c.legend(title='Categoria Clínica dbSNP', loc='upper right', frameon=True, fontsize=8, title_fontsize=9)

# ---------------------------------------------------------
# PAINEL D: Top Variantes Humanas Mais Frequentes em Outras Espécies
# ---------------------------------------------------------
ax_d = axes[1, 1]

# Selecionar as 12 variantes humanas com maior número de organismos no Alinhamento 4668
top_vars = df_4668.sort_values(by='num_organismos', ascending=False).head(12).iloc[::-1]

colors_top = [COLOR_PALETTE[cat] for cat in top_vars['categoria_clinica']]
bars_d = ax_d.barh(top_vars['AAposAA'], top_vars['num_organismos'], color=colors_top, alpha=0.85, edgecolor='black', height=0.6)

for bar, row in zip(bars_d, top_vars.itertuples()):
    xval = bar.get_width()
    ax_d.text(xval + 4, bar.get_y() + bar.get_height()/2.0, 
              f"{row.num_organismos} orgs ({row.freq_organismos*100:.1f}%) | {row.categoria_clinica}", 
              ha='left', va='center', fontsize=8, fontweight='bold')

ax_d.set_title('D. Top 12 Variantes Humanas Mais Abundantes em Outros Organismos (Aln 4668)', fontsize=12, fontweight='bold', pad=10)
ax_d.set_xlabel('Número de Organismos Não-Humanos com a Variante Solicitada', fontsize=10, fontweight='bold')
ax_d.set_xlim(0, 430)

# Ajuste fino do layout
plt.tight_layout(rect=[0, 0, 1, 0.96])

# Salvar figura em alta resolução
output_png = 'variantes_humanas_gatm_dbsnp_em_outros_organismos.png'
plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.close()

print(f"Figura gerada com sucesso: {output_png}")
