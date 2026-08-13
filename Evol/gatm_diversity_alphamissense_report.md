# Análise de Diversidade de Resíduos vs. Predições do AlphaMissense (GATM - P50440)

Esta análise mapeou cada coluna dos 4 alinhamentos `.fasta` da pasta `GATM` para as posições da sequência de referência humana **P50440** (GATM). Para cada posição, calculou-se a **diversidade de resíduos** (número de aminoácidos únicos no alinhamento, ignorando lacunas/gaps) e correlacionou-se esse valor com as contagens das classes de predição do **AlphaMissense** (*Pathogenic*, *Benign* e *Ambiguous*) dos arquivos `.tsv`.

---

## 📊 Principais Descobertas Biológicas em GATM

1. **Forte Associação entre Conservação Evolutiva ($D=1$) e Patogenicidade**:
   - Assim como em SLC6A8, os sítios totalmente conservados no alinhamento ($D=1$) concentram a maior quantidade de predições **Patogênicas** do AlphaMissense (média de **~12 a 13.7** predições patogênicas por sítio de 19 possíveis).

2. **Crescimento de Variantes Benignas em Posições Variáveis ($D > 1$)**:
   - Conforme a diversidade de aminoácidos aumenta nas colunas do alinhamento ($D \ge 2$), a proporção de mutações classificadas como **Benignas** aumenta progressivamente (passando de ~20-30% em $D=1$ para mais de 50-60% em sítios mais diversos).

---

## 🖼️ Arquivos de Gráficos Gerados

Todos os gráficos foram gerados em alta resolução (300 DPI) e salvos na pasta `/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /GATM/`:

1. **Painel Comparativo Geral**:
   - `GATM_all_alignments_diversity_vs_alphamissense.png`: Compara a distribuição percentual das classes do AlphaMissense por nível de diversidade nos 4 alinhamentos side-by-side (Grid 2x2).

2. **Gráficos Detalhados Individuais** (Média de Predições por Sítio, Proporção % e Track Sequencial de P50440):
   - `guest_4660_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4668_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4700_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4725_seq_aligned_trimmed_diversity_vs_alphamissense.png`

---

## 🛠️ Script de Geração
O script em Python responsável pelo mapeamento das posições da P50440, cruzamento de dados e plotagem para a proteína GATM foi salvo em:
`/media/jpmslima/home2/coding/creatine-1/plot_gatm_diversity_am.py`
