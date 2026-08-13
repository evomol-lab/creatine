# Análise de Diversidade de Resíduos vs. Predições do AlphaMissense (SLC6A8 - P48029)

Esta análise mapeou cada coluna dos 6 alinhamentos `.fasta` da pasta `SLC6A8` para as posições da sequência de referência humana **P48029** (SLC6A8). Para cada posição, calculou-se a **diversidade de resíduos** (número de aminoácidos únicos no alinhamento, ignorando lacunas/gaps) e correlacionou-se esse valor com as contagens das classes de predição do **AlphaMissense** (*Pathogenic*, *Benign* e *Ambiguous*) obtidas a partir dos arquivos `.tsv`.

---

## 📊 Principais Descobertas Biológicas

1. **Alta Patogenicidade em Sítios Conservados (Menor Diversidade - $D=1$)**:
   - Os sítios com menor diversidade ($D=1$, 100% conservados no alinhamento) possuem uma proporção esmagadora de predições **Patogênicas** pelo AlphaMissense.
   - No alinhamento `guest_4682` (200 sequências), a média de predições patogênicas para posições com $D=1$ é de **18.03** (de 19 possíveis), com menos de **0.34** predições benignas por sítio.

2. **Aumento de Variantes Benignas em Sítios Diversos ($D > 1$)**:
   - À medida que a diversidade de resíduos aumenta ($D \ge 2$), a proporção de mutações classificadas como **Benignas** cresce drasticamente, enquanto as mutações patogênicas diminuem.

---

## 🖼️ Arquivos de Gráficos Gerados

Todos os gráficos foram gerados em alta resolução (300 DPI) e salvos na pasta `/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /SLC6A8/`:

1. **Painel Comparativo Geral**:
   - `SLC6A8_all_alignments_diversity_vs_alphamissense.png`: Compara a distribuição percentual das classes do AlphaMissense por nível de diversidade nos 6 alinhamentos side-by-side.

2. **Gráficos Detalhados Individuais** (Contagem Média, Proporção % e Track Sequencial de P48029):
   - `guest_4682_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4754_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4755_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4756_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4757_seq_aligned_trimmed_diversity_vs_alphamissense.png`
   - `guest_4758_seq_aligned_trimmed_diversity_vs_alphamissense.png`

---

## 🛠️ Script de Geração
O script em Python responsável pelo mapeamento das posições da P48029, cruzamento de dados e plotagem está salvo em:
`/media/jpmslima/home2/coding/creatine-1/plot_slc6a8_diversity_am.py`
