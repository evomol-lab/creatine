import requests
import pandas as pd

# Mapeamento para código de 3 letras (opcional)
AA_3LETTER = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'
}

def obter_variacoes_dbsnp_uniprot(uniprot_id="P48029"):
    """
    Obtém as variações missense do UniProt/dbSNP para uma proteína específica
    e gera arquivos de saída com as mutações no formato AAposAA (ex: M1V, R390W).
    """
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{uniprot_id}"
    headers = {"Accept": "application/json"}
    
    print(f"Buscando variações missense na UniProt API para {uniprot_id}...")
    response = requests.get(url, headers=headers)
    
    if not response.ok:
        print(f"Erro ao acessar API: {response.status_code}")
        return
        
    data = response.json()
    features = data.get("features", [])
    records = []
    
    for f in features:
        wild = f.get("wildType")
        mut = f.get("alternativeSequence")
        pos = f.get("begin")
        
        # Filtra apenas substituições missense (1 aminoácido)
        if wild and mut and len(wild) == 1 and len(mut) == 1 and wild != mut:
            # Extrai os IDs do dbSNP (rsIDs)
            xrefs = f.get("xrefs", [])
            dbsnp_ids = [x["id"] for x in xrefs if x.get("name", "").lower() == "dbsnp"]
            
            # Formato 1 letra: AAposAA (ex: M1V)
            aa_pos_aa = f"{wild}{pos}{mut}"
            
            # Formato 3 letras: AAposAA (ex: Met1Val)
            wild_3 = AA_3LETTER.get(wild, wild)
            mut_3 = AA_3LETTER.get(mut, mut)
            aa_pos_aa_3 = f"{wild_3}{pos}{mut_3}"
            
            # Significância clínica (ClinVar / UniProt)
            clin_sig = []
            if f.get("clinicalSignificances"):
                clin_sig = [cs.get("type", "") for cs in f.get("clinicalSignificances", [])]
            
            records.append({
                "uniprot_id": uniprot_id,
                "gene": data.get("geneName", "SLC6A8"),
                "posicao": int(pos),
                "aa_wild": wild,
                "aa_mut": mut,
                "AAposAA": aa_pos_aa,
                "AAposAA_3letras": aa_pos_aa_3,
                "dbSNP_id": ", ".join(dbsnp_ids) if dbsnp_ids else "N/A",
                "significancia_clinica": ", ".join(clin_sig) if clin_sig else "N/A"
            })
            
    df = pd.DataFrame(records)
    
    # 1. Filtra apenas as que possuem ID no dbSNP (rsID)
    df_dbsnp = df[df["dbSNP_id"] != "N/A"].copy()
    
    # Salva arquivo CSV completo com metadados e dbSNP ID
    file_csv = f"{uniprot_id}_missense_dbSNP_AAposAA.csv"
    df_dbsnp.to_csv(file_csv, index=False)
    
    # 2. Salva uma lista simples em .txt apenas com os códigos no formato AAposAA
    file_txt = f"{uniprot_id}_missense_AAposAA_list.txt"
    with open(file_txt, "w") as f_out:
        for val in df_dbsnp["AAposAA"]:
            f_out.write(f"{val}\n")
            
    print(f"\n--- Processamento Concluído ---")
    print(f"Total de mutações missense encontradas: {len(df)}")
    print(f"Total de mutações missense associadas ao dbSNP: {len(df_dbsnp)}")
    print(f"Arquivo CSV gerado: {file_csv}")
    print(f"Arquivo TXT (lista simples) gerado: {file_txt}\n")
    print("Exemplo das 10 primeiras linhas:")
    print(df_dbsnp[["AAposAA", "AAposAA_3letras", "dbSNP_id", "significancia_clinica"]].head(10))

if __name__ == "__main__":
    obter_variacoes_dbsnp_uniprot("P48029")