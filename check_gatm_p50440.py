import glob
import os

gatm_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /GATM"
fasta_files = sorted(glob.glob(os.path.join(gatm_dir, "*.fasta")))

for f in fasta_files:
    fname = os.path.basename(f)
    print(f"\n=================== {fname} ===================")
    
    headers = []
    p50440_seq = ""
    curr_header = ""
    curr_seq = []
    
    with open(f, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if curr_header:
                    if "P50440" in curr_header:
                        p50440_seq = "".join(curr_seq)
                curr_header = line
                headers.append(line)
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_header and "P50440" in curr_header:
            p50440_seq = "".join(curr_seq)
            
    print(f" Total sequences: {len(headers)}")
    print(f" Contains P50440: {'P50440' in ''.join(headers)}")
    if p50440_seq:
        unaligned_len = len(p50440_seq.replace('-', ''))
        aligned_len = len(p50440_seq)
        print(f" P50440 aligned length: {aligned_len}, unaligned length: {unaligned_len}")
    else:
        # Print sample headers to see how UniProt IDs are formatted
        print(" Sample headers:")
        for h in headers[:5]:
            print("  ", h)
