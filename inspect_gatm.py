import glob
import os

gatm_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /GATM"

fasta_files = glob.glob(os.path.join(gatm_dir, "*.fasta"))
tsv_files = glob.glob(os.path.join(gatm_dir, "*.tsv"))

print("FASTA files in GATM:")
for f in sorted(fasta_files):
    print(" ", os.path.basename(f))

print("\nTSV files in GATM:")
for t in sorted(tsv_files):
    print(" ", os.path.basename(t))

def check_gatm_fasta(fasta_file):
    headers = []
    p00390_seq = ""
    curr_header = ""
    curr_seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if curr_header:
                    if "P00390" in curr_header:
                        p00390_seq = "".join(curr_seq)
                curr_header = line
                headers.append(line)
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_header and "P00390" in curr_header:
            p00390_seq = "".join(curr_seq)
            
    print(f"\nFile: {os.path.basename(fasta_file)}")
    print(f" Total sequences: {len(headers)}")
    print(f" Contains P00390: {'P00390' in ''.join(headers)}")
    if p00390_seq:
        unaligned_len = len(p00390_seq.replace('-', ''))
        aligned_len = len(p00390_seq)
        print(f" P00390 aligned length: {aligned_len}, unaligned length: {unaligned_len}")

for f in sorted(fasta_files):
    check_gatm_fasta(f)
