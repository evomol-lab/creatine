import glob
import os

slc6a8_dir = "/media/jpmslima/home2/coding/creatine-1/Evol/GATM - SLC6A8 /SLC6A8"

fasta_files = glob.glob(os.path.join(slc6a8_dir, "*.fasta"))
tsv_files = glob.glob(os.path.join(slc6a8_dir, "*.tsv"))

print("FASTA files:")
for f in sorted(fasta_files):
    print(" ", os.path.basename(f))

print("\nTSV files:")
for t in sorted(tsv_files):
    print(" ", os.path.basename(t))

def check_p48029(fasta_file):
    headers = []
    has_p48029 = False
    p48029_seq = ""
    curr_header = ""
    curr_seq = []
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if curr_header:
                    if "P48029" in curr_header:
                        p48029_seq = "".join(curr_seq)
                curr_header = line
                headers.append(line)
                curr_seq = []
            else:
                curr_seq.append(line)
        if curr_header and "P48029" in curr_header:
            p48029_seq = "".join(curr_seq)
            
    print(f"\nFile: {os.path.basename(fasta_file)}")
    print(f" Total sequences: {len(headers)}")
    print(f" Contains P48029: {'P48029' in ''.join(headers)}")
    if p48029_seq:
        unaligned_len = len(p48029_seq.replace('-', ''))
        aligned_len = len(p48029_seq)
        print(f" P48029 aligned length: {aligned_len}, unaligned length: {unaligned_len}")

for f in sorted(fasta_files):
    check_p48029(f)
