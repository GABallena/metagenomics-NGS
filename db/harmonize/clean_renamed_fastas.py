# Portfolio-safe copy: paths/identifiers generalized; databases not included.
import os
import re

def is_valid_seq(seq):
    return bool(re.search(r'[ACGTNacgtn]', seq))

def clean_seq(seq):
    return re.sub(r'[^ACGTNacgtn]', '', seq.upper())

def clean_fasta_streamed(input_path, output_path, log_path):
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as infile, \
         open(output_path, 'w') as outfile, \
         open(log_path, 'w') as logfile:

        header = None
        seq_lines = []

        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequence = clean_seq(''.join(seq_lines))
                    if sequence:
                        outfile.write(f"{header}\n")
                        for i in range(0, len(sequence), 70):
                            outfile.write(sequence[i:i+70] + '\n')
                    else:
                        logfile.write(f"Removed: {header[1:]}\n")
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # final entry
        if header:
            sequence = clean_seq(''.join(seq_lines))
            if sequence:
                outfile.write(f"{header}\n")
                for i in range(0, len(sequence), 70):
                    outfile.write(sequence[i:i+70] + '\n')
            else:
                logfile.write(f"Removed: {header[1:]}\n")

def main():
    for file in os.listdir():
        if file.endswith("_renamed.fasta"):
            out = file.replace("_renamed.fasta", "_cleaned.fasta")
            log = file.replace("_renamed.fasta", "_cleaned.log")
            print(f"Cleaning: {file}")
            clean_fasta_streamed(file, out, log)
    print("Done.")

if __name__ == "__main__":
    main()
