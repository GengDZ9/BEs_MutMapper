import tkinter as tk
from tkinter import messagebox
import pandas as pd
from itertools import combinations, product
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Define the base list
bases = ['A', 'T', 'C', 'G']

# Define the codon table
codon_table = {
    'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
    'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
    'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stop', 'TAG': 'Stop',
    'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stop', 'TGG': 'Trp',
    'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
    'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
    'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
    'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

# Define the mapping table from three-letter to single-letter amino acids
single_letter_amino_acid = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D',
    'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Stop': '*'
}

def generate_mutations(codon, editor_type):
    mutations = []
    for num_mutations in range(1, 4):
        for positions in combinations(range(3), num_mutations):
            new_codons = [codon]
            for pos in positions:
                if editor_type == 'ABE':
                    new_codons = [
                        new_codon[:pos] + ('G' if new_codon[pos] == 'A' else new_codon[pos]) + new_codon[pos + 1:]
                        for new_codon in new_codons
                    ]
                elif editor_type == 'CBE':
                    new_codons = [
                        new_codon[:pos] + ('T' if new_codon[pos] == 'C' else new_codon[pos]) + new_codon[pos + 1:]
                        for new_codon in new_codons
                    ]
                elif editor_type == 'DFBE':
                    new_codons = [
                        new_codon[:pos] + ('G' if new_codon[pos] == 'A' else 'T' if new_codon[pos] == 'C' else new_codon[pos]) + new_codon[pos + 1:]
                        for new_codon in new_codons
                    ]
            mutations.extend(new_codons)
    return mutations

def calculate():
    input_seq = entry.get()
    editor_type = editor_var.get()
    if len(input_seq) <= 51 and len(input_seq) % 3 == 0 and all(base in bases for base in input_seq):
        codons = [input_seq[i:i + 3] for i in range(0, len(input_seq), 3)]
        all_mutation_results = []
        all_mutated_codons_combinations = product(*[generate_mutations(codon, editor_type) for codon in codons])
        for mutated_codons in all_mutated_codons_combinations:
            original_aa_seq = ''.join([single_letter_amino_acid.get(codon_table.get(codon, 'Unknown'), 'X') for codon in codons])
            mutated_aa_seq = ""
            for mutated_codon in mutated_codons:
                amino_acid = single_letter_amino_acid.get(codon_table.get(mutated_codon, 'Unknown'), 'X')
                mutated_aa_seq += amino_acid
                if amino_acid == '*':
                    break
            total_mutsite_num = sum([sum([1 for i in range(3) if codon[i] != mutated_codon[i]]) for codon, mutated_codon in zip(codons, mutated_codons)])
            mutated_seq = ''.join(mutated_codons)

            # Calculate the number of amino acid changes
            min_len = min(len(original_aa_seq), len(mutated_aa_seq))
            aa_change_num = sum([1 for i in range(min_len) if original_aa_seq[i] != mutated_aa_seq[i]])
            aa_change_num += abs(len(original_aa_seq) - len(mutated_aa_seq))

            result = {
                'original_seq': input_seq,
                'original_aa_seq': original_aa_seq,
                'mutated_seq': mutated_seq,
                'mutated_aa_seq': mutated_aa_seq,
                'total_mutsite_num': total_mutsite_num,
                'aa_change_num': aa_change_num
            }
            all_mutation_results.append(result)
        all_mutation_results = [result for result in all_mutation_results if result['total_mutsite_num'] > 0]
        mutation_df = pd.DataFrame(all_mutation_results)
        mutation_df.to_excel('user_input_mutation_results.xlsx', index=False)
        messagebox.showinfo("Result", "The results have been saved to user_input_mutation_results.xlsx. Data graphs have also been generated.")

        # Draw the distribution histogram of aa_change_num
        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        # Adjust the layout of the bar chart by setting the align parameter to 'left''
        n, bins, patches = ax1.hist(mutation_df['aa_change_num'], bins=range(int(mutation_df['aa_change_num'].min()), int(mutation_df['aa_change_num'].max()) + 2), edgecolor='black', align='left')
        ax1.set_title('Distribution Histogram of Amino Acid Change Number')
        ax1.set_xlabel('Number of Amino Acid Changes')
        ax1.set_xticks(range(int(mutation_df['aa_change_num'].min()), int(mutation_df['aa_change_num'].max()) + 1))
        ax1.set_ylabel('Frequency')

        # Add bar chart value labels
        for rect, num in zip(patches, n):
            height = rect.get_height()
            ax1.text(rect.get_x() + rect.get_width()/2, height, str(int(num)), ha='center', va='bottom')

        # Draw the proportion chart of aa_change_num
        ax2 = fig.add_subplot(1, 2, 2)
        counts = mutation_df['aa_change_num'].value_counts()
        percentages = counts / counts.sum() * 100
        ax2.pie(percentages, labels=percentages.index, autopct='%1.1f%%')
        ax2.set_title('Proportion Chart of Amino Acid Change Number')

        plt.tight_layout()

        # Save the figure as a PDF with 300 PPI
        with PdfPages('amino_acid_mutation_graphs.pdf') as pdf:
            pdf.savefig(fig, dpi=300)

        plt.show()

    else:
        messagebox.showerror("Input Error", "The input sequence length must be a multiple of 3 and not exceed 51 bases, and can only contain A, T, C, G. Please re-enter.")

root = tk.Tk()
root.title("Sequence Input")

# Create a frame to place the base editor selection section
editor_frame = tk.Frame(root)
editor_frame.pack(pady=10)

editor_label = tk.Label(editor_frame, text="Select Base Editor:")
editor_label.pack(side=tk.LEFT, padx=10)

editor_var = tk.StringVar()
editor_var.set('DFBE')

abe_radio = tk.Radiobutton(editor_frame, text="ABE (A to G)", variable=editor_var, value='ABE', foreground='blue')
abe_radio.pack(side=tk.LEFT, padx=5)

cbe_radio = tk.Radiobutton(editor_frame, text="CBE (C to T)", variable=editor_var, value='CBE', foreground='green')
cbe_radio.pack(side=tk.LEFT, padx=5)

dfbe_radio = tk.Radiobutton(editor_frame, text="DFBE (A to G and C to T)", variable=editor_var, value='DFBE', foreground='orange')
dfbe_radio.pack(side=tk.LEFT, padx=5)

# Create a frame to place the sequence input section
sequence_frame = tk.Frame(root)
sequence_frame.pack(pady=20)

label = tk.Label(sequence_frame, text="Please enter a sequence with a length that is a multiple of 3 and does not exceed 51 bases:")
label.pack(pady=5)

entry = tk.Entry(sequence_frame)
entry.pack(pady=5)

button = tk.Button(root, text="Submit", command=calculate)
button.pack(pady=20)

root.mainloop()
    