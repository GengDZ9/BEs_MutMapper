import pandas as pd
from itertools import combinations


bases = ['A', 'T', 'C', 'G']

codon_library = []

for first in bases:
    for second in bases:
        for third in bases:
            codon = first + second + third
            codon_library.append(codon)


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


stop_codons = {'TAA', 'TAG', 'TGA'}


mutation_results = []
for codon in codon_library:
    original_amino_acid = codon_table.get(codon, 'Unknown')
    original_is_stop = codon in stop_codons
    for num_mutations in range(1, 4):
        for positions in combinations(range(3), num_mutations):
            new_codons = [codon]
            for pos in positions:
                new_codons = [
                    new_codon[:pos] + ('G' if new_codon[pos] == 'A' else 'T' if new_codon[pos] == 'C' else new_codon[pos]) + new_codon[pos + 1:]
                    for new_codon in new_codons
                ]
            for mutated_codon in new_codons:
                mutated_amino_acid = codon_table.get(mutated_codon, 'Unknown')
                mutated_is_stop = mutated_codon in stop_codons
                result = {
                    'ori_codon': codon,
                    'ori_aa': original_amino_acid,
                    'mut_conda': mutated_codon,
                    'mut_aa': mutated_amino_acid,
                    'aa_change': 'yes' if original_amino_acid != mutated_amino_acid else 'no',
                    'mutsite_num': num_mutations,
                    'new_stopconda': 'yes' if not original_is_stop and mutated_is_stop else 'no',
                    'del_stopconda': 'yes' if original_is_stop and not mutated_is_stop else 'no'
                }
                mutation_results.append(result)

codon_df = pd.DataFrame(codon_library, columns=['codon'])
amino_acid_df = pd.DataFrame(list(codon_table.items()), columns=['codon', 'amino_acid'])
mutation_df = pd.DataFrame(mutation_results)
ori_mut_df = mutation_df[['ori_aa', 'mut_aa']]

codon_df.to_excel('codon_database.xlsx', index=False)
amino_acid_df.to_excel('amino_database.xlsx', index=False)
mutation_df.to_excel('mutation_results.xlsx', index=False)
ori_mut_df.to_excel('ori_mut.xlsx', index=False)

print("results have been saved to codon_database.xlsx, amino_database.xlsx, mutation_results.xlsx,ori_mut.xlsx") 