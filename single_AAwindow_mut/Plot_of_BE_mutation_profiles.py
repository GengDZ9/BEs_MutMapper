import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
import numpy as np

hydrophobic_amino_acids = ['Ala', 'Val', 'Leu', 'Ile', 'Pro', 'Phe', 'Trp', 'Met']
hydrophilic_amino_acids = ['Gly', 'Ser', 'Thr', 'Cys', 'Tyr', 'Asn', 'Gln']
positive_amino_acids = ['Lys', 'Arg', 'His']
negative_amino_acids = ['Asp', 'Glu']

amino_acids = hydrophobic_amino_acids + hydrophilic_amino_acids + positive_amino_acids + negative_amino_acids

try:
    ori_mut_df = pd.read_excel('ori_mut.xlsx')
except FileNotFoundError:
    print("error：no 'ori_mut.xlsx'")
else:
    mutation_counts = pd.DataFrame(0, index=amino_acids, columns=amino_acids)
    for _, row in ori_mut_df.iterrows():
        original_amino_acid = row['ori_aa']
        mutated_amino_acid = row['mut_aa']
        if original_amino_acid in amino_acids and mutated_amino_acid in amino_acids:
            mutation_counts.loc[original_amino_acid, mutated_amino_acid] += 1

    plt.rcParams['figure.dpi'] = 300

    plt.figure(figsize=(12, 10))

    ax = sns.heatmap(mutation_counts, annot=True, fmt='d', cmap='YlGnBu',
                xticklabels=amino_acids, yticklabels=amino_acids)
    cbar = ax.collections[0].colorbar        

    plt.title('Amino Acid Mution Heatmap',fontweight='bold')
    plt.xlabel('Mution',fontweight='bold')
    plt.ylabel('Wild Type',fontweight='bold')
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    cbar.ax.tick_params(labelsize=12)


    plt.savefig('amino_acid_mutation_heatmap.pdf')
    
    print("plot has been saved in amino_acid_mutation_heatmap.pdf")


mutation_df = pd.read_excel('mutation_results.xlsx')
synonymous_mutation_df = mutation_df[(mutation_df['mut_aa'] != 'Stop') & (mutation_df['aa_change'] == 'yes')]
nonsynonymous_mutation_df = mutation_df[(mutation_df['mut_aa'] != 'Stop') & (mutation_df['aa_change'] == 'no')]



synonymous_pivot_table = pd.pivot_table(synonymous_mutation_df, values='ori_codon', index='mut_aa', columns='mutsite_num', aggfunc='count', fill_value=0)
synonymous_pivot_table = synonymous_pivot_table.reindex(amino_acids)
synonymous_pivot_table = synonymous_pivot_table.fillna(0)


nonsynonymous_pivot_table = pd.pivot_table(nonsynonymous_mutation_df, values='ori_codon', index='mut_aa', columns='mutsite_num', aggfunc='count', fill_value=0)
nonsynonymous_pivot_table = nonsynonymous_pivot_table.reindex(amino_acids)
nonsynonymous_pivot_table = nonsynonymous_pivot_table.fillna(0)


combined_pivot_table = pd.concat([synonymous_pivot_table.add_suffix('_sense'), nonsynonymous_pivot_table.add_suffix('_nosense')], axis=1)


plt.figure(figsize=(12, 8))
max_value = combined_pivot_table.max().max()
step = 5
if max_value > 20:
    step = 2 if max_value <= 50 else 5
tick_values = np.arange(0, int(max_value) + step, step=step)
cbar_kws = {
    "ticks": tick_values,
    "format": lambda x, pos: f'{int(x)}'
}


ax = sns.heatmap(combined_pivot_table, annot=True, fmt='.0f', cmap='YlGnBu', cbar_kws=cbar_kws)
cbar = ax.collections[0].colorbar


plt.xlabel('Number of mutation sites and Mutation Type', fontsize=12)
plt.ylabel('Amino acid', fontsize=12)
plt.title('Amino acid mutations at different number of base mutation sites (Synonymous and Nonsynonymous Mutations)', fontsize=16)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

cbar.ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig('Cluster_amino_acid_mutation_heatmap.pdf')

print("plot has been saved in Cluster_amino_acid_mutation_heatmap.pdf。")
