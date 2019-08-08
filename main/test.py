from project import *
from genetics import *
import re
from Bio.Seq import Seq

full_dict_proteins = Protein.full_protein_dict()
gene_list = mutated_gene_list(score_threshold=True)
for gene in gene_list:
    if gene.get_gene_name() == 'TRPS1':
        break
protein = gene.protein_list()[0]
description = full_dict_proteins[protein.get_protein_id()].description
position = re.search(r'GRCh38:8:(\d+):(\d+)', description)
start_pos = int(position.group(1))
end_pos = int(position.group(2))
# Use full dict proteins to get the exact allelic position of the protein.
adict = chromo_dict('test_chromosomes')
mut_pos_arr = [2560, 2686, 976, 164, 2627]
reference = str(Seq(
    'gattttaattctctcctaccgggcgcagcggaggagctggaggtggttggacccgagggg'.upper()).reverse_complement())
index = adict[8].find(reference)
print('Index    ', index)
print('End pos  ', end_pos)
print('Differenc', end_pos - adict[8].find(reference))
for mut_pos in mut_pos_arr:
    print(adict[8][index + mut_pos - 3:index + mut_pos + 2])
