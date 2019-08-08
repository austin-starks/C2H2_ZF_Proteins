"""
A program to extract, store, and transform biological data.

This program is used to extract protein IDs from HMMER output.From the output,
it extracts the gene IDs and gene names for each gene. Finally, it searches
through a variety of gene-disease assoociation databases and determines if a
gene is associated with pathological disorder.

Afterwards, this program searches through databases such as ClinVar and UniProt
that have information on the specific allelic locations and amino acid
substition. Finally it verifies those substitions by looking at a chromosome
dictionary, and removes diseases that have impossible substitions.

Finally, this program filters out the proteins and removes all except the
longest trandscript. Then it looks at mutations in that domain

This program was built to observe Zinc Finger Proteins. However, the program was
built to be able to be used with any HMMER output and any FASTA file. Additionally,
there are many useful standalone functions that can be used to study other
biological questions.

Requirements: Python 3 and HMMER 2.3.2

Author: Austin Starks
Data: June 14, 2019
"""

import genetics
import output
import database
import helper
import count
#### Temp ####


def execute_program(HMMER=False, filter=True, score_threshold=False):
    """
    Executes the program and creates a gene list.

    This function is used if this is the first time executing the progam.
    This fuction creates a gene list and stores all of the relevant information
    from a variety of databases and HMMER output.

    Parameter HMMER: says whether to create a new HMMER file.
    Preconditon: HMMER is a bool

    Parameter filter: says whether to filter out other proteins except the longest
    transcript
    Preconditon: filter is a bool

    Parameter score_threshold: says whether to get only output proteins with a
    score above 17.7
    Preconditon: score_threshold is a bool
    """
    data = genetics.hmmer_output(create_new=HMMER)
    full_dict_proteins = genetics.Protein.full_protein_dict()
    protein_list = genetics.Protein.create_protein_set(data)
    gene_list = genetics.Gene.create_gene_list(
        protein_list, full_dict_proteins)
    genetics.Protein.get_zf_protein_info(data, gene_list)
    genetics.Protein.gene_list_set_sequence(gene_list, full_dict_proteins)
    database.search(gene_list)
    database.search_mutation_info(gene_list)
    genetics.Protein.remove_proteins(gene_list, score_threshold)
    output.info(gene_list, create_txt=True, create_tab=True,
                create_chart=True, create_tab_dis=True)
    ######################################################################


def create_gene_list(get_domain_info=False):
    """
    Creates a gene list.

    This method creates a list of genes with all of the needed information in
    them. This includes protein information, disease information, and domain
    information.

    The gene list is created from existing tab files.

    Parameter get_domain_info: boolean that says whether to get domain info.
    Preconditon: get_domain_info is a boolean

    Returns: a list of genes
    """
    assert type(get_domain_info) == bool, 'zf_info must be a boolean'
    try:
        data = genetics.hmmer_output()
        gene_list = genetics.Gene.create_gene_list_tab(
            zf_domain_info=get_domain_info)
        protein_dict = genetics.Protein.full_protein_dict()
        genetics.Protein.get_proteins(gene_list, data)
        genetics.Protein.gene_list_set_sequence(gene_list, protein_dict)
        genetics.Domain.gene_list_set_domain_positions_sequence(gene_list)
    except helper.DoesNotExistError:
        print('Error, tab files do not exist. Creating new list of genes')
        gene_list = execute_program()
    return gene_list


def unconfident_domains():
    """
    Looks at the removed ZF proteins.

    This function looks at the list of proteins that have a low score theshold
    and outputs statistics about them.

    Parameter remove_all: says whether to remove proteins with a low score
    threshold even if it is a known ZF protein.
    Preconditon: remove_all is a bool
    """
    gene_list = genetics.Gene.diseased_gene_list()
    gene_list = genetics.Protein.remove_protein_info(
        gene_list, True)
    removal = genetics.verify_protein_sequence(gene_list)
    genetics.remove_bad_diseases(gene_list, removal)
    genetics.gene_list_process_mutations(gene_list, False)
    output.info_mutation_list(
        gene_list, filename_mut="non_zf_info.txt", create_mut=True)
    adict = count.protein_and_disease_stats(gene_list)
    for i in range(len(adict)):
        tmp = adict[i]
        for j in tmp:
            print(j + ': ' + str(tmp[j]))
    return gene_list


def mutated_gene_list(score_threshold=True):
    """
    This function returns a list of genes with the mutations in the correct
    position.

    Parameter score_threshold: says whether to get only output proteins with a
    score above 17.7
    Preconditon: score_threshold is a bool

    Parameter remove_all: says whether to remove proteins with a low score
    threshold even if it is a known ZF protein.
    Preconditon: remove_all is a bool
    """
    ori_list = genetics.Gene.diseased_gene_list()
    genetics.Protein.remove_proteins(ori_list, score_threshold)
    removal = genetics.verify_protein_sequence(ori_list)
    gene_list = genetics.remove_bad_diseases(ori_list, removal)
    genetics.gene_list_process_mutations(gene_list)
    if score_threshold:
        genetics.Disease.create_GRCh38_file(gene_list)
        genetics.Disease.convert_to_GRCh37_from_file(gene_list)
    output.mutation_info(gene_list, ori_list, removal)
    output.chart_binding_specifity('pearsons')
    return gene_list


def population_variation_info():
    """
    This function returns a list of Genes that has detailed information
    for each protein.
    """
    data = genetics.hmmer_output(
        create_new=False, fasta="Homo_sapiens.GRCh37.pep.all.fa", filename='hmmer37.txt')
    full_dict_proteins = genetics.Protein.full_protein_dict(
        "Homo_sapiens.GRCh37.pep.all.fa")
    protein_list = genetics.Protein.create_protein_set(data)
    gene_list = genetics.Gene.create_gene_list(
        protein_list, full_dict_proteins)
    genetics.Protein.get_zf_protein_info(data, gene_list)
    genetics.Protein.gene_list_set_sequence(gene_list, full_dict_proteins)
    genetics.Domain.gene_list_set_domain_positions_sequence(gene_list)
    genetics.Domain.store_key_positions(gene_list)
    output.population_info(gene_list, full_dict_proteins)


if __name__ == '__main__':
    population_variation_info()
