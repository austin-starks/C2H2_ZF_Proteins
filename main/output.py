"""
A file used to output information on a file.

This file stores functions that outputs information in a list of genes.

Author: Austin Starks
Data: June 25, 2019
"""

from datetime import date
import os.path
from texttable import Texttable  # pip install texttable
import matplotlib.pyplot as plt
import numpy
import database
import genetics
import helper
import count
import re


def info(alist, filename_txt="zf_protein_info.txt",
         filename_chart="zf_protein_chart.png", filename_tab="zf_protein_tabs.tsv",
         filename_tab_dis="zf_genes_w_disease.tsv", create_txt=True,
         create_tab=True, create_chart=True, create_tab_dis=True):
    """
    Creates a text file containing the information from a list of genes,
    a tab-delimited file of all of the info in genes, and a chart visualizing
    the number of proteins that have a certain number of domains.

    Parameter alist: a list of genes to create the file from.
    Precondtion: alist is a list of genes

    Parameter filename_txt: a string to name the info file ending in .txt
    Preconditon: filename_txt is a string ending in .txt

    Parameter filename_chart: a string to name the graph ending in .png
    Preconditon: filename_chart is a string ending in .png

    Parameter filename_tab: a string to name the tab delimited file ending in .txt
    Preconditon: filename_tab is a string ending in .txt

    Parameter create_txt: A boolean that says whether or not to create a new
    text file.
    Preconditon: create_txt is a bool

    Parameter create_tab: a boolean that says whether or not to create a new
    tab delimited file.
    Preconditon: create_tab is a bool

    Parameter create_chart: a boolean that says whether or not to create a new
    tab delimited file.
    Preconditon: create_chart is a bool
    """
    info_txt(alist, filename_txt, create_txt)
    info_tab(alist, filename_tab, create_tab)
    info_chart(alist, filename_chart, create_chart)
    info_tab_diseased(alist, filename_tab_dis, create_tab_dis)


def mutation_info(gene_list, ori_list, removal,
                  filename_mut="mutation_list_info.txt",
                  filename_DNA="DNA_binding_region_info.txt",
                  analysis='euclidean',
                  create_mut=True, create_DNA=True, create_spec=True):
    """
    Creates a text file containing the information from a list of genes with
    mutations, removed disease files, and a chart of the ED of the binding specifity

    Parameter gene_list: a list of genes to create the file from.
    Precondtion: alist is a list of genes

    Parameter ori_list: the original list of genes that has all mutations
    (including the removed mutations)
    Preconditon: ori_list is a superset of gene_list

    Parameter removal: a dictionary of diseases that were removed for gene_list
    Preconditon: removal is a dictionary

    Parameter filename_mut: a string to name the info file ending in .txt
    Preconditon: filename_mut is a string ending in .txt

    Parameter filename_DNA: a string to name the DNA-binding-spec info text file
    Preconditon: filename_DNA is a string ending in .txt

    Parameter filename_specifity: a string to name the graph ending in .png
    Preconditon: filename_specifity is a string ending in .png

    Parameter create_mut: A boolean that says whether or not to create a new
    text file about mutation info.
    Preconditon: create_mut is a bool

    Parameter create_DNA: a boolean that says whether or not to create a new
    file about the DNA-binding Specificity
    Preconditon: create_DNA is a bool

    Parameter create_spec: a boolean that says whether or not to create a new
    chart about the ED
    Preconditon: create_spec is a bool
    """
    removed_disease_list(ori_list, removal, create_dis_list=True)
    info_mutation_list(gene_list, filename_mut, create_mut)
    DNA_binding_list(gene_list, filename_DNA, create_DNA)
    # chart_binding_specifity(analysis, create_spec)


def info_txt(alist, filename_txt="zf_protein_info.txt", create_txt=True):
    """
    Creates a text file containing the information from a list of genes.

    Parameter alist: the list to create the file from.
    Precondtion: alist is a list

    Parameter full_dict: a dictionary of the full list of proteins. Used by
    statistics to gather information about the proteins
    Preconditon: full_dict is a dictionary of protein IDs

    Parameter filename_txt: a string to name the info file ending in .txt
    Preconditon: filename_txt is a string ending in .txt

    Parameter create_txt: A boolean that says whether or not to create a new
    text file.
    Preconditon: create_txt is a bool
    """
    filename = helper.edit_filename(
        filename_txt, "sibling", foldername='output')
    if not os.path.exists(filename) or create_txt:
        full_dict = genetics.Protein.full_protein_dict()
        x = Texttable()
        today = date.today()
        month = str(today.month)
        day = str(today.day)
        if len(month) < 2:
            month = "0" + month
        if len(day) < 2:
            day = "0" + day
        with open(filename, "w+", encoding='utf-8') as f:
            stats = count.zf_protein_stats(alist, full_dict)
            f.write("This text file was generated by the Python module " +
                    os.path.basename(__file__) + "." + "\n\n" + "This file contains information " +
                    "about genes, including their \ngene IDs, their common names, and the list of " +
                    "proteins they \ncode for. This file also contains information about the \n" +
                    "proteins including the number of ZF domains in each protein, \ntheir core " +
                    "sequences, and statistics about the proteins." +
                    "\n\nAuthor: Austin Starks\nDate: " + month + "/" + day + "/" + str(today.year))
            f.write(stats + "Gene info: \n\nNote: A core sequence of '????' " +
                    "indicates that domain is likely invalid.\n\n")
            categories = ["Gene name", "Protein ID",
                          "Core Seqs", "# ZFs", "Diseases"]
            x.add_row(categories)
            for gene in alist:
                protein_list = gene.protein_list()
                for protein in protein_list:
                    row = [gene.get_gene_name(), protein.get_protein_id(), protein.core_sequences(),
                           protein.num_domains(), gene.disease_list_string()]
                    x.add_row(row)
            f.write(x.draw())
        print(filename_txt + " created successfully!")


def info_tab(alist, filename_tab="zf_protein_tabs.tsv", create_tab=True):
    """
    Creates a tab delimited file containing the information from a list of genes.

    Parameter alist: the list to create the file from.
    Precondtion: alist is a list

    Parameter filename_tab: a string to name the tab delimmited file ending in .txt
    Preconditon: filename_tab is a string ending in .txt

    Parameter create_tab: a boolean that says whether or not to create a new
    tab delimited file.
    Preconditon: create_tab is a bool
    """
    filename = helper.edit_filename(
        filename_tab, "sibling", foldername='output')
    if not os.path.exists(filename) or create_tab:
        with open(filename, "w+", encoding='utf-8') as g:
            head = ['Gene Name', 'Gene ID', 'Protein ID', 'Core',
                    'Num Domains', 'Disease Info']
            g.write('#' + '\t'.join(str(x) for x in head) + '\n')
            for gene in alist:
                protein_list = gene.protein_list()
                for protein in protein_list:
                    row = [gene.get_gene_name(), gene.get_gene_id(), protein.get_protein_id(),
                           protein.core_sequences(
                               False), protein.num_domains(),
                           gene.disease_str_without_full_info()]
                    g.write('\t'.join(str(x) for x in row) + '\n')
        print(filename_tab + " created successfully!")


def info_chart(alist, filename_chart="zf_protein_chart.png", create_chart=True):
    """
    Creates a chart visualizing the number of proteins that have a certain
    number of domains.

    Parameter alist: a list of genes
    Precondtion: alist is a list of genes

    Parameter filename_tab: a string to name the chart ending in png
    Preconditon: filename_tab is a string ending in .png

    Parameter create_chart: a boolean that says whether or not to create a chart
    Preconditon: create_chart is a bool
    """
    filename = helper.edit_filename(
        filename_chart, "sibling", foldername='output')
    if not os.path.exists(filename) or create_chart:
        dict_domains = count.dict_domains(alist)
        left = []
        height = []
        for key in dict_domains:
            left.append(key)
            height.append(dict_domains[key])
        tick_label = [str(i) for i in left]
        plt.bar(left, height, tick_label=tick_label, width=0.4, color=["blue"])
        plt.xlabel("Number of proteins")
        plt.ylabel("Number of ZFs")
        plt.title("Number of ZFs in ZF proteins")
        figure = plt.gcf()
        figure.set_size_inches(18, 9)
        plt.rc('font', size=8)
    #    plt.show() #Use this to see the chart.
        plt.savefig(filename, format="png", dpi=100)
        print(filename_chart + " created successfully!")


def info_tab_diseased(alist, filename_tab_dis="zf_genes_w_disease.tsv",
                      create_tab_dis=True):
    """
    Creates a tab delimited file containing the disease information from a list
    of genes with disease information

    Parameter alist: the list to create the file from.
    Precondtion: alist is a list

    Parameter filename_tab:_dis a string to name the tab delimmited file ending in .txt
    Preconditon: filename_tab_dis is a string ending in .txt

    Parameter create_tab_dis: a boolean that says whether or not to create a new
    tab delimited file.
    Preconditon: create_tab_dis is a bool
    """
    filename = helper.edit_filename(
        filename_tab_dis, "sibling", foldername='output')
    if not os.path.exists(filename) or create_tab_dis:
        with open(filename, "w+", encoding='utf-8') as g:
            head = ['Gene Name', 'Gene ID', 'Chromosome', 'Position of SNV',
                    'Allele Change', 'Amino Change', 'Disease Desc', 'SNP', 'Disease Source']
            g.write('#' + '\t'.join(str(x) for x in head) + '\n')
            for gene in alist:
                disease_list = gene.disease_list()
                for disease in disease_list:
                    if disease.has_full_info():
                        row = [gene.get_gene_name(), gene.get_gene_id(),
                               disease.get_chromosome(), disease.get_string_full_position(),
                               disease.get_allele_change(), disease.get_amino_change(),
                               disease.get_disease_name(), disease.get_dbSNP(),
                               disease.get_source()]
                        g.write('\t'.join(str(x) for x in row) + '\n')
        print(filename_tab_dis + " created successfully!")


def removed_disease_list(gene_list, removal, filename="removed_disease_info.txt",
                         create_dis_list=False):
    """
    Outputs info about removed diseases.

    This procedure outputs information about diseases that are to be removed.
    Information includes the protein

    Parameter gene_list: a list of genes.
    Preconditon: gene_list is a list of Genes.

    Parameter removal: a set of diseases to remove
    Preconditon removal is a set of Diseases

    Parameter filename: a string to name the info file ending in .txt
    Preconditon: filename_txt is a string ending in .txt

    Parameter create_dis_list: A boolean that says whether or not to create a new
    text file.
    Preconditon: create_dis_list is a bool
    """
    assert type(gene_list) == list, "Gene list must be a list."
    assert type(removal) == dict, "Removal must be a dict of diseases."
    assert type(filename) == str, "Filename must be a string."
    file = helper.edit_filename(filename, "sibling", "output")
    if not os.path.exists(file) or create_dis_list:
        today = date.today()
        month = str(today.month)
        day = str(today.day)
        abool = False
        if len(month) < 2:
            month = "0" + month
        if len(day) < 2:
            day = "0" + day
        with open(file, "w+", encoding='utf-8') as f:
            f.write("This file was generated by the Python module " +
                    os.path.basename(__file__) + '\n\n' +
                    "This file contains information about the diseases that\n" +
                    "were removed from the disease list. Specifically, this\n" +
                    "file outputs the gene associated with the disease,\n" +
                    "the list of proteins that the gene codes for, the allelic\n" +
                    "and amino acid substitution, the position of the mutation,\n" +
                    "and other important info.\n\n" +
                    "Author: Austin Starks\nDate: " + month + "/" + day + "/" +
                    str(today.year) + '\n\n')
            f.write("Number of diseases removed: " + str(len(removal)))
            f.write(
                '\n----------------------------------------------\n\n')
            for gene in gene_list:
                diseases = gene.disease_list()
                for disease in diseases:
                    if disease in removal:
                        abool = True
                        f.write('\nGene name: ' + gene.get_gene_name() +
                                '\nGene ID:' + gene.get_gene_id() + '\n')
                        f.write('Disease description: ' +
                                disease.get_disease_name() + '\n')
                        f.write('Disease source: ' +
                                disease.get_source() + '\n')
                        f.write('Disease allelic position: ' +
                                disease.get_string_full_position() + '\n')
                        f.write('Disease allelic change: ' +
                                disease.get_allele_change() + '\n')
                        f.write('Disease amino change: ' +
                                disease.get_amino_change() + '\n')
                        f.write('Notes: ' + removal[disease] + '\n')
                if abool:
                    abool = False
                    proteins = gene.protein_list()
                    for protein in proteins:
                        f.write('\nProtein ID: ' +
                                protein.get_protein_id() + '\n\n')
                        f.write('Protein Seq:\n' +
                                str(protein.prosequence()) + '\n')
                    f.write(
                        '\n----------------------------------------------\n')
    if create_dis_list:
        print("Text file about removed diseases created successfully.")
    else:
        print("Text file about removed diseases is already there. \n  Change" +
              " create_dis_list to be True.")


def info_mutation_list(mut_list, filename_mut="mutation_list_info.txt",
                       create_mut=True):
    """
    Parameter alist: a list of genes.
    Preconditon: alist is a list of Genes.

    Parameter filename_mut: a string to name the info file ending in .txt
    Preconditon: filename_mut is a string ending in .txt

    Parameter create_mut: A boolean that says whether or not to create a new
    text file.
    Preconditon: create_mut is a bool
    """
    assert type(mut_list) == list, "Gene list must be a list."
    assert type(filename_mut) == str, "Filename must be a string."
    assert type(create_mut) == bool, "Create mut must be a string."
    print("Gathering statistics about mutations.")
    file = helper.edit_filename(filename_mut, "sibling", "output")
    dna_count, ion_count, zf_muts = count.mutation_gene_stats(mut_list)
    diseased_list = genetics.Gene.diseased_gene_list(print_=False)
    initial = count.total_mutations(diseased_list)
    total = count.total_mutations(mut_list)
    stats_pro_dis = count.protein_and_disease_stats(mut_list)
    missense = count.missense_mutation_stats(mut_list)
    if not os.path.exists(file) or create_mut:
        today = date.today()
        month = str(today.month)
        day = str(today.day)
        if len(month) < 2:
            month = "0" + month
        if len(day) < 2:
            day = "0" + day
        with open(file, "w+", encoding='utf-8') as f:
            f.write("This file was generated by the Python module " +
                    os.path.basename(__file__) + '\n\n' +
                    "This file contains information about the effects of\n" +
                    "mutations on human zinc finger proteins. Specifically,\n" +
                    "it lists the original and mutated zinc finger domains, and \n" +
                    "describes if that mutation is on a zinc/dna-binding region.\n" +
                    "This file also gives statistics about the region.\n\n"
                    "Author: Austin Starks\nDate: " + month + "/" + day + "/" +
                    str(today.year) + '\n\n')
            f.write(
                '\n----------------------------------------------\n\n')
            f.write('Mutation Stats:\n')
            f.write('The following are statistics about mutations in ZF domains.\n')
            f.write('The diseases were "filtered" to remove mutations that were\n')
            f.write('faulty.\n\nFor example:\nHis184Pro\n')
            f.write("But when you check the 184th protein, it isn't Histidine.\n\n")
            # f.write("Initial number of mutations           :  " +
            #        str(initial) + '\n')
            f.write("# of mutations after filtering        :  " +
                    str(total) + '\n')
            f.write("# of mutations in DNA-binding region  :  " +
                    str(dna_count) + '\n')
            f.write("# of mutations in ion-binding region  :  " +
                    str(ion_count) + '\n')
            f.write("# of mutations in zinc-finger domains :  " +
                    str(zf_muts) + '\n')
            f.write("# of missense mutations               :  %s\n" %
                    missense['Total missense mutations'])
            f.write("# of nonsense mutations               :  %s\n" %
                    (missense['# of nonsense mutations']))
            f.write("# of nonsense mut in section of 1 domain   :  %s\n" %
                    (missense['# of nonsense mutations eliminating a part of 1 domain']))
            f.write("# of nonsense mut eliminating entire domain:  %s\n" %
                    (missense['# of nonsense mutations eliminating an entire domain']))
            f.write(
                '\n----------------------------------------------\n')
            f.write('\nProtein Stats:\n')
            f.write('The following are statistics about the amino acids. Note the\n')
            f.write('apparent overlap between these stats and the mutation stats.\n')
            f.write("However, these stats are important for performing a Fisher's\n")
            f.write("Exact test.\n\n")
            for i in range(len(stats_pro_dis)):
                tmp = stats_pro_dis[i]
                for j in tmp:
                    msg = j
                    while len(msg) < 52:
                        msg += ' '
                    f.write(msg + ': ' + str(tmp[j]) + '\n')
            f.write(
                '\n----------------------------------------------\n')
            f.write("\nMissense Stats:\n")
            f.write("The following are statistics about missense mutations. These\n")
            f.write("stats are used to get information on the number of missense\n")
            f.write(
                "mutations (point mutations that don't end with termination of\n")
            f.write("the protein transcript).\n\n")
            for ele in missense:
                msg = ele
                while len(msg) < 40:
                    msg += ' '
                f.write(msg + ": " + str(missense[ele]) + '\n')
            f.write(
                '\n----------------------------------------------\n')
            for gene in mut_list:
                proteins = gene.protein_list()
                f.write("\nGene Name: " + gene.get_gene_name() + '\n')
                f.write("Gene ID: " + gene.get_gene_id() + '\n')
                f.write("Chromosome: %s" % gene.disease_list()
                        [0].get_chromosome() + '\n')
                for protein in proteins:
                    f.write("Protein ID: " + protein.get_protein_id() + '\n')
                    domains = protein.domain_list()
                    for domain in sorted(domains):
                        mutation_dict = domain.get_mutation_dictionary()
                        ori_sequence = domain.domsequence()
                        f.write("\n---\nDomain number: " +
                                str(domain.get_domain_number()) + ' | ')
                        f.write(
                            "From " + str(domain.get_start_position()) + " to ")
                        f.write(str(domain.get_end_position()) + '\n')
                        f.write("Score            : " +
                                str(domain.get_score()) + '\n')
                        if domain.is_valid():
                            f.write("Original Sequence: " +
                                    ori_sequence + '\n\n')
                            for mut_seq in mutation_dict:
                                notes = mutation_dict[mut_seq][0]
                                disease = mutation_dict[mut_seq][4]
                                if notes != "none":
                                    notes = notes + " region"
                                f.write("Mutated Sequence : " + mut_seq + '\n')
                                f.write("Original Amino   : " +
                                        mutation_dict[mut_seq][1] + '\n')
                                f.write("Mutated Amino    : " +
                                        mutation_dict[mut_seq][2] + '\n')
                                f.write("Mutation Index   : " +
                                        str(mutation_dict[mut_seq][3]) + '\n')
                                f.write("Amino Acid Change: " +
                                        mutation_dict[mut_seq][4].get_amino_change() + '\n')
                                f.write("Allelic Change   : " +
                                        mutation_dict[mut_seq][4].get_allele_change() + '\n')
                                f.write("Allele Position  : %s\n" %
                                        mutation_dict[mut_seq][4].get_string_full_position())
                                f.write("Disease Name     : " +
                                        disease.get_disease_name() + '\n')
                                f.write("Source           : " +
                                        disease.get_source() + '\n')
                                f.write("Notes            : " + notes + '\n')
                                f.write('\n')
                        else:
                            f.write("Invalid C2H2 ZF Domain\n\n")
                f.write(
                    '\n----------------------------------------------\n\n')
            f.write(
                '\n\n----------------------------------------------\n\n')
            f.write("The following diseases did not affect a zinc finger domain:\n")
            for gene in mut_list:
                diseases = gene.disease_list()
                for disease in diseases:
                    if not disease.get_hit_domain():
                        f.write("Gene name: " + gene.get_gene_name() + ' | ')
                        f.write("Gene ID: " + gene.get_gene_id() + '\n')
                        f.write("Disease: " + disease.get_disease_name() + '\n')
                        f.write("Allelic position: " +
                                disease.get_string_full_position() + ' | ')
                        f.write("Allelic substitution: " +
                                disease.get_allele_change() + ' | ')
                        f.write("Chromosome: " +
                                str(disease.get_chromosome()) + '\n')
                        f.write("Amino acid substitution: " +
                                disease.get_amino_change() + ' | ')
                        f.write("Source: " + disease.get_source() + '\n')
                        f.write('\n')
                f.write(
                    '\n----------------------------------------------\n')
    print(filename_mut + " was created successfully!")


def DNA_binding_list(DNA_list, filename_DNA="DNA_binding_region_info.txt",
                     create_DNA=True):
    """
    Parameter DNA_list: a list of genes.
    Preconditon: DNA_list is a list of Genes.

    Parameter filename_DNA: a string to name the info file ending in .txt
    Preconditon: filename_DNA is a string ending in .txt

    Parameter create_DNA: A boolean that says whether or not to create a new
    text file.
    Preconditon: create_DNA is a bool
    """
    assert type(DNA_list) == list, "Gene list must be a list."
    assert type(filename_DNA) == str, "Filename must be a string."
    assert type(create_DNA) == bool, "Create mut must be a string."
    print("Outputting DNA-binding region info")
    file = helper.edit_filename(filename_DNA, "sibling", "output")
    if not os.path.exists(file) or create_DNA:
        today = date.today()
        month = str(today.month)
        day = str(today.day)
        if len(month) < 2:
            month = "0" + month
        if len(day) < 2:
            day = "0" + day
        with open(file, 'w+', encoding='utf-8') as f:
            f.write("This file was generated by the Python module " +
                    os.path.basename(__file__) + '\n\n')
            f.write(
                "This file contains information about missense mutations in the\n")
            f.write("DNA-binding region. Specifically, this file gives information\n")
            f.write("on the gene name, the protein ID, the domain number, the\n")
            f.write("HMMER score, the disease information, as well as other\n")
            f.write(
                "relevant information about mutations in the DNA-binding region.\n\n")
            f.write("Author: Austin Starks\nDate: " + month + "/" + day + "/" +
                    str(today.year) + '\n\n')
            for gene in DNA_list:
                proteins = gene.protein_list()
                for protein in proteins:
                    domains = protein.domain_list()
                    for domain in sorted(domains):
                        mutation_dict = domain.get_mutation_dictionary()
                        ori_sequence = domain.domsequence()
                        if domain.is_valid:
                            for mut_seq in mutation_dict:
                                notes = mutation_dict[mut_seq][0]
                                disease = mutation_dict[mut_seq][4]
                                amino_change = mutation_dict[mut_seq][4].get_amino_change(
                                )
                                allele_change = mutation_dict[mut_seq][4].get_allele_change(
                                )
                                position = mutation_dict[mut_seq][4].get_string_full_position(
                                )
                                Ter = amino_change[-3:]
                                if 'dna' in notes and Ter != 'Ter':
                                    f.write("\nGene Name: " +
                                            gene.get_gene_name() + '\n')
                                    f.write("Protein ID: %s\n" %
                                            protein.get_protein_id())
                                    f.write("Domain number: " +
                                            str(domain.get_domain_number()) + ' | ')
                                    f.write(
                                        "From " + str(domain.get_start_position()) + " to ")
                                    f.write(
                                        str(domain.get_end_position()) + '\n')
                                    f.write("Score            : " +
                                            str(domain.get_score()) + '\n')
                                    f.write("Original Sequence: " +
                                            ori_sequence + '\n')
                                    f.write("Mutated Sequence : " +
                                            mut_seq + '\n')
                                    f.write("Original Amino   : " +
                                            mutation_dict[mut_seq][1] + '\n')
                                    f.write("Mutated Amino    : " +
                                            mutation_dict[mut_seq][2] + '\n')
                                    f.write("Mutation Index   : " +
                                            str(mutation_dict[mut_seq][3]) + '\n')
                                    f.write("Amino Acid Change: " +
                                            amino_change + '\n')
                                    f.write("Allelic Change   : " +
                                            allele_change + '\n')
                                    f.write("Allele Position  : " +
                                            position + '\n')
                                    f.write("Disease Name     : " +
                                            disease.get_disease_name() + '\n')
                                    f.write("Source           : " +
                                            disease.get_source() + '\n')
                                    f.write("Sequence: \n%s\n\n" %
                                            protein.prosequence())
                                    f.write(
                                        '\n-----------------------------------------------\n')
    print(filename_DNA + " was created successfully!")


def chart_binding_specifity(analysis='euclidean', create_spec=True):
    """
    Creates a histogram displaying the change in binding specifity caused by
    missense mutations in the DNA-binding region.

    This function scans a list of PWM (obtained manually from zf.princeton.edu)
    and outputs a charge that displays the euclidean distance of each mutation in
    the DNA-binding region.

    Parameter analysis: the method of analyzing
    Preconditon: analysis is pearsons or Euclidean

    Parameter create_spec:  a boolean that says whether or not to create a new
    chart
    Preconditon: create_spec is a bool
    """
    assert 'pearson' in analysis or 'euclidean' in analysis
    filename_specifity = analysis + '.png'
    filename = helper.edit_filename(
        filename_specifity, "sibling", foldername='output')
    if not os.path.exists(filename) or create_spec:
        # with open(filename, 'w+') as f:
        #     for ele in adict:
        #         f.write('%s\t%s\n' % (ele, adict[ele]))
        if 'euclidean' in analysis:
            xlabel = numpy.arange(0.0, 2.6, 0.2)
        else:
            xlabel = numpy.arange(0.0, 1.0, 0.05)
        arr = count.quantify_dna_specifity(
            foldername='analysis', analysis=analysis)
        plt.hist(arr, color='lightgrey', ec='black', bins=xlabel)
        plt.xticks(xlabel, fontsize=18)
        plt.yticks(range(0, 10), fontsize=18)
        plt.ylabel("Frequency", fontsize=24)
        if 'euclidean' in analysis:
            plt.xlabel("Euclidean Distance", fontsize=24)
        else:
            plt.xlabel("Pearson Correlation Coefficient", fontsize=24)
        plt.title("Predicted Change in DNA Binding Specificity", fontsize=36)
        figure = plt.gcf()
        figure.set_size_inches(18, 9)
    #    plt.show()  # Use this to see the chart.
        plt.savefig(filename, format="png", dpi=100)
        print(filename_specifity + " created successfully!")


def population_info(pop_list, full_dict_proteins=dict(),
                    filename_pop="population_info.tsv", create_pop=True):
    """
    Outputs population info.

    This function outputs a tab file for population information. This file can
    then be searched through ExAC to get population information.

    Parameter pop_list: a list of Genes
    Preconditon: pop_list is a list of Genes

    Parameter filename_pop: the name of the file
    Preconditon: filename_pop is a string

    Parameter create_pop: says whether to create a new file
    Preconditon: create_pop is a boolean
    """
    print("Creating population info file.")
    if full_dict_proteins == dict():
        full_dict_proteins = genetics.Protein.full_protein_dict(
            "Homo_sapiens.GRCh37.pep.all.fa")
    filename = helper.edit_filename(
        filename_pop, "sibling", foldername='output')
    if not os.path.exists(filename) or create_pop:
        with open(filename, "w+", encoding='utf-8') as g:
            head = ['Gene ID', 'Protein ID', 'Domain Number',
                    'Chromosome-Start_Pos-End_Pos', 'Core Indices',
                    'Ion-Binding Indices', 'Hydrophobic Indices']
            g.write('#' + '\t'.join(str(x) for x in head) + '\n')
            for gene in pop_list:
                proteins = gene.protein_list()
                for protein in proteins:
                    description = full_dict_proteins[protein.get_protein_id()]
                    match = re.search(
                        r"GRCh37:(\w+):(\w+):(\w+)", str(description))
                    try:
                        chromosome = match.group(1)
                        if len(chromosome) > 2:
                            int(X)
                        start_pos = match.group(2)
                        end_pos = match.group(3)
                        domains = protein.domain_list()
                        for domain in domains:
                            key_pos = domain.get_key_positions()
                            if len(key_pos) == 0:
                                pass
                            else:
                                g.write(gene.get_gene_id() + '\t')
                                g.write(protein.get_protein_id() + '\t')
                                g.write(str(domain.get_domain_number()) + '\t')
                                g.write(str(chromosome) + '-' +
                                        str(start_pos) + '-' + str(end_pos) + '\t')
                                g.write(str(key_pos['Core Indices']) + '\t')
                                g.write(
                                    str(key_pos['Ion-Binding Indices']) + '\t')
                                g.write(str(key_pos['Hydrophobic Indices']))
                                g.write('\n')
                    except:
                        pass
    print(filename_pop + " created successfully!")


# chart_binding_specifity()
