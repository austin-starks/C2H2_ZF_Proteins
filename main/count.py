"""
A file used to store counting functions.

This file stores functions static functions that have the job of counting
something. This file (obviously) doesn't store functions that count things
within an object; it is up to the object to do that.

Aithor: Austin Starks
Date: July 8, 2019
"""
import collections
import genetics
import helper
import re
import os.path
import glob
import statistics
import numpy
import scipy.stats
import scipy.spatial


def dict_domains(alist):
    """
    This function counts how many proteins have x domains from a list of
    genes.

    For example, this function will get how many proteins have 1 domain, how
    many have 2 domains, etc.

    Parameter alist: a list of genes.
    Precondition: alist is a list of genes.

    Returns: a dictionary mapping number of domains to count.
    """
    assert type(alist) == list, "alist must be a list"
    assert type(alist[0]) == genetics.Gene, "elements in alist must be Genes"
    assert type(alist[-1]) == genetics.Gene, "elements in alist must be Genes"
    adict = dict()
    for gene in alist:
        pro_list = gene.protein_list()
        for protein in pro_list:
            num_domains = protein.num_domains()
            if num_domains not in adict:
                adict[num_domains] = 1
            else:
                tmp = adict[num_domains]
                adict[num_domains] = tmp + 1
    return dict(sorted(adict.items()))


def zf_protein_stats(alist, pro_dict, fasta="Homo_sapiens.GRCh38.pep.all.fa"):
    """
    This function takes the data in a list of genes and extracts important
    statistics about the list. It also outputs a chart corresponding the
    number of proteins that have a corresponding number of domains.

    Parameter alist: a list of genes to get statistics from
    Preconditon: alist is a list of genes

    Parameter pro_dict: a full list of proteins that information was obtained from.
    Preconditon: pro_dict is a dictionary

    Parameter filename_chart: a string to name the graph ending in .png
    Preconditon: filename_chart is a string ending in .png

    Returns: a string of statistics
    """
    msg = "\n\n-------------------------------------\n\nStatistics: \n\n"
    # Get dictionay of domains and display it.
    dict_domain = dict_domains(alist)
    numproteins = 0
    for gene in alist:
        numproteins += gene.numproteins()
    numdomains = total_num_domains(alist)
    total_pro_processed = len(pro_dict)
    total_genes_processed = len(genetics.Gene.full_gene_id_set(fasta))
    numgenes = len(alist)
    num_diseased_zf_genes = count_genes_w_disease(alist)
    num_diseased_zf_genes_w_full_dis = genes_w_full_dis_info(alist)
    percent_tot_zf_genes = str(
        round(numgenes / total_genes_processed * 100, 1))
    percent_zf_genes_dis = str(
        round(num_diseased_zf_genes / numgenes * 100, 1))
    percent_zf_g_full_dis = str(
        round(num_diseased_zf_genes_w_full_dis / numgenes * 100, 1))
    msg = msg + "Total genes processed: " + str(total_genes_processed) + "\n" +\
                "Total ZF genes found: " + str(numgenes) + "\n" + \
                "Percent ZF genes in the human genome: " + percent_tot_zf_genes + "%\n" + \
                "Percent ZF genes associated with disease: " + percent_zf_genes_dis + \
                "%\nPercent ZF genes with allele information: " + \
                percent_zf_g_full_dis + "%\n\n" + \
                "The following lists the number of proteins that have\n" + \
                "the corresponding number of zinc-fingers. \n\n" + \
                '(Note: One "=" is equal to 2 proteins.)\n\n'
    msg += "Zfs: Proteins\n"
    for key in dict_domain:
        strvalue = str(dict_domain[key])
        while len(strvalue) < 3:
            strvalue = " " + strvalue
        strkey = str(key)
        while len(strkey) < 3:
            strkey = " " + strkey
        equals = "=" * (dict_domain[key] // 2)
        msg = msg + strkey + ": " + strvalue + "|" + equals + "\n"
    msg = msg + "\nThe average number of zinc-fingers per ZF-protein is " + \
        str(round(numdomains / numproteins, 2)) + \
        ".\n\n-------------------------------------\n\n"
    return msg


def basic_stats(gene_list):
    """
    Outputs basic statistics. File can be edited to write and do more complex
    statistics.

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes
    """
    genes_w_dis = genes_w_disease(gene_list)
    genes_w_full_dis = genes_w_full_dis_info(gene_list)
    len_gene_list = len(gene_list)
    print("Stats:")
    print("  " + str(genes_w_dis) + " ZF genes have diseases. ")
    print("  " + str(genes_w_full_dis) + " ZF genes w/ disease have nucleotide " +
          "information.")
    print("  There are " + str(len_gene_list) + " total ZF genes.")
    print("  " + str(round(genes_w_dis / len_gene_list * 100, 2)) +
          "% of ZF genes are associated with disease.")


def special_mutations(gene_list):
    """
    This procedure counts the number of special mutations.

    This function counts the number of terminations, ion-binding, and dna-binding
    mutations from a list of genes. One domain can have multiple mutations
    associated with it, so this function is used to make sure every mutation
    is counted

    Parameter gene_list: a list of Genes
    Preconditon: gene_list is a list of Genes
    """
    cnt = collections.Counter()
    for gene in gene_list:
        proteins = gene.protein_list()
        for protein in proteins:
            domains = protein.domain_list()
            for domain in domains:
                if domain.is_mutated:
                    dom_dict = domain.get_mutation_dictionary()
                    for mut_seq in dom_dict:
                        mut_type = dom_dict[mut_seq][0]
                        cnt[mut_type] += 1
                else:
                    raise helper.IncorrectError
    return cnt


def total_mutations(gene_list):
    """
    This procedure counts the total number of mutations.

    This procedure counts the total number of diseases with full disease info
    in a list of genes.

    Parameter gene_list: a list of Genes
    Preconditon: gene_list is a list of Genes
    """
    i = 0
    for gene in gene_list:
        diseases = gene.disease_list()
        for disease in diseases:
            if disease.has_full_info():
                i += 1
    return i


def count_genes_w_disease(gene_list):
    """
    This function takes in a list of genes and returns the number of genes
    that has a disease

    Parameter gene_list: a list of genes to get the IDs from.
    Precondtion gene_list is a list of Genes.

    Returns: the number of genes with disease[int]
    """
    assert type(gene_list) == list or type(gene_list) == set
    assert type(gene_list[0]) == genetics.Gene
    i = 0
    for gene in gene_list:
        if gene.has_disease():
            i += 1
    return i


def genes_w_full_dis_info(gene_list):
    """
    This function takes in a list of genes and returns the number of genes
    that have nucleotide information, amino acid change, and other information
    as well as the names of the disease.

    Parameter gene_list: a list of genes to get the IDs from.
    Precondtion gene_list is a list of Genes.

    Returns: the number of genes with full disease info[int]
    """
    assert type(gene_list) == list or type(gene_list) == set
    assert type(gene_list[0]) == genetics.Gene
    i = 0
    for gene in gene_list:
        if gene.has_disease_w_full_info():
            i += 1
    return i


def total_num_domains(alist):
    """
    This function calculates the total number of protein domains in a list
    of genes.

    Parameter alist: a list of genes.
    Precondition: alist is a list of genes.

    Returns: the total number of protein domains[int]
    """
    assert type(alist) == list
    i = 0
    for gene in alist:
        pro_list = gene.protein_list()
        for protein in pro_list:
            i += protein.num_domains()
    return i


def mutation_gene_stats(gene_list):
    """
    Outputs statistics about mutated gene list.

    This function gathers up information about a mutated list of genes. This
    info includes the number of mutations, the number in ion-binding region,
    the number in DNA-binding region, and other important statistics.
    """
    binding_region_count = special_mutations(gene_list)
    dna_count = binding_region_count['dna-binding']
    ion_count = binding_region_count['ion-binding']
    total = dna_count + ion_count + binding_region_count['none']
    return (dna_count, ion_count, total)


def protein_and_disease_stats(gene_list):
    """
    Outputs statistics about proteins and diseases.

    This function outputs the length of each ZF domain and the length of each
    protein sequence. It also counts the number of amino acids that mutations
    cause disease in. Used to perform statistical analysis.
    """
    pro = collections.Counter()
    dis = collections.Counter()
    mis = collections.Counter()
    amin_set = set()
    for gene in gene_list:
        proteins = gene.protein_list()
        diseases = gene.disease_list()
        for disease in diseases:
            amino = re.search(
                r"\D+(\d+)", disease.get_amino_change()).group(0)
            gene_name = disease.get_gene_name()
            if (gene_name, amino) not in amin_set:
                amin_set.add((gene_name, amino))
                dis['# of mutated amino acids'] += 1
        for protein in proteins:
            pro['# of amino acids'] += len(
                protein.prosequence())
            pro['# of proteins'] += 1
            domains = protein.domain_list()
            for domain in domains:
                if domain.is_valid():
                    pro['# of amino acids in ZF domains'] += len(
                        domain.domsequence())
                    pro['# of domains'] += 1
                    mut_dict = domain.get_mutation_dictionary()
                    mut_set = set()
                    for mut in mut_dict:
                        notes = mut_dict[mut][0]
                        new_amino = mut_dict[mut][2]
                        ind = mut_dict[mut][3]
                        if ind not in mut_set:
                            mut_set.add(ind)
                            dis['# of mutated amino acids in domains'] += 1
                            if new_amino != '*':
                                mis['# of MISSENSE amino acids in domains'] += 1
                            if 'binding' in notes:
                                dis['# of mutated amino acids in DNA/ion-binding regions'] += 1
                            if 'dna' in notes:
                                dis['# of mutated amino acids in DNA-binding regions'] += 1
                            if 'ion' in notes:
                                dis['# of mutated amino acids in ion-binding regions'] += 1
                            if 'hydrophobic' in notes:
                                dis['# of mutated amino acids in hydrophobic regions'] += 1
                            if 'binding' in notes and new_amino != '*':
                                mis['# of MISSENSE amino acids in DNA/ion-binding regions'] += 1
                            if 'dna' in notes and new_amino != '*':
                                mis['# of MISSENSE amino acids in DNA-binding regions'] += 1
                            if 'ion' in notes and new_amino != '*':
                                mis['# of MISSENSE amino acids in ion-binding regions'] += 1
                            if 'hydrophobic' in notes and new_amino != '*':
                                mis['# of MISSENSE amino acids in hydrophobic regions'] += 1
    return (dis, pro, mis)


def missense_mutation_stats(gene_list):
    """
    Outputs statistics about point mutations.

    This function outputs a collections Counter of mutation stats to later
    output.
    """
    cnt = collections.Counter()
    for gene in gene_list:
        diseases = gene.disease_list()
        protein = gene.protein_list()[-1]
        for disease in diseases:
            if disease.has_full_info():
                amino_change = disease.get_amino_change()
                change = re.search(r"\D{3}(\d+)(\D{3})", amino_change)
                Ter = change.group(2)
                cnt['Total mutations'] += 1
                if Ter != 'Ter':
                    cnt['Total missense mutations'] += 1
                else:
                    cnt['# of nonsense mutations'] += 1
                    domains = protein.domain_list()
                    mut_ind = int(change.group(1))
                    for domain in domains:
                        end_pos = domain.get_end_position()
                        if mut_ind <= end_pos - 1:
                            # print("Mutation Index: ", mut_ind)
                            # print("Domain End pos: ", end_pos)
                            # print('End Amino Acid: ',
                            #       protein.prosequence()[int(end_pos - 1)])
                            cnt['# of nonsense mutations eliminating a part of 1 domain'] += 1
                            break
                    for domain in domains:
                        start_pos = domain.get_start_position()
                        if mut_ind <= start_pos - 1:
                            cnt['# of nonsense mutations eliminating an entire domain'] += 1
                            break

    cnt['Missense mutations in DNA-binding region'] += 0
    for gene in gene_list:
        proteins = gene.protein_list()
        for protein in proteins:
            domains = protein.domain_list()
            for domain in domains:
                if domain.is_valid():
                    mut_dict = domain.get_mutation_dictionary()
                    for mutation in mut_dict:
                        notes = mut_dict[mutation][0]
                        disease = mut_dict[mutation][-1]
                        amino_change = disease.get_amino_change()
                        change = re.search(
                            r"\D{3}\d+(\D{3})", amino_change).group(1)
                        if 'dna' in notes and change != "Ter":
                            cnt['Missense mutations in DNA-binding region'] += 1
                        if 'ion' in notes and change != "Ter":
                            cnt['Missense mutations in ion-binding region'] += 1
                        if 'hydrophobic' in notes and change != "Ter":
                            cnt['Missense mutations in hydrophobic region'] += 1
    return cnt


def quantify_dna_specifity(foldername='analysis', analysis='euclidean'):
    """
    Quantifies the change in DNA-binding specficity.

    This function quantifies the change in DNA-binding specficity and calculates
    pearson's R for a list of mutations that affect the DNA-binding regions
    of zinc finger proteins.

    Parameter foldername: the folder the file is located in
    Preconditon: foldername is a string

    Parameter analysis: the measure used to analyze the data.
    Preconditon: analysis is 'euclidean' or 'pearsons'
    """
    assert type(analysis) == str and type(
        analysis) == str, 'Invalid Parameters'
    distance = ''.join([i for i in analysis if i.isalpha()])
    assert analysis == 'euclidean' or analysis == 'pearsons', "Distance can " + \
        'only be euclidean or pearsons.'
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    fileDir = os.path.join(fileDir, '../' + foldername)
    ori_list = glob.glob(fileDir + "/*o.txt")
    iteration = 0
    total_sum = 0
    adict = dict()
    arr = []
    for filename in sorted(ori_list):
        with open(filename, 'r') as f:
            original = f.read()
        a_ori = re.search(
            r"a\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", original)
        c_ori = re.search(
            r"c\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", original)
        g_ori = re.search(
            r"g\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", original)
        t_ori = re.search(
            r"t\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", original)
        mutation_list = glob.glob(filename[0:-5] + 'm*')
        for mutation_file in sorted(mutation_list):
            with open(mutation_file, 'r') as g:
                mutated = g.read()
            # print(mutation_file)
            a_mut = re.search(
                r"a\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", mutated)
            c_mut = re.search(
                r"c\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", mutated)
            g_mut = re.search(
                r"g\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", mutated)
            t_mut = re.search(
                r"t\s+([\d\.]+)\s+([\.\d]+)\s+([\.\d]+)\s+([\d\.]+)", mutated)
            a1_ori = float(a_ori.group(1))
            c1_ori = float(c_ori.group(1))
            g1_ori = float(g_ori.group(1))
            t1_ori = float(t_ori.group(1))
            a1_mut = float(a_mut.group(1))
            c1_mut = float(c_mut.group(1))
            g1_mut = float(g_mut.group(1))
            t1_mut = float(t_mut.group(1))

            a2_ori = float(a_ori.group(2))
            c2_ori = float(c_ori.group(2))
            g2_ori = float(g_ori.group(2))
            t2_ori = float(t_ori.group(2))
            a2_mut = float(a_mut.group(2))
            c2_mut = float(c_mut.group(2))
            g2_mut = float(g_mut.group(2))
            t2_mut = float(t_mut.group(2))

            a3_ori = float(a_ori.group(3))
            c3_ori = float(c_ori.group(3))
            g3_ori = float(g_ori.group(3))
            t3_ori = float(t_ori.group(3))
            a3_mut = float(a_mut.group(3))
            c3_mut = float(c_mut.group(3))
            g3_mut = float(g_mut.group(3))
            t3_mut = float(t_mut.group(3))

            a4_ori = float(a_ori.group(4))
            c4_ori = float(c_ori.group(4))
            g4_ori = float(g_ori.group(4))
            t4_ori = float(t_ori.group(4))
            a4_mut = float(a_mut.group(4))
            c4_mut = float(c_mut.group(4))
            g4_mut = float(g_mut.group(4))
            t4_mut = float(t_mut.group(4))

            matrix_ori = numpy.array([[a1_ori, c1_ori, g1_ori, t1_ori], [a2_ori, c2_ori, g2_ori, t2_ori],
                                      [a3_ori, c3_ori, g3_ori, t3_ori], [a4_ori, c4_ori, g4_ori, t4_ori]])
            matrix_mut = numpy.array([[a1_mut, c1_mut, g1_mut, t1_mut], [a2_mut, c2_mut, g2_mut, t2_mut],
                                      [a3_mut, c3_mut, g3_mut, t3_mut], [a4_mut, c4_mut, g4_mut, t4_mut]])

            if distance == 'euclidean':
                ED = scipy.spatial.distance.euclidean(
                    matrix_ori.flatten(), matrix_mut.flatten())
                arr.append(ED)
                #helper.easy_print("ED", ED)
            else:
                # print(mutation_file)
                ###############
                correlation_coefficient = scipy.stats.pearsonr(
                    matrix_ori.flatten(), matrix_mut.flatten())
                arr.append(1 - correlation_coefficient[0])
                #print(1 - correlation_coefficient[0])
                # if correlation_coefficient[0] < 0:
                #     print('negative')
    # print('...............')
    method = ''
    if distance == 'euclidean':
        method = 'ED'
    else:
        method == 'PCC'
    helper.easy_print("Mean " + method, statistics.mean(arr))
    helper.easy_print("Median " + method, statistics.median(arr))
    return arr


if __name__ == '__main__':
    quantify_dna_specifity(foldername='analysis', analysis='pearsons')
