"""
A file used to store databases.

This file stores functions that search a gene list through a database.

This function is used to build functio

Author: Austin Starks
Data: June 25, 2019
"""

import helper
import genetics
import csv
import re
import time
import random
from string import ascii_lowercase
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from bs4 import BeautifulSoup as soup
from fake_useragent import UserAgent


def search_mutation_info(gene_list):
    """
    This function searches databases to get information about diseases in a
    list of genes. It will only search databases that contain information on
    nucleotide positions

    Parameter gene_list: a list of gene.
    Preconditon: gene_list is a list of genes.
    """
    gene_name_list = genetics.Gene.gene_name_list(gene_list)
    gene_name_set = set(gene_name_list)
    CLINVAR_nuc(gene_list, gene_name_list, gene_name_set)
    UniProt(gene_list, gene_name_list, True)
    paper_search(gene_list, gene_name_list, gene_name_set)
    genetics.Disease.delete_repeated_info(gene_list)


def search(gene_list):
    """
    This function searches databases to get information about diseases in a
    list of genes.

    Parameter gene_list: a list of gene.
    Preconditon: gene_list is a list of genes.

    Parameter gene_name_list: a list of gene names.
    Preconditon: gene_name_list is a list of gene names.
    """
    gene_name_list = genetics.Gene.gene_name_list(gene_list)
    gene_name_set = set(gene_name_list)
    CGD(gene_list, gene_name_list, gene_name_set)
    CLINVAR(gene_list, gene_name_list, gene_name_set)
    OMIM(gene_list, gene_name_list, gene_name_set)
    NLM(gene_list, gene_name_list, gene_name_set)
    DisGeNet(gene_list, gene_name_list, gene_name_set)
    UniProt(gene_list, gene_name_list, False)


def CGD(gene_list, gene_name_list, gene_name_set):
    """
    This procedure searches the Clinical Genome Database (download to computer
    and convert to tab-delimited file) for information about disease

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names
    """
    print("Searching the Clinical genetics Database")
    filename = helper.edit_filename(
        "CGD.txt", "sibling", foldername='databases')
    with open(filename, "r", encoding='utf-8') as f:
        reader = csv.reader(f, dialect="excel", delimiter='\t')
        for row in reader:
            gene_name = row[0]
            if gene_name in gene_name_set:
                ind = gene_name_list.index(gene_name)
                gene = gene_list[ind]
                gene.set_disease(row[3], gene_name, "CGD")
    print("CGD database complete")


def CLINVAR(gene_list, gene_name_list, gene_name_set):
    """
    Searches the gene_condition_id database that is downloaded to computer.
    From ClinVar

    Source file from /pub/clinvar/

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names
    """
    print("Searching ClinVar")
    filename = helper.edit_filename(
        "gene_condition_source_id", "sibling", foldername='databases')
    with open(filename, "r", encoding='utf-8') as f:
        reader = csv.reader(f, dialect="excel", delimiter='\t')
        for row in reader:
            gene_name = row[1] + row[2]
            if gene_name in gene_name_set:
                ind = gene_name_list.index(gene_name)
                gene = gene_list[ind]
                gene.set_disease(row[4], gene_name, "ClinVar")
    print("ClinVar database complete")


def NLM(gene_list, gene_name_list, gene_name_set):
    """
    This procedure searches the GHR NLM database of genes and finds all of the
    genes that are in a list of genes.

    The database can be found here: https://ghr.nlm.nih.gov/gene
    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names

    Parameter start_time: the time this function began. This is used to calculate
    the amount of time this function takes to execute and print it on the screen.
    Preconditon: start_time is a time object.
    """
    url1 = "https://ghr.nlm.nih.gov/gene?initial="
    alphabet = list(ascii_lowercase)  # list of letters
    ua = UserAgent()  # Needed so the website won't block my acesss; generates
    # a thingy that looks just look regular chrome.
    this_list = []
    print("Searching NLM website. This takes a few seconds.")
    start_time = time.time()
    for let in alphabet:
        headers = {"User-Agent": ua.random}  # uses the thingy I generated
        request = Request(url1 + let, headers=headers)  # request website
        response = urlopen(request)  # open the website
        respData = response.read()  # read contents on the website
        response.close()  # close the website
        page_soup = soup(respData, "html.parser")  # get HTML from website
        raw_info = page_soup.findAll("ul", {"class": "browse-results"})[0].text
        # Useful function that I need to learn to use better.
        this_list = this_list + \
            re.findall("(" + let.upper() + "[A-Z\d]+): \w", raw_info)
        helper.time_elapsed(start_time)
    for gene_name in this_list:
        helper.time_elapsed(start_time)
        url2 = "https://ghr.nlm.nih.gov/gene/"
        if gene_name in gene_name_set:
            try:
                ind = gene_name_list.index(gene_name)
                gene = gene_list[ind]
                headers = {"User-Agent": ua.random}
                request = Request(url2 + gene_name +
                                  "#conditions", headers=headers)
                response = urlopen(request)
                respData = response.read()
                response.close()
                page_soup = soup(respData, "html.parser")
                dis1 = re.findall(
                    "Health Conditions Related to Genetic Changes\s+(([\w+,-]+ )+)", page_soup.text)
                disease = dis1[0][0]
                gene.set_disease(disease, gene_name, "NLM")
                diseases = re.findall("More About This Health Condition\s+(([\w+,-]+ )+)",
                                      page_soup.text)[0:-1]
                for disease in diseases:
                    gene.set_disease(disease[0], gene_name, "NLM")
            except IndexError:
                pass
    print("\nNLM database complete.")


def OMIM(gene_list, gene_name_list, gene_name_set):
    """
    This procedure searches the OMIM database (downloaded to computer)

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names
    """

    print("Searching the OMIM Database")
    filename = helper.edit_filename(
        "morbidmap.txt", "sibling", foldername='databases')
    with open(filename, "r", encoding='utf-8') as f:
        reader = csv.reader(f, dialect="excel", delimiter='\t')
        for row in reader:
            try:
                gene_row_list = row[1].split(", ")
                for gene_name in gene_row_list:
                    if gene_name in gene_name_set:
                        ind = gene_name_list.index(gene_name)
                        gene = gene_list[ind]
                        disease = row[0]
                        gene.set_disease(disease, gene_name, "OMIM")
            except IndexError:
                pass
    print("OMIM database complete")


def DisGeNet(gene_list, gene_name_list, gene_name_set):
    """
    This procedure searches DisGeNet (download to computer and
    and convert to tab-delimited file) for information about disease

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names
    """
    print("Searching DisGesNet")
    filename = helper.edit_filename("all_gene_disease_associations.tsv", "sibling",
                                    foldername='databases')
    with open(filename, "r", encoding='utf-8') as f:
        reader = csv.reader(f, dialect="excel", delimiter='\t')
        for row in reader:
            gene_name = row[1]
            if gene_name in gene_name_set:
                ind = gene_name_list.index(gene_name)
                gene = gene_list[ind]
                gene.set_disease(row[5], gene_name, "DisGeNet")
    print("DisGeNet database complete")


def UniProt(gene_list, gene_name_list, get_info=True):
    """
    This procedure searches UniProt human polymorphism and disease mutations
    (download to computer) for information about disease.

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names

    Parameter get_info: a boolean that says whether UniProt should search the
    web for dbSNP information
    Preconditon: get_info is a bool
    """
    if not get_info:
        print("Searching UniProt. This takes a few seconds.")
    else:
        print("Searching UniProt for mutation info. This takes a few minutes.")
    filename = helper.edit_filename(
        "humsavar.txt", "sibling", foldername='databases')
    with open(filename, "r") as g:
        data = g.read()
    start_time = time.time()
    i = 0
    for gene_name in gene_name_list:
        if '\n' + gene_name + " " in data:
            row_list = re.findall('\n' + gene_name + r" .+", data)
            helper.time_elapsed(start_time)
            for row in row_list:
                if "Disease" in row:
                    disease_name = re.search(r"(rs\w+|-)\s+(.+)", row).group(2)
                    ind = gene_name_list.index(gene_name)
                    gene = gene_list[ind]
                    gene.set_disease(disease_name, gene_name,
                                     "UniProt", get_info)
                    if get_info:
                        disease = gene.disease_list()[-1]
                        dbSNP = re.search(r"(rs\w+|-)", row).group(1)
                        aa_change = re.search(r"(p.\w+)", row).group(0)
                        disease.set_amino_change(aa_change)
                        disease.set_SNP(dbSNP)
                        disease.get_SNP_info()
                        if disease.has_full_info():
                            i = i + 1
                        if i == 6:
                            i = 0
                            time.sleep(random.uniform(0.5, 1.25))
                            helper.time_elapsed(start_time)
                            time.sleep(random.uniform(0.5, 1.25))
    print("\nUniProt database complete")


def CLINVAR_nuc(gene_list, gene_name_list, gene_name_set):
    """
    Searches the variant_summary database that is downloaded to computer.
    From ClinVar

    Source file from /pub/clinvar/

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names
    """
    print("Searching ClinVar")
    filename = helper.edit_filename(
        "variant_summary.txt", "sibling", foldername='databases')
    with open(filename, "r", encoding='utf-8') as f:
        reader = csv.reader(f, dialect="excel", delimiter='\t')
        for row in reader:
            gene_name = row[4]
            gene_name_in_set = gene_name in gene_name_set
            if gene_name_in_set:
                snv = row[1] == 'single nucleotide variant'
                pathogenic = row[6] == 'Pathogenic' or row[6] == 'Likely pathogenic'
                GRCh38 = "GRCh38" == row[16]
                germline = "germline" == row[15] or "germline/somatic" == row[15]
                if snv and pathogenic and GRCh38 and germline:
                    disease_name = row[13]
                    chromosome = row[18]
                    position = int(row[19])
                    nucleotide = re.search(r":c.+(\w>\w)", row[2]).group(1)
                    try:
                        protein = re.search(r"\((p.\w+).", row[2]).group(1)
                    except:
                        protein = ''
                    if protein != '':
                        ind = gene_name_list.index(gene_name)
                        gene = gene_list[ind]
                        gene.set_disease(
                            disease_name, gene_name, "ClinVar", True)
                        disease = gene.disease_list()[-1]
                        disease.set_amino_change(protein)
                        disease.set_allele_change(nucleotide)
                        disease.set_chromosome(chromosome)
                        disease.set_position(position, 'GRCh38')
    print("ClinVar database complete")


def paper_search(gene_list, gene_name_list, gene_name_set):
    """
    This procedure searches the database 'zf_database.txt' for mutations that
    affect ZF proteins.

    This database was created by looking at the list of ZF proteins that had
    mutation info, and looking for papers about those proteins. Each mutation
    is sourced to know what paper it came from.

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of genes

    Parameter gene_name_list: a list of gene names
    Preconditon: gene_name_list is a list of Gene names [str]

    Parameter gene_name_set: a set of gene names (used for searching)
    Preconditon: gene_name_set is a set of gene names
    """
    print("Searching created tab file of mutations")
    filename = helper.edit_filename('zf_database.txt', "sibling",
                                    foldername='databases')
    with open(filename, 'r', encoding='ISO-8859-1') as f:
        next(f)
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        for row in reader:
            gene_name = row[0]
            if gene_name in gene_name_set:
                ind = gene_name_list.index(gene_name)
                gene = gene_list[ind]
                gene.set_disease(row[3], gene_name, row[8], True)
                disease = gene.disease_list()[-1]
                disease.set_chromosome(row[4])
                position = int(row[5])
                notes = row[9]
                if 'Nucleotide # based on:' in notes:
                    type_pos = re.search(
                        r'Nucleotide # based on: ([^\s;]+)', notes).group(1)
                elif position < 20000:
                    type_pos = "RELATIVE"
                else:
                    type_pos = 'CHROMOSOME'
                disease.set_position(position, type_pos)
                disease.set_allele_change(row[6])
                disease.set_amino_change(row[7])
                # notes = row[9]
                # if 'NM' in notes:
                # allele_change = convert_to_chrom(allele_change)
    print("Created tab file database complete")
