"""
A file used to maintain classes and functions related to genetic data.

This file stores all of the functions that have to deal with biological data
such as information about disease, genes, and proteins. The bulk of these functions
have to deal with storing infomration in list of genes, proteins they encode for,
and diseases associated with each gene. Disease, Gene, and Protein are their own
classes and have their own static functions.

Author: Austin Starks
Data: June 17, 2019
"""

import copy
import re
import glob
import os.path
import subprocess
import csv
import time
import helper
import database
import collections
from ast import literal_eval
from mygene import MyGeneInfo  # pip install mygene
from Bio.SeqIO import parse  # Need to install BioPython
from Bio.SearchIO import to_dict  # Need to install BioPython
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from bs4 import BeautifulSoup as soup
from fake_useragent import UserAgent


def hmmer_output(command="hmmsearch", hmm="zf_C2H2.ls.hmm",
                 fasta="Homo_sapiens.GRCh38.pep.all.fa", filename='hmmer.txt', create_new=False):
    """
    This method runs HMMER and creates a textfile of its output. Then, it
    returns that text as a string.

    Parameter command: a command that outputs HMMER output as a text file
    Precondition: command is a string include a filename.txt in the output

    Parameter create_new: a boolean that will create a new HMMER output if
    the file exists already (if True). False otherwise

    Returns: the HMMER output as a string.
    """
    command = helper.edit_filename(
        "hmmer-2.3.2/src/" + command, "sibling", "hmmer")
    # This assumes you have hmmer in the same folder as the hmmer folder
    # Code may need to be edited on a different computer to get hmmer to run.
    # Just need to put hmmer files in that folder, and it should execute.
    hmm = helper.edit_filename(hmm, "sibling", "hmmer")
    fasta = helper.edit_filename(fasta, "sibling", "hmmer")
    filename = helper.edit_filename(filename, "sibling", "hmmer")
    if not os.path.exists(filename) or create_new:
        subprocess.run(command + " " + hmm + " " +
                       fasta + "> " + filename, shell=True)
    # open the file and turn the data into a string.
    with open(filename, 'r') as myfile:
        data = myfile.read()
    if data == "":
        raise AttributeError("hmmer.txt should not be blank. Check to be sure" +
                             " that HMMER is installed properly. You can do this by adding hmmer-2.3.2" +
                             "to the HMMER directory")
    if create_new:
        print("hmmer.txt created successfully.")
    return data


def amino_map():
    """
    Creates a dictionary that maps the 3 letter amino acid code to the 1 letter
    code (and vice-versa)

    Returns: a amino acid dictionary (map)
    """
    amino_dict = dict()
    # Add all of the nucleotides
    amino_dict['F'] = 'Phe'
    amino_dict['L'] = 'Leu'
    amino_dict['I'] = 'Ile'
    amino_dict['M'] = 'Met'
    amino_dict['V'] = 'Val'
    amino_dict['S'] = 'Ser'
    amino_dict['P'] = 'Pro'
    amino_dict['T'] = 'Thr'
    amino_dict['A'] = 'Ala'
    amino_dict['Y'] = 'Tyr'
    amino_dict['H'] = 'His'
    amino_dict['Q'] = 'Gln'
    amino_dict['N'] = 'Asn'
    amino_dict['K'] = 'Lys'
    amino_dict['D'] = 'Asp'
    amino_dict['E'] = 'Glu'
    amino_dict['C'] = 'Cys'
    amino_dict['W'] = 'Trp'
    amino_dict['R'] = 'Arg'
    amino_dict['S'] = 'Ser'
    amino_dict['G'] = 'Gly'
    amino_dict['Ter'] = "*"
    tmp_dict = dict()
    # Reiterate through the list to get the opposite
    for one in amino_dict:
        three = amino_dict[one]  # save the 3 letter code
        tmp_dict[three] = one  # 3 letter code maps to one letter code
        # lower case one letter code maps to three
        tmp_dict[one.lower()] = three
        # lower case three letter code maps to one
        tmp_dict[three.lower()] = one
    amino_dict['X'] = 'Ter'
    return {**amino_dict, **tmp_dict}  # return both dicts combined


def chromo_dict(foldername="chromosomes"):  # "test_chromosomes"
    """
    Creates a dictionary containing the nucleotide information in each chromosome.

    This function goes into the chromosome folder and reads all of the files.
    Then, it maps each chromosome number to the text in each file.

    Returns: a chromosome dictionary
    """
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    fileDir = os.path.join(fileDir, '../' + foldername)
    file_list = glob.glob(fileDir + "/*.fa")
    chromo_dict = dict()
    i = 0
    print("Creating a chromosome dictionary for the human body.")
    for filename in file_list:
        with open(filename, "r") as f:
            data = f.read()
            first_new_line = data.find('\n')
            num = re.search("(\d+|X|Y)", data[0:first_new_line]).group(0)
            try:
                num = int(num)
            except:
                pass
            data = data[first_new_line + 1:]
            data = data.replace("\n", "")
            chromo_dict[num] = data
            print("Have " + str(i) + " chromosome(s).")
            i = i + 1
    print("Chromosome dictionary created successfully")
    return chromo_dict


def opposite_nucleotide_dict():
    """
    This function returns a dictionary of the opposing nucleotide. A pairs with
    T and G pairs with C.

    Returns: a dictionary mapping of nucleotides
    """
    opp_nuc = dict()
    opp_nuc["A"] = "T"
    opp_nuc["a"] = "T"
    opp_nuc["T"] = "A"
    opp_nuc["t"] = "A"
    opp_nuc["G"] = "C"
    opp_nuc["g"] = "C"
    opp_nuc["C"] = "G"
    opp_nuc["c"] = "G"
    return opp_nuc


def verify_nucleotide_and_amino(gene_list):
    """
    Verifies that the nucleotide position is correct.

    This function takes in a list of genes and checks that the position in each
    disease is correct.

    Throws IncorrectError if the nucleotide position cannot be confirmed. If
    the amino acid position can't be confirmed, it removes that disease from
    the disease list. However, most that are removed are amino acids that have
    'Ter' in it, so that explains why that position is flagged.

    Parameter gene_list: a list of genes
    Preconditon: gene_list is a list of Genes.

    Returns: a list of diseases to remove.
    """
    chrom_dict = chromo_dict()
    opp_dict = opposite_nucleotide_dict()
    amap = amino_map()
    adict = dict()
    print("Verifying nucleotide and amino acid positions")
    for gene in gene_list:
        diseases = gene.disease_list()
        for disease in diseases:
            position = disease.get_position()
            chromosome = disease.get_chromosome()
            allele_change = disease.get_allele_change()
            direction = gene.get_direction()
            amino_change = disease.get_amino_change()
            verify_nuc = verify_nucleotide(chromosome, position, allele_change,
                                           direction, chrom_dict, opp_dict)
            verify_amino, codon = verify_amino_acid(chromosome, position, allele_change,
                                                    amino_change, direction, chrom_dict, amap, adict)
            if not verify_amino:
                adict[disease] = 'DNA -> RNA -> AA translation error\n' +\
                    "DNA sequence: 5'-(position-" + str(position - 2) + ") " + codon +\
                    "-(position " + str(position + 2) + ")-3'"
    return adict
    print("Nucleotide Positions verfied.")


def verify_nucleotide(chromo, position, allele_change, direction, chrom_dict,
                      nucleotide_pair):
    """
    Returns: True if the amino acid can be coded by the sequence. False otherwise.

    Parameter chromo: the chromosome number
    Preconditon: chromosome is an int between 1 and 22 or X or Y

    Parameter position: the position of the nucleotide change
    Preconditon: position is an int

    Parameter allele_change: the change in alleles
    Preconditon: allele_change is a 3 character string

    Parameter direction: the direction of the transcript
    Preconditon: direction is 1 or -1

    Parameter chrom_dict: the dictionary of chromosomes
    Preconditon: chrom_dict is a dictionary mapping of chromosomes

    Parameter nucleotide_pair: a dictionary mapping A to T and G to C (etc.)
    Preconditon: nucleotide_pair must be a dict
    """
    assert chromo is "X" or chromo is "Y" or type(chromo) == int, \
        "Chromosome should be an int, X, or Y."
    if type(chromo) == int:
        assert chromo >= 1 and chromo <= 22
    assert direction is 1 or direction is -1, "Invalid direction"
    assert type(position) == int, "Position must be an integer"
    assert type(nucleotide_pair) == dict, "Nucleotide Dict must be a dict"
    assert len(allele_change) == 3, "Allele change should be 3 characters"
    test_nuc = allele_change[0]
    correct_nuc = chrom_dict[chromo][position - 1]
    if direction == 1:
        verify_nuc = test_nuc == correct_nuc
        if not verify_nuc:
            print("Correct nuc: " + correct_nuc)
            print("Test nuc: " + test_nuc)
            print("Gene Direction: " + str(direction))
            print("Gene ID: " + gene.get_gene_id())
            raise helper.IncorrectError()
    elif direction == -1:
        test_nuc = nucleotide_pair[test_nuc]
        verify_nuc = correct_nuc == test_nuc
        if not verify_nuc:
            print("Correct nuc: " + correct_nuc)
            print("Test nuc: " + test_nuc)
            print("Gene Direction: " + str(direction))
            print("Gene ID: " + gene.get_gene_id())
            raise helper.IncorrectError()
    else:
        raise helper.IncorrectError("Test nuc: " + test_nuc +
                                    "Correct nuc: " + correct_nuc)
    return True


def verify_amino_acid(chromo, position, allele_change, amino_change, direction,
                      chrom_dict, amino_map, removal_dict):
    """
    Returns: True if the amino acid can be coded by the sequence. False otherwise.

    Parameter chromo: the chromosome number
    Preconditon: chromosome is an int between 1 and 22 or X or Y

    Parameter position: the position of the nucleotide change
    Preconditon: position is an int

    Parameter allele_change: the change in alleles
    Preconditon: allele_change is a 3 character string

    Parameter amino_change: the change in AA
    Preconditon: amino_change is a str

    Parameter direction: the direction of the transcript
    Preconditon: direction is 1 or -1

    Parameter chrom_dict: the dictionary of chromosomes
    Preconditon: chrom_dict is a dictionary mapping of chromosomes

    Parameter amino_map: a dictionary mapping amino acids to their 1 letter code
    Preconditon: amino_map must not be empty if correct_amino is length 3.

    Parameter removal_dict: a dictionary of mutations to remove.
    Preconditon: removal_dict is a dict
    """
    assert chromo is "X" or chromo is "Y" or type(chromo) == int, \
        "Chromosome should be an int, X, or Y."
    if type(chromo) == int:
        assert chromo >= 1 and chromo <= 22
    assert direction is 1 or direction is -1, "Invalid direction"
    assert type(chrom_dict) == dict, "chrom_dict must be a dictionary."
    assert type(amino_map) == dict, "amino_map must be a dictionary."
    assert type(position) == int, "Position must be an integer"
    seq_five = Seq(chrom_dict[chromo][position - 3:position + 2], generic_dna)
    if direction == -1:
        seq_five = seq_five.reverse_complement()
    seq_list = [seq_five[0:3]] + [seq_five[1:4]] + [seq_five[2:5]]
    amino_process = re.search(r"p.(\D{3})(\d+)(\D{3})", amino_change)
    if amino_process == None:
        print(amino_change)
    correct_amino = amino_map[amino_process.group(1)]
    target_amino = amino_map[amino_process.group(3)]
    new_allele = allele_change[-1]
    i = 0
    for sequence in seq_list:
        if sequence.translate()[0] == correct_amino:
            codon = str(seq_five[i:i + 3])
            mutated_codon = Seq(codon[0] + new_allele + codon[2])
            mutated_amino = mutated_codon.translate()
            if mutated_amino == target_amino:
                return (True, str(codon))
            mutated_codon = Seq(codon[0] + codon[1] + new_allele)
            mutated_amino = mutated_codon.translate()
            if mutated_amino == target_amino:
                return (True, str(codon))
            mutated_codon = Seq(new_allele + codon[1] + codon[2])
            mutated_amino = mutated_codon.translate()
            if mutated_amino == target_amino:
                return (True, str(codon))
        i = i + 1
    seq_five = seq_five.reverse_complement()
    seq_list = [seq_five[0:3]] + [seq_five[1:4]] + [seq_five[2:5]]
    i = 0
    for sequence in seq_list:
        if sequence.translate()[0] == correct_amino:
            codon = str(seq_five[i:i + 3])
            mutated_codon = Seq(codon[0] + new_allele + codon[2])
            mutated_amino = mutated_codon.translate()
            if mutated_amino == target_amino:
                return (True, str(codon))
            mutated_codon = Seq(codon[0] + codon[1] + new_allele)
            mutated_amino = mutated_codon.translate()
            if mutated_amino == target_amino:
                return (True, str(codon))
            mutated_codon = Seq(new_allele + codon[1] + codon[2])
            mutated_amino = mutated_codon.translate()
            if mutated_amino == target_amino:
                return (True, str(codon))
        i = i + 1
        # print("...")
        # helper.easy_print("Gene name", gene.get_gene_name())
        # helper.easy_print("Codon", codon)
        # helper.easy_print("target_amino   ", target_amino)
        # mutated_amino = Seq(codon[0] + new_allele + codon[2]).translate()
        # helper.easy_print("mutated_amino 1", mutated_amino)
        # mutated_amino = Seq(codon[0] + codon[1] + new_allele).translate()
        # helper.easy_print("mutated_amino 2", mutated_amino)
        # mutated_amino = Seq(new_allele + codon[1] + codon[2]).translate()
        # helper.easy_print("mutated_amino 3", mutated_amino)
        # print("...")
    return (False, str(seq_five))


def verify_protein_sequence(gene_list):
    """
    Verifies that the protein sequence amino acid change is possible.

    This function takes in a list of genes, check its list of proteins, and then
    verifies that the sequence is possible to mutate.

    Parameter gene_list: a list of genes.
    Preconditon: gene_list is a list of Genes.

    Example of removal:
    zf_genes_w_disease.tsv contains an amino acid change (p.Arg112Ter)
    If the protein sequence doesn't have Arg at its 112th position, then that
    disease is removed.

    Returns: a set of diseases (to remove and get info from)
    """
    remove_dict = dict()
    amap = amino_map()
    for gene in gene_list:
        proteins = gene.protein_list()
        diseases = gene.disease_list()
        for disease in diseases:
            amino_change = re.search(
                r"p.(\D{3})(\d+)(\D{3})", disease.get_amino_change())
            amino_ind = int(amino_change.group(2))
            check_amino = amino_change.group(1)
            resulting_amino = amino_change.group(3)
            for protein in proteins:
                seq = protein.prosequence()
                try:
                    target_amino = amap[seq[amino_ind - 1]]
                except IndexError:
                    target_amino = "No amino acid"
                try:
                    pos_min_1 = amap[seq[amino_ind - 2]]
                except IndexError:
                    pos_min_1 = "No amino acid"
                try:
                    pos_plus_1 = amap[seq[amino_ind]]
                except IndexError:
                    pos_plus_1 = "No amino acid"
                # if target_amino != check_amino and pos_plus_1 == check_amino:
                #     new_change = 'p.' + check_amino + str(amino_ind + 1) + \
                #         resulting_amino + '?'
                #     disease.set_amino_change(new_change)
                if target_amino != check_amino:
                    remove_dict[disease] = 'Amino acid position mismatch' + \
                        "\n       " + pos_min_1 + \
                        " is at position " + str(amino_ind - 1) + \
                        "\n       " + target_amino + \
                        " is at position " + str(amino_ind) + \
                        "\n       " + pos_plus_1 + \
                        " is at position " + str(amino_ind + 1)
    return remove_dict


def remove_bad_diseases(gene_list, removal):
    """
    This procedure removes the diseases in removal from the gene_list.

    Parameter gene_list: a list of genes.
    Preconditon: gene_list is a list of Genes.

    Parameter removal: a set of diseases to remove
    Preconditon removal is a set of Diseases
    """
    length = str(len(removal))
    new_list = copy.deepcopy(gene_list)
    print("Removing unverifiable diseases.")
    for gene in new_list:
        diseases = gene.disease_list()
        for disease in diseases:
            if disease in removal:
                gene.remove_disease(disease)
    print("   A total of " + length + " diseases were removed.")
    return new_list


def gene_list_process_mutations(gene_list, console=False):
    """
    Processes the mutation stored in disease list

    This function takes the amino acid change in the list of diseases, and
    processes them. It returns a list of genes and notes in each gene
    if a mutation hit a core amino acid, ion-binding domain, or neither.

    Parameter gene_list: a list of Genes
    Preconditon: gene_list is a list of Genes

    Parameter console: a boolean that says whether to print the output to the
    console.
    Preconditon: console is a boolean

    Parameter adict: a dictionary of possible notes
    Preconditon: adict is a dictionary
    """
    print('Mutating list of genes')
    amap = amino_map()
    for gene in gene_list:
        proteins = gene.protein_list()
        diseases = gene.disease_list()
        for protein in proteins:
            domains = protein.domain_list()
            for domain in domains:
                start_pos = domain.get_start_position()
                end_pos = domain.get_end_position()
                new_seq = domain.domsequence()
                seq = domain.domsequence()
                for disease in diseases:
                    new_seq = domain.domsequence()
                    amino_change = disease.get_amino_change()
                    process_change = re.search(
                        r'p.(\w{3})(\d{0,4})(\w{3})', amino_change)
                    mutation_ind = process_change.group(2)
                    mutation_ind = int(mutation_ind)
                    mutant_amino = amap[process_change.group(3)]
                    if mutation_ind >= start_pos and mutation_ind <= \
                            end_pos and domain.is_valid():
                        current_amino_ind = mutation_ind - start_pos
                        current_amino = seq[current_amino_ind]
                        three_amino = process_change.group(1)
                        verify_amino = amap[three_amino]
                        if current_amino.upper() == verify_amino.upper():
                            new_seq = seq[0:current_amino_ind] + mutant_amino + \
                                seq[current_amino_ind + 1:]
                            if console:
                                helper.easy_print(
                                    "\nProtein ID  ", protein.get_protein_id())
                                helper.easy_print(
                                    "Domain num  ", domain.get_domain_number())
                                helper.easy_print(
                                    "Domain Score", domain.get_score())
                                helper.easy_print("Start pos   ", start_pos)
                                helper.easy_print("Mutation Ind", mutation_ind)
                                helper.easy_print("End pos     ", end_pos)
                                helper.easy_print("Original seq", seq)
                                helper.easy_print("Mutant seq  ", new_seq)
                                helper.easy_print(
                                    "orig amino  ", current_amino)
                                helper.easy_print("Mutant amino", mutant_amino)
                                print('\n')
                            if len(new_seq) != len(seq):
                                raise helper.IncorrectError
                        else:
                            print("\nProtein: " + protein.get_protein_id())
                            print("Gene: " + gene.get_gene_name())
                            print("Mutation: " + str(mutation_ind))
                            print("Index: " + str(current_amino_ind))
                            print("Start position: " + str(start_pos))
                            print("End position: " + str(end_pos))
                            print("Sequence: " + seq)
                            print("Current amino: " + current_amino)
                            print("Verify amino: " + verify_amino)
                            print("Amino change: " + amino_change)
                            print("Disease source: " + disease.get_source())
                            raise helper.IncorrectError("There " +
                                                        "was an error" +
                                                        " verifying " +
                                                        "amino.")
                    if seq != new_seq:
                        disease.set_hit_domain(True)
                        domain.set_mutated()
                        domain.get_mutated_region_info(
                            disease, new_seq, "none")
                        # gene.remove_disease(disease)
    print("List of genes mutated.")


class Protein(object):
    """
    A class representing a protein.

    This class maintains all information about a protein, including the domains
    the protein has.

    INSTANCE ATTRIBUTES:
        _id:        the protein ID [str]
        _info:      information about each domain in the protein [str arr]
        _domains:   a list of domains [Domain arr]
        _prosequence: the amino acid code for this protein.
    """

    def __init__(self, id, domains=[]):
        """
        Initializes a protein.

        Parameter protein_id: the ID of the protein starting with ENSP
        Precondition: protein_id is a valid ID for a protein.

        Parameter domains: a list of protein domains
        Precondition: domains is a list
        """
        self._id = id
        self._info = []
        self._domains = []
        if domains != []:
            self._domains = domains
        self._sequence = ''

    def __eq__(self, other):
        """
        Two proteins are equal if they have the same ID.

        Returns: True if this protein is equal to another protein. False
        otherwise.
        """

        return len(self._sequence) == len(other._sequence)

    def __hash__(self):
        """
        Ensures two equivalent proteins hash to the same value.
        """
        return hash(len(self._sequence)) + hash(self._id)

    def __lt__(self, other):
        """
        One protein is less than another gene if it hashes to a lower value.
        Used for sorting list of genes.
        """
        return len(self._sequence) < len(other._sequence)

    def get_protein_id(self):
        """
        Returns: the protein ID for this protein.
        """
        return self._id

    def info(self):
        """
        Returns: the protein info as a list
        """
        return self._info

    def prosequence(self):
        """
        Returns: the fasta information for this protein including the position
        and the amino acid sequence
        """
        return self._sequence

    def domain_list(self):
        """
        Returns: the protein domains in this protein
        """
        return self._domains

    def num_domains(self):
        """
        Returns: the number of protein domains this protein has.
        """
        i = 0
        for domain in self._domains:
            if domain.is_valid():
                i = i + 1
        return i

    def set_info(self, info):
        """
        Stores raw information about the protein domains as a string
        """
        self._info = info

    def set_prosequence(self, sequence):
        """
        Stores raw information about the protein domains as a string
        """
        self._sequence = sequence

    def core_sequences(self, newline=True):
        """
        Returns: the core sequences seperated by a "+" as a string.
        """
        core = []
        for domain in self._domains:
            if domain.is_valid():
                core.append(domain.core())
            else:
                core.append("????")
        core_str = "+".join(core)
        if newline:
            i = 0
            j = 10
            core_sep_num = len(core_str) // j
            # Insert a \n at pos 40, 81, 122...
            while core_sep_num >= 0:
                core_str = core_str[:j + i] + "\n" + core_str[j + i:]
                i = i + 1
                j = j + 10
                core_sep_num = core_sep_num - 1
        return core_str

    def set_domain(self, start_time):
        """
        Takes the data in info and stores the information in the domains

        Parameter start_time: the time this function began. This is used to calculate
        the amount of time this function takes to execute and print it on the screen.
        Preconditon: start_time is a time object.
        """
        for info in self._info:
            helper.time_elapsed(start_time)
            match = re.search(r"domain (\d+)", info)
            dom_num = int(match.group(1))
            check = re.search(r"\*->(.+)<-\*", info)
            ind_h1 = check.group(1).find("Hlrt")
            seq = re.search("\d{1,4}\s{4}([\w-]+)\s{4}\d{2,4}", info).group(1)
            c1 = seq[ind_h1 - 7]
            c2 = seq[ind_h1 - 5]
            c3 = seq[ind_h1 - 4]
            c4 = seq[ind_h1 - 1]
            core = c1 + c2 + c3 + c4
            check = re.search(r"\*->(.+)<-\*", info)
            ind_c1 = re.search(r"C\.*p\.*f", check.group(1)).span()[0]
            first_cys = seq[ind_c1]
            # helper.easy_print("First Cys", first_cys)
            first_cys_bool = first_cys == 'C'
            ind_c2 = re.search(
                r"C\.*g\.*k\.*s", check.group(1)).span()[0]
            second_cys = seq[ind_c2]
            # helper.easy_print("Sec Cys", second_cys)
            second_cys_bool = second_cys == 'C'
            ind_h2 = ind_h1 + 4
            second_his_bool = seq[ind_h2]
            # helper.easy_print("Second Cys is C", second_cys_bool)
            cys_is_present = first_cys_bool and second_cys_bool
            his_is_present = seq[ind_h1] == "H" and second_his_bool
            if his_is_present and cys_is_present:
                self._domains.append(
                    Domain(self._id, dom_num, seq, core))
            else:
                self._domains.append(
                    Domain(self._id, dom_num, seq, "????"))
        self._domains = sorted(self._domains)

    def set_domain_positions_sequence(self):
        """
        This function gets the starting and ending position for each domain in
        the list of domains for this protein.
        """
        for domain in self._domains:
            for info in self._info:
                from_ = re.search(
                    r", from (\d{0,4}) to (\d{0,4}): score (\d{1,2}\.\d)", info)
                start_pos = from_.group(1)
                end_pos = from_.group(2)
                dom_num = int(re.search(r": domain (\d{1,2})", info).group(1))
                score = float(from_.group(3))
                seq = re.search(
                    "\d{1,4}\s{4}([\w-]+)\s{4}\d{2,4}", info).group(1)
                if domain.get_domain_number() == dom_num:
                    domain.set_start(int(start_pos))
                    domain.set_end(int(end_pos))
                    domain.set_domsequence(seq)
                    domain.set_score(score)

    def all_DNA_ion_bind_regions(self):
        """
        This procedure creates a mapping of DNA and ion-binding region for all
        domains in this protein. This is used to later get information on
        population variatio of mutations in the DNA/ion-binding region.
        """
        for domain in self._domains:
            for info in self._info:
                from_ = re.search(
                    r", from (\d{0,4}) to (\d{0,4}): score (\d{1,2}\.\d)", info)

    @staticmethod
    def create_protein_set(data):
        """
        This function uses regular expression search to find all of the proteins IDs
        in a string. Proteins begin with ENSP and are followed by a series of
        digits, then a period, then 1 or 2 digits.

        Parameter data: a string that contains protein IDs
        Precondition: data is a string

        Returns: a set of protein IDs.
        """
        # First characters are 'ENSP'
        # Followed by at least 7 digits
        # Followed by more digits or a period.
        # Followed by one or more digits
        protein_set = set(re.findall("ENSP\d{7}[\d\.]+", data))
        print("Protein ID set created successfully")
        return protein_set

    @staticmethod
    def gene_list_set_sequence(gene_list, full_dict, fasta="Homo_sapiens.GRCh38.pep.all.fa"):
        """
        This procedure takes in a list of genes and looks at the proteins in the
        list. It then stores information in a fasta file for each protein.

        Parameter gene_list: a list of genes
        Preconditon: gene_list is a list of Genes

        Parameter full_dict: a full dictionary of proteins
        Preconditon: full_dict is a dictionary of proteins

        Parameter fasta: a fasta file
        Preconditon: fasta must be a fasta file
        """
        assert type(gene_list) == list, "Gene list must be a list"
        assert type(gene_list[0]) == Gene, "Gene list must be a list of Genes"
        assert type(fasta) == str, "fasta file must be a fasta file"
        assert "fa" in fasta or "fasta" in fasta, "fasta must be a fasta file."
        assert type(full_dict) == dict, "Full dict must be a dict"
        for gene in gene_list:
            proteins = gene.protein_list()
            for protein in proteins:
                sequence = full_dict[protein.get_protein_id()].seq
                protein.set_prosequence(sequence)

    @staticmethod
    def full_protein_dict(fasta="Homo_sapiens.GRCh38.pep.all.fa", print_=True):
        """
        This procedure parses a FASTA file

        Parameter fasta: the fasta file
        Precondtion: fasta is a string name for a fasta file

        Parameter print_: says whether to print to console
        Preconditon: print_ is a bool

        Returns: a dictionary of the protein ID and the decription mapped to a
        description
        """
        fasta = helper.edit_filename(fasta, "sibling", "hmmer")
        results = parse(fasta, "fasta")
        adict = to_dict(results)
        if print_:
            print("Dictionary of proteins created from FASTA file.")
        return adict

    @staticmethod
    def get_zf_protein_info(data, alist, set_domain=True, print_=True):
        """
        This function takes in a list of genes and from those genes, it gets the
        list of proteins, and finds important information about the proteins
        and stores the info in each protein.

        Parameter data: the data to find the information from.
        Precondition: data is HMMER output as a string.

        Parameter alist: a list of genes.
        Precondition: alist is a list of genes.

        Parameter set_domain: whether or not to execute set_domain (True)
        Precondtion: set_domain is a boolean

        Parameter print_: says whether to print to console.
        Preconditon: print_ is a bool
        """
        if print_:
            print("Creating the protein domains. This takes a minute.")
        start_time = time.time()
        for gene in alist:
            pro_list = gene.protein_list()
            for pro in pro_list:
                info = re.findall(
                    pro.get_protein_id() + r": domain.+\n.+\n.+\n.+\n.+\d+", data)
                pro.set_info(info)
                if print_:
                    helper.time_elapsed(start_time)
                if set_domain:
                    pro.set_domain(start_time)
        if print_:
            print("")
            print("Protein domains stored in proteins.")

    @staticmethod
    def get_proteins(gene_list, data, get_domain_info=True,
                     file="zf_protein_tabs.tsv", print_=True):
        """
        This function gets the proteins that a list of genes code for.

        Parameter data: the data to find the information from.
        Precondition: data is HMMER output as a string.

        Parameter gene_list: a list of genes.
        Preconditon: gene_list is a list of genes.

        Parameter file: the file containing the proteins.
        Preconditon: file is a filename

        Parameter print_: says whether to print to console
        Preconditon: print_ is a bool

        Returns: a list of genes with proteins encoded.
        """
        if print_:
            print("Getting proteins for a list of genes.")
        file = helper.edit_filename(file, "sibling", "output")
        gene_id_list = Gene.gene_id_list(gene_list)
        gene_id_set = set(gene_id_list)
        if not os.path.exists(file):
            raise helper.DoesNotExistError(
                "Gene info tab with protein info does not exist.")
        with open(file, "r", encoding='utf-8') as f:
            reader = csv.reader(f, dialect="excel", delimiter='\t')
            for row in reader:
                gene_id = row[1]
                if gene_id in gene_id_set:
                    gene = Gene.create_gene_from_row(row)
                    ind = gene_list.index(gene)
                    gene_list[ind].addprotein(gene.protein_list()[0])
        if get_domain_info:
            Protein.get_zf_protein_info(
                data, gene_list, set_domain=False, print_=print_)
        if print_:
            print("Proteins obtained successfully.")

    @staticmethod
    def remove_threshold(gene_list):
        """
        Filters proteins below the score threshold out of a gene list.

        This function is used to remove proteins from a gene list. It is used
        when the gene list is storing more than one protein ID.

        Parameter gene_list: a list of Genes
        Preconditon: gene_list is a list of Genes.
        """
        assert type(gene_list) == list, "Gene list must be a list"
        assert len(gene_list) > 0, "Gene list must be a list of Genes."
        assert type(gene_list[0]) == Gene, "Gene list must be a list of genes"
        print("Filtering low-scoring proteins.")
        remove_gene_list = []
        for gene in gene_list:
            remove_protein_list = []
            proteins = gene.protein_list()
            for protein in proteins:
                domains = protein.domain_list()
                abool = True
                for domain in domains:
                    if domain.is_valid(True):
                        abool = False
                        break
                if abool:
                    remove_protein_list.append(protein)
            for pro in remove_protein_list:
                proteins.remove(pro)
            if len(proteins) == 0:
                remove_gene_list.append(gene)
        for gene in remove_gene_list:
            gene_list.remove(gene)
        print("Low scoring proteins removed.")

    @staticmethod
    def remove_proteins(gene_list, score_threshold=False):
        """
        Filters proteins from a gene list.

        This function removes all of the proteins in a gene except the largest protein
        This is used to make data processing easier.

        Parameter gene_list: a list of genes
        Preconditon: gene_list is a list of Genes.

        Parameter score_threshold: says whether to get only output proteins with a
        score above 17.7
        Preconditon: score_threshold is a bool

        Parameter remove_all: says whether to remove proteins with a low score
        threshold even if it is a known ZF protein.
        Preconditon: remove_all is a bool
        """
        assert type(gene_list) == list, "Gene list must be a list"
        assert len(gene_list) > 0, "Gene list must be a list of Genes."
        assert type(gene_list[0]) == Gene, "Gene list must be a list of genes"
        print("Removing all except target protein from genes")
        amap = amino_map()
        c = collections.Counter()
        remove_list = []
        for gene in gene_list:
            proteins = sorted(gene.protein_list())
            pro_choice = proteins[-1]
            diseases = gene.disease_list()
            for protein in proteins:
                seq = protein.prosequence()
                for disease in diseases:
                    try:
                        amino_change = re.search(
                            r"p.(\D{3})(\d+)(\D{3})", disease.get_amino_change())
                        amino_ind = int(amino_change.group(2))
                        amino_check = amap[amino_change.group(1)]
                        real_amino = seq[amino_ind - 1]
                        try:
                            if real_amino == amino_check:
                                # helper.easy_print("Success", i)
                                c[protein] += 1
                        except IndexError:
                            pass
                    except:
                        pass
            if len(c) > 0:
                pro_choice = c.most_common(1)[0][0]
                c.clear()
            domains = pro_choice.domain_list()
            abool = False
            for domain in domains:
                if domain.is_valid(score_threshold):
                    abool = True
                    break
            if abool:
                gene.clear_proteins(pro_choice)
            else:
                remove_list.append(gene)
        for gene in remove_list:
            gene_list.remove(gene)
        print("Proteins removed successfully")

    @staticmethod
    def remove_protein_info(gene_list, score_threshold=False):
        """
        Returns filtered out proteins from a gene list.

        This function returns the gene list of the proteins that would be removed.
        This is done to gather statistics about the proteins

        Parameter gene_list: a list of genes
        Preconditon: gene_list is a list of Genes.

        Parameter score_threshold: says whether to get only output proteins with a
        score above 17.7
        Preconditon: score_threshold is a bool
        """
        assert type(gene_list) == list, "Gene list must be a list"
        assert len(gene_list) > 0, "Gene list must be a list of Genes."
        assert type(gene_list[0]) == Gene, "Gene list must be a list of genes"
        print("Filtering the proteins")
        amap = amino_map()
        c = collections.Counter()
        remove_list = []
        for gene in gene_list:
            proteins = sorted(gene.protein_list())
            pro_choice = proteins[-1]
            diseases = gene.disease_list()
            for protein in proteins:
                seq = protein.prosequence()
                for disease in diseases:
                    try:
                        amino_change = re.search(
                            r"p.(\D{3})(\d+)(\D{3})", disease.get_amino_change())
                        amino_ind = int(amino_change.group(2))
                        amino_check = amap[amino_change.group(1)]
                        real_amino = seq[amino_ind - 1]
                        try:
                            if real_amino == amino_check:
                                # helper.easy_print("Success", i)
                                c[protein] += 1
                        except IndexError:
                            pass
                    except:
                        pass
            if len(c) > 0:
                pro_choice = c.most_common(1)[0][0]
                c.clear()
            domains = pro_choice.domain_list()
            abool = False
            for protein in proteins:
                for domain in domains:
                    if domain.is_valid(score_threshold):
                        abool = True
                        break
            if abool:
                gene.clear_proteins(pro_choice)
            else:
                remove_list.append(gene)
        return remove_list


class Domain(object):
    """
    A class representing a protein domain.

    This class maintains all information about a protein domain. Each domain is
    either a valid or invalid ZF domain, and if it is valid, it has a core
    sequence. Each domain also has a number.

    INSTANCE ATTRIBUTES:
        _pro_id:    the protein ID for which this domain belongs to.
        _num:       the domain number [int]
        _valid:     the domain is valid if it has the standard C2H2 region.
                    invalid otherwise [bool]
        _seq:       the amino acid sequence coding this domain [str]
        _core:      the 4 amino acid core sequence of a ZF protein [str]
        _start_pos: the starting position of this domain [int]
        _end_pos:   the end position of this domain [int]
        _mutated:   a bool that says whether this domain has been mutated.
        _score:     the Domain score from HMMER output [float]
        _mut_seq:   an array of mutated sequences [str arr]
        _mut_dict:  a dictionary mapping each mutation to other important info.
        _key_positions: a dictionary mapping the key position indices

    """

    def __init__(self, pro_id, num, seq, core, valid=True):
        """
        Initializes a domain.

        Parameter num: the domain number
        Precondition: num is an int

        Parameter seq: the AAs coding for this domain
        Precondition: seq is a string

        Parameter valid: says whether this is a valid protein domain
        Precondition: valid is a boolean
        """
        assert type(pro_id) == str, "pro_id must be a string"
        assert type(num) == int, "num must be an int"
        assert type(seq) == str, "seq must be a string"
        assert type(core) == str, "core must be a string"
        assert type(valid) == bool, "valid must be a boolean"
        assert len(core) == 4, "length of core must be 4"
        self._pro_id = pro_id
        self._num = num
        self._valid = valid
        self._seq = seq
        self._core = core
        self._start_pos = -1
        self._end_pos = -1
        self._mutated = False
        self._score = -1
        self._mut_dict = dict()
        self._key_positions = dict()

    def __eq__(self, other):
        """
        Two domains are equal if they have the same domain number and are within
        the same protein.

        Returns: True if this domain is equal to another domain. False
        otherwise.
        """
        return self._num == other._num and self._pro_id == other._pro_id

    def __lt__(self, other):
        """
        This domain is less than another domain if this num is less than the
        other domain's num. Two domains can only be compared if they are
        domains of the same protein ID.

        Returns: True if this domain is less than another domain. False
        if this domain is greater than another domain. NotImplemented if the
        domains cannot be compared.
        """
        if self._pro_id == other._pro_id:
            return self._num < other._num
        else:
            return NotImplemented

    def __hash__(self):
        """
        Ensures two equivalent domains hash to the same value.
        """
        return hash(self._pro_id + str(self._num))

    def get_domain_number(self):
        """
        Returns: the domain number as an int.
        """
        return self._num

    def domsequence(self):
        """
        Returns: the domain sequence
        """
        return self._seq.replace("-", "")

    def core(self):
        """
        Returns: the core sequence of this domain as a string. If this domain
        is not valid, returns "????"
        """
        if self._valid:
            return self._core
        else:
            return "????"

    def get_score(self):
        """
        Returns: the domain score. Confident scores are above 17.7
        """
        return self._score

    def set_score(self, score):
        """
        Stores the score for this domain. Confident domains have a score of above 17.7.
        """
        assert type(score) == float, "Score must be a float."
        self._score = score

    def set_mutated(self, mutated=True):
        """
        Sets this domain to become mutated (or not)
        """
        assert type(mutated) == bool
        self._mutated = mutated

    def is_mutated(self):
        """
        Returns: True if this is a valid ZF protein domain. False otherwise.
        """
        return self._mutated

    def is_valid(self, score_threshold=False):
        """
        Returns: True if this is a valid ZF protein domain. False otherwise.

        Parameter score_threshold: says whether to implement a minimum score
        threshold. In other words, if True, will only return domains that are
        confident ZF domains. False will return all ZF domains from HMMER.

        Preconditon: score_threshold is a boolean
        """
        assert type(score_threshold) == bool, "Score theshold must be a bool"
        if score_threshold:
            return self._core != "????" and self._score >= 17.7
        else:
            return self._core != "????"

    def get_start_position(self):
        """
        Returns: the starting position for this domain
        """
        return self._start_pos

    def get_end_position(self):
        """
        Returns: the ending position for this domain
        """
        return self._end_pos

    def set_domsequence(self, seq):
        """
        This procedure sets this domain's sequence to be seq

        Parameter seq: the sequence to set this domain to.
        Precondtion: seq is a string
        """
        self._seq = seq

    def set_start(self, start_pos):
        """
        This procedure adds the starting position to this domain.

        Parameter start_pos: the starting position
        Precondition: start_pos is an int.
        """
        self._start_pos = start_pos

    def set_end(self, end_pos):
        """
        This procedure adds the ending position to this domain.

        Parameter end_pos: the end position
        Precondition: end_pos is an int.
        """
        self._end_pos = end_pos

    def get_mutated_region_info(self, disease, mut_seq, notes="none"):
        """
        This function takes in a list of Genes and gets information about the
        mutated regions (such as if it's in an ion-binding region).

        Parameter disease: a disese to store in the mutation dictionary
        Preconditon disease is a disease

        Parameter mut_seq: a string that display the mutated sequence
        Preconditon: mut_seq is a string
        """
        assert isinstance(disease, Disease), "disease must be a Disease"
        assert type(mut_seq) == str, "Mutated sequence must be a string"
        if self._mutated:
            core = self._core
            sequence = self.domsequence()
            minus_one = core[0]
            two = core[1]
            three = core[2]
            six = core[3]
            match = re.search(minus_one + r"\w" + two +
                              three + r"\w{2}" + six, sequence)
            match_span = match.span()
            match_group = match.group(0)
            h1_ind = match_span[1]
            #    helper.easy_print('match', match)
            h2_ind = sequence.rfind('H')
            minus_one_ind = match_span[0]
            two_ind = minus_one_ind + 2
            three_ind = two_ind + 1
            six_ind = three_ind + 3
            new_match = re.search(
                r"(C)\w+(C)\w+" + match_group, sequence)
            c1_ind = new_match.span()[0]
            pattern = re.compile(r"(C)\w+")
            third_match = pattern.search(sequence, c1_ind + 1)
            c2_ind = third_match.span()[0]
            hydro_1_ind = six_ind - 2
            hydro_2_ind = minus_one_ind - 2
            for z in range(len(sequence)):
                if sequence[z] != mut_seq[z]:
                    mut_ind = z
            # print("...")
            # helper.easy_print("Pattern", pattern)
            # helper.easy_print("third match", third_match)
            # print("...")
            # helper.easy_print("C1 ind", c1_ind)
            # helper.easy_print("C1", sequence[c1_ind])
            # print("...")
            # helper.easy_print("C2 ind", c2_ind)
            # helper.easy_print("C2", sequence[c2_ind])
            ion_binding = c1_ind == mut_ind or c2_ind == mut_ind or h1_ind \
                == mut_ind or h2_ind == mut_ind
            dna_binding = minus_one_ind == mut_ind or two_ind == mut_ind or \
                three_ind == mut_ind or six_ind == mut_ind
            termination = mut_seq[mut_ind] == "*"
            other_structure = mut_ind == hydro_1_ind or mut_ind == hydro_2_ind
            if ion_binding:
                # helper.easy_print("Ori Seq", sequence)
                # helper.easy_print("Mut Seq", mut_seq)
                # helper.easy_print("Index", mut_ind)
                # helper.easy_print("Ion-binding", mut_seq[mut_ind])
                # print('')
                self._mut_dict[mut_seq] = ("ion-binding", sequence[mut_ind],
                                           mut_seq[mut_ind], mut_ind, disease)
            elif dna_binding:
                # helper.easy_print("Ori Seq", sequence)
                # helper.easy_print("Mut Seq", mut_seq)
                # helper.easy_print("Index", mut_ind)
                # helper.easy_print('DNA-binding', mut_seq[mut_ind])
                # print('')
                self._mut_dict[mut_seq] = ("dna-binding", sequence[mut_ind],
                                           mut_seq[mut_ind], mut_ind, disease)
            elif other_structure:
                self._mut_dict[mut_seq] = ("hydrophobic", sequence[mut_ind],
                                           mut_seq[mut_ind], mut_ind, disease)
            else:
                self._mut_dict[mut_seq] = (notes, sequence[mut_ind],
                                           mut_seq[mut_ind], mut_ind, disease)
            # elif termination:
                # helper.easy_print("Ori Seq", sequence)
                # helper.easy_print("Mut Seq", mut_seq)
                # helper.easy_print("Index", mut_ind)
                # helper.easy_print("Termination", mut_seq[mut_ind])
                # print('')
            #    self._mut_dict[mut_seq] = ("termination", sequence[mut_ind],
            #                               mut_seq[mut_ind], mut_ind)
            # else:
            #    self._mut_dict[mut_seq] = ("none", sequence[mut_ind],
            #                               mut_seq[mut_ind], mut_ind)
            # if not (ion_binding or dna_binding or termination):
            #     helper.easy_print("Mut Seq", mut_seq)
            #     helper.easy_print("Ori Seq", sequence)
            #     helper.easy_print("Index", mut_ind)
            #     helper.easy_print("None", mut_seq[mut_ind])
            #     print('')
        else:
            raise helper.IncorrectError("Domain is not mutated")

    def get_mutation_dictionary(self):
        """
        Returns: the mutation dictionary
        """
        return self._mut_dict

    def set_key_positions(self, full_core_ind, ion_ind, hydro_ind):
        """
        Creates an orderd dictionary of key positions and stores it in the domain.

        Parameter full_core_ind: the index of the core region
        Precondtion: full_core_ind is a list of indices

        Parameter ion_ind: the index of the core region
        Precondtion: fion_ind is a list of indices

        Parameter hydro_ind: the index of the core region
        Precondtion: hydro_ind is a list of indices
        """
        assert type(full_core_ind) == list, "full_core_ind must be a list"
        assert type(ion_ind) == list, "ion_ind must be a list"
        assert type(hydro_ind) == list, "hydro_ind must be a list"
        adict = collections.OrderedDict()
        adict['Core Indices'] = full_core_ind
        adict['Ion-Binding Indices'] = ion_ind
        adict['Hydrophobic Indices'] = hydro_ind
        self._key_positions = adict

    def get_key_positions(self):
        """
        Returns: the positions of the DNA-binding region, ion-binding region, and
        hydrophobic region.
        """
        return self._key_positions

    @staticmethod
    def gene_list_set_domain_positions_sequence(gene_list):
        """
        This function gets domain position for a list of genes.

        This function looks at self._info for each protein in a list of genes.
        The proteins in the gene must have self._info as not blank.

        Parameter gene_list: a list of Genes
        Preconditon: gene_list is a list of genes.
        """
        assert type(gene_list) == list, "Gene list must be a list of genes"
        assert len(gene_list) > 0, "Gene list must be a list of genes"
        assert type(gene_list[0]) == Gene, "Gene list must be a list of genes"
        assert len(gene_list[0].protein_list()[0].info()
                   ) > 0, "Info can not be empty."
        for gene in gene_list:
            proteins = gene.protein_list()
            for protein in proteins:
                protein.set_domain_positions_sequence()

    @staticmethod
    def store_key_positions(gene_list):
        """
        This function stores key position information in a list of genes
        to the domains.

        This function finds the absolute positions for the DNA-binding core, ion-
        binding, and hydrophobic regions of ZF proteins. Primarily used to get
        population information.

        Parameter gene_list: a list of Genes.
        Preconditon: gene_list is a list of Genes
        """
        aminomap = amino_map()
        for gene in gene_list:
            proteins = gene.protein_list()
            for protein in proteins:
                domains = protein.domain_list()
                info_list = protein.info()
                sequence = protein.prosequence()
                for domain in domains:
                    #print('domain.is_valid(True)', domain.is_valid(True))
                    #print("Score", domain.get_score())
                    if domain.is_valid(True):
                        for info in info_list:
                            # print("Test 2")
                            from_ = re.search(
                                r", from (\d{0,4}) to (\d{0,4}): score (\d{1,2}\.\d)", info)
                            dom_num = int(
                                re.search(r": domain (\d{1,2})", info).group(1))
                            if domain.get_domain_number() == dom_num:
                                start_pos = domain.get_start_position()
                                check = re.search(
                                    r"\*->(.+)<-\*", info).group(1)
                                # print(check)
                                c_1_ind = re.search(
                                    r'C\.*p\.*f\.*', check).span()[0]
                                c_2_ind = re.search(
                                    r'C\.*g\.*k\.*', check).span()[0]
                                h_1_ind = re.search(
                                    r'H\.*l\.*r\.*t', check).span()[0]
                                h_2_ind = check.rfind('H')
                                ion_ind = [c_1_ind] + [c_2_ind] + \
                                    [h_1_ind] + [h_2_ind]
                                two_three = check.find('snL')
                                seq = re.search(
                                    "\d{1,4}\s{4}([\w-]+)\s{4}\d{2,4}", info).group(1)
                                # print(seq)
                                core_ind = [two_three - 2, two_three,
                                            two_three + 1, h_1_ind - 1]
                                hydro_ind = [h_1_ind - 9, h_1_ind - 3]
                                num_minus_core = seq[0:-4].count('-')
                                num_minus_c2 = seq[0:c_2_ind].count('-')
                                num_minus_h1 = seq[0:h_1_ind].count('-')
                                num_minus_h2 = seq[0:h_2_ind].count('-')
                                num_minus_hydro1 = seq[0:h_1_ind -
                                                       9].count('-')
                                full_core_ind = [start_pos + i -
                                                 1 - num_minus_core for i in core_ind]
                                full_ion_ind = [start_pos + k -
                                                1 for k in ion_ind]
                                full_ion_ind[0] = start_pos + ion_ind[0] - 1
                                full_ion_ind[1] -= num_minus_c2
                                full_ion_ind[2] -= num_minus_h1
                                full_ion_ind[3] -= num_minus_h2
                                full_hydro_ind = [
                                    start_pos - 1 + h for h in hydro_ind]
                                full_hydro_ind[0] -= num_minus_hydro1
                                full_hydro_ind[1] -= num_minus_h1
                                if [seq[i] for i in core_ind] == [sequence[j]
                                                                  for j in full_core_ind]:
                                    pass
                                else:
                                    # print('Full Core')
                                    # print([sequence[j] for j in full_core_ind])
                                    # print(full_core_ind)
                                    # print(check)
                                    # print(seq)
                                    # print('...')
                                    pass
                                if "CCH" in "".join([seq[i] for i in ion_ind]) and "CCH" in "".join([sequence[i] for i in full_ion_ind]):
                                    pass
                                else:
                                    # print('Ion')
                                    # print([sequence[j] for j in full_ion_ind])
                                    # print(full_ion_ind)
                                    # print([seq[a] for a in ion_ind])
                                    # # print(ion_ind)
                                    # print(check)
                                    # print(seq)
                                    # print('...')
                                    pass
                                if [seq[a] for a in hydro_ind] == [sequence[b] for b in full_hydro_ind]:
                                    pass
                                else:
                                    # print("Hydro")
                                    # print('full', [sequence[b]
                                    #                for b in full_hydro_ind])
                                    # print('partial', [seq[a] for a in hydro_ind])
                                    # print(check)
                                    # print(seq)
                                    # print('...')
                                    pass
                                full_core_ind = [aminomap[sequence[full_core_ind[i]]] + str(full_core_ind[i] + 1)
                                                 for i in range(len(full_core_ind))]
                                full_hydro_ind = [aminomap[sequence[full_hydro_ind[i]]] + str(full_hydro_ind[i] + 1)
                                                  for i in range(len(full_hydro_ind))]
                                full_ion_ind = [aminomap[sequence[full_ion_ind[i]]] + str(full_ion_ind[i] + 1)
                                                for i in range(len(full_ion_ind))]

                                domain.set_key_positions(
                                    full_core_ind, full_ion_ind, full_hydro_ind)


class Gene(object):
    """
    A class representing a gene.

    This class maintains all information about a gene, including the gene's
    common name, it's ID, and the proteins that code for this gene.

    INSTANCE ATTRIBUTES:
        _gene_name:     the common name of the gene. If it doesn't have a common
                        name, the name is the gene_id. [str]
        _gene_id:       the gene ID. [str]
        _protein_list:  the list of proteins this gene codes for [Protein arr]
        _disease_list:  the list of diseases that mutations in this gene causes
                        [str array]
        _direction:     the direction of the transcript
    """

    def __init__(self, gene_name, gene_id='UNKOWNN', protein=''):
        """
        Initializes a gene.

        Parameter gene_name: the name of the gene
        Precondition: gene_name is a valid name for a gene.

        Parameter gene_id: the ID of the gene starting with ENSG
        Precondition: gene_id is a valid ID for a gene.

        Parameter protein: the ID of the protein corresponding to this gene or
        a list of proteins coded by this gene.
        Precondition: protein is a valid ID for a protein or a Protein
        """
        assert type(gene_name) == str, "gene name must be a string"
        assert type(gene_id) == str, "gene ID must be a string"
        assert type(protein) == str or type(protein) == Protein, "protein must \
        be a string or a protein"
        self._gene_name = gene_name
        if type(protein) == str and len(protein) > 0:
            self._protein_list = [Protein(protein)]
        elif type(protein) == Protein:
            self._protein_list = [protein]
        elif type(protein) == str and len(protein) == 0:
            self._protein_list = []
        self._gene_id = gene_id
        self._disease_list = []
        self._direction = 0

    def __eq__(self, other):
        """
        Two genes are equal if they have the same ID.

        Returns: True if this gene is equal to another gene. False otherwise
        """
        return self._gene_id == other._gene_id

    def __lt__(self, other):
        """
        One gene is less than another gene if it hashes to a lower value.
        Used for sorting list of genes.
        """
        return float(self._gene_id[5:]) < float(other._gene_id[5:])

    def __hash__(self):
        """
        Ensures two equivalent genes hash to the same value.
        """
        return hash(self._gene_id)

    def get_gene_name(self):
        """
        Returns: the common name for this gene.

        If this gene doesn't have a common name, it returns the gene ID.
        """
        return self._gene_name

    def get_gene_id(self):
        """
        Returns: the gene ID for this gene.
        """
        return self._gene_id

    def protein_list(self):
        """
        Returns: the list of proteins this gene codes for (arr of Proteins)
        """
        return self._protein_list

    def protein_ids(self):
        """
        Returns: the list of protein IDs this gene codes for (arr of strings)
        """
        prot_list = []
        for protein in self._protein_list:
            prot_list.append(protein.get_protein_id())
        return prot_list

    def protein_str(self):
        """
        Returns: the list of protein IDs this gene codes for as a multi - line string.
        """
        prot_list = self.protein_ids()
        return "\n".join(prot_list)

    def numproteins(self):
        """
        Returns: the number of proteins this gene codes for.
        """
        return len(self._protein_list)

    def disease_list(self):
        """
        Returns: the list of diseases in this gene.
        """
        return self._disease_list

    def disease_list_string(self):
        """
        Returns: the list of diseases mutations in this gene causes
        """
        disease_list = Disease.disease_list_no_rep(self._disease_list)
        length_disease = len(disease_list)
        if length_disease == 0:
            return "No Disease Found"
        elif length_disease == 1:
            return disease_list[0].get_disease_name() + \
                " [from " + disease_list[0].get_source() + "]"
        else:
            sources = Disease.source_list(self._disease_list)
            if len(sources) < 2:
                str_sources = ", ".join(sources)
            if len(sources) >= 2:
                str_sources = ", ".join(sources[0:2])
            if len(sources) > 2:
                str_sources = str_sources + ", etc."
            dis_arr = disease_list[0:2]
            length_arr = len(dis_arr)
            str_disease = ""
            for i in range(length_arr):
                disease = dis_arr[i]
                str_disease = str_disease + disease.get_disease_name()
                if i < length_arr - 1:
                    str_disease = str_disease + ", "
            if length_disease == 2:
                return str_disease + " [from " + str_sources + "]"
            else:
                return str_disease + ", and more [from " + str_sources + "]"

    def disease_str_without_full_info(self):
        """
        Returns: the string of diseases and the source they came from. It doesn't
        count diseases with full info. Used to generate tab - delimited file.
        """
        disease_list = []
        for disease in self._disease_list:
            if not disease.has_full_info():
                disease_list.append(disease)
        length_disease = len(disease_list)
        if length_disease == 0:
            return "No Disease Found (No Source Found)"
        str = ""
        for i in range(length_disease):
            disease = disease_list[i]
            str = str + disease.get_disease_name() + \
                " [" + disease.get_source() + "]"
            if i < length_disease - 1:
                str = str + " | "
        return str

    def get_direction(self):
        """
        Returns: 1 if the direction is forward. -1 if it is reverse
        """
        return self._direction

    def set_direction(self, direction):
        """
        Adds the direction of this protein's transcript
        """
        if direction != 1 or direction != -1:
            assert helper.IncorrectError
        if self._direction == 0 and direction == 0:
            assert helper.IncorrectError
        elif self._direction == 1 and direction == -1:
            assert helper.IncorrectError
        elif self._direction == -1 and direction == 1:
            assert helper.IncorrectError
        self._direction = direction

    def has_disease_w_full_info(self):
        """
        Returns: True if gene has a disease with full disease info. False
        otherwise.
        """
        i = 0
        for disease in self._disease_list:
            if disease.has_full_info():
                i = i + 1
                break
        return i > 0

    def has_disease(self):
        """
        Returns: True if this gene is associated with diseases. False otherwise.
        """
        return len(self._disease_list) > 0

    def set_disease(self, disease, gene_name, source='', mutated=False):
        """
        This procedures adds the string to the list of diseases that mutations
        in this gene can cause.

        Parameter disease: the disease description of the disease
        Preconditon: disease is a string

        Parameter gene_name: the name of the gene this disease affects
        Preconditon: gene_name is a string

        Parameter source: the source the disease came from.
        Preconditon: source is a string

        Parameter mutated: says whether the disease has mutation info.
        Preconditon: mutated is a bool
        """
        assert type(disease) == str or isinstance(disease, Disease), \
            "Disease must be a string or disease"
        assert type(gene_name) == str, "gene name must be a string"
        assert type(source) == str, "Source name must be a string"
        assert type(mutated) == bool, "Mutated must be a bool"
        if type(disease) == str and not mutated:
            disease = Disease(disease, gene_name, source)
        elif type(disease) == str and mutated:
            disease = MutatedDisease(disease, gene_name, source)
        if isinstance(disease, Disease):
            self._disease_list.append(disease)
        else:
            raise helper.IncorrectError("Disease must be a string or disease")

    def set_full_disease_list(self, new_dis_list):
        """
        This procedure sets new_dis_list to be the new disease list for this gene.

        Parameter new_dis_list: a list of Diseases
        Preconditon: new_dis_list is a list of diseases
        """
        self._disease_list = new_dis_list

    def remove_disease(self, disease):
        """
        This procedure removes the disease from the disease list
        """
        ind = self._disease_list.index(disease)
        alist = []
        for i in range(len(self._disease_list)):
            if i != ind:
                alist.append(self._disease_list[i])
        self._disease_list = alist

    def addprotein(self, protein):
        """
        This procedure adds a protein to the list of proteins that this gene
        codes for.
        """
        if isinstance(protein, Protein):
            self._protein_list.append(protein)
        elif isinstance(protein, str):
            self._protein_list.append(Protein(protein))
        else:
            raise ValueError("Protein must be a string or protein")

    def clear_proteins(self, protein=None):
        """
        This procedure clears the protein list. If protein is not none, it will
        set the protein list to be just protein.

        Parameter protein: a Protein to set the gene list to have(or None)
        Preconditon: protein is a Protein or none.
        """
        assert protein is None or type(protein) == Protein, "Protein must " +\
            "be None or a Protein"
        if protein is not None:
            self._protein_list = [protein]
        else:
            self._protein_list = []

    @staticmethod
    def create_gene_list(protein_list, all_proteins):
        """
        This function searches the protein dictionary for each protein, and gets
        the corresponding gene IDs. Then it creates a gene object and adds it to a
        dictionary mapping of genes to proteins

        Parameter protein_list: a list of proteins
        Precondition: protein_list is a list of valid protein IDs

        Parameter all_proteins: a dictionary of proteins that contains information
        including the genes that code for those proteins.
        Precondition: all_proteins is a dictionary of FASTA results of polypeptides

        Returns: a list of genes
        """
        assert type(protein_list) == list or type(protein_list) or set
        assert type(all_proteins) == dict
        # For each protein in the list,
        gene_list = list()
        gene_set = set()
        print("Creating a list of genes.")
        for protein_id in protein_list:
            # Get the description from the corresponding protein
            description = all_proteins[protein_id].description
            # Get the gene ID
            gene_id = re.search("ENSG\d+", description).group(0)
            gene_name = re.search("gene_symbol:(\S+)", description)
            try:
                gene_name = gene_name.group(1)
            except AttributeError:
                gene_name = gene_id
            if gene_name == "":
                gene_name = gene_id
            # print("Gene name: " + gene_name + " | Gene ID: " + gene_id)
            gene = Gene(gene_name, gene_id, protein_id)
            if gene in gene_set:
                gene_index = gene_list.index(gene)
                gene_list[gene_index].addprotein(protein_id)
            else:
                gene_list.append(gene)
                gene_set.add(gene)
        print("List of genes created successfully.")
        return gene_list

    @staticmethod
    def diseased_gene_list(gene_list=[], get_domain_info=True, print_=True):
        """
        This function creates a list of genes that have disease and mutation
        information. If a list is provided, it will create the diseased gene
        list from that list.

        Parameter gene_list: a list of Genes(or an empty list)
        Preconditon: gene_list is a list of Genes or empty

        Parameter get_domain_info: says whether to get the domain info.
        Preconditon: get_domain_info is abool

        Parameter print_: says whether to print to console
        Preconditon: print_ is a bool

        Returns: a list of genes with diseases.
        """
        assert type(gene_list) == list, "Gene list must be a list"
        assert type(get_domain_info) == bool, "Get domain info must be a bool"
        data = hmmer_output()
        protein_dict = Protein.full_protein_dict(print_=print_)
        if gene_list == []:
            try:
                # file = 'test_file_w_disease.tsv'
                gene_list = Gene.gene_list_w_full_dis_info_tab(
                    data, print_=print_)
                Protein.get_proteins(
                    gene_list, data, get_domain_info=get_domain_info, print_=print_)
            except helper.DoesNotExistError:
                print('Error, tab files do not exist. Creating new list of genes.')
                execute_program()
                gene_list = Gene.gene_list_w_full_dis_info_tab(data)
                Protein.get_proteins(
                    gene_list, data, get_domain_info=get_domain_info)
        Protein.gene_list_set_sequence(gene_list, protein_dict)
        Domain.gene_list_set_domain_positions_sequence(gene_list)
        return gene_list

    @staticmethod
    def create_gene_list_tab(file_protein="zf_protein_tabs.tsv",
                             file_disease="zf_genes_w_disease.tsv",
                             nucleotide_info=True, zf_domain_info=True):
        """
        This function searches a tab - seperated file and creates a list of genes
        (with all of the information about proteins and genes stored in each
        gene).

        Parameter file: a file to search for information about the genes.
        Precondtion: file is the name of a file[str]

        Parameter gene_list: the list of genes to build upon(if any)
        Preconditon: gene_list is an array of genes(possible empty)

        Parameter nucleotide_info: a boolean that states whether to get the
        nucleotide info
        Preconditon: nucleotide_info is a boolean

        Parameter nucleotide_info: a boolean that states whether to get the
        zf protein domain info and store it in each domain
        Preconditon: nucleotide_info is a boolean

        Returns: a list of genes
        """
        assert type(file_protein) == str, "Protein file must be a string."
        assert type(file_disease) == str, "Disease file must be a string."
        assert type(
            nucleotide_info) == bool, "Nucleotide info must be a boolean"
        assert type(zf_domain_info) == bool, "Nucleotide info must be a boolean"
        file = helper.edit_filename(file_protein, "sibling", "output")
        if not os.path.exists(file):
            raise helper.DoesNotExistError(
                "Gene info tab with protein info does not exist.")
        data = hmmer_output(create_new=False)
        if nucleotide_info:
            gene_list = Gene.gene_list_w_full_dis_info_tab(data)
        else:
            gene_list = []
        gene_set = set(gene_list)
        print("Scanning " + file_protein + " and building set of genes.")
        with open(file, "r", encoding='utf-8') as f:
            next(f)
            reader = csv.reader(f, dialect="excel", delimiter='\t')
            for row in reader:
                # domain = Domain(pro_id, num, seq, core, valid = True):
                # Create a list of domains
                gene = Gene.create_gene_from_row(row)
                # Add that gene to an array
                if gene not in gene_set:
                    gene_list.append(gene)
                    gene_set.add(gene)
                else:
                    gene_index = gene_list.index(gene)
                    gene_list[gene_index].addprotein(gene.protein_list()[0])
        print("Set of genes created successfully.")
        if zf_domain_info:
            Protein.get_zf_protein_info(data, gene_list, set_domain=False)
        return gene_list

    @staticmethod
    def create_gene_from_row(row):
        """
        This function creates a gene, protein tuple from a row of a tab file.

        This helper function creates a gene and a protein from a row of a tab
        file. It is up to the program that called it to process this gene and
        protein.

        Parameter row: a row to process
        Preconditon: row is a row from csv.reader()
        """
        protein_id = row[2]
        domain_list = []
        core_str = row[3]
        core_list = core_str.split("+")
        num = 1
        for core in core_list:
            domain = Domain(protein_id, num, "UNKNOWN", core)
            domain_list.append(domain)
            num = num + 1
        # Add that list of domains to a protein
        # protein = Protein(id, domains_list)
        protein = Protein(protein_id, domain_list)
        domain_list = []
        # Then create a list of proteins
        # Add that list of proteins to a gene.
        # gene = Gene(name, gene_id, protein)
        gene = Gene(row[0], row[1], protein)
        disease_info = row[5]
        if disease_info != "No Disease Found (No Source Found)":
            disease_info = disease_info.split(" | ")
            for element in disease_info:
                parenthesis_index = element.find("[")
                disease = element[:parenthesis_index - 1]
                source = element[parenthesis_index + 1: -1]
                #helper.easy_print("Element", element)
                gene.set_disease(disease, gene.get_gene_name(), source)
        return gene

    @staticmethod
    def full_gene_id_set(fasta="Homo_sapiens.GRCh38.pep.all.fa"):
        """
        This function takes a fasta file and retrieves all of the genes in it.

        Returns: a set of gene IDs
        """
        assert type(
            fasta) == str, "Fasta must be the name of a fasta file [str]"
        assert "fa" in fasta or "fasta" in fasta, "fasta must be a fasta file."
        fasta = helper.edit_filename(fasta, "sibling", "hmmer")
        with open(fasta) as myfile:
            data = myfile.read()
        return set(re.findall("ENSG\d+", data))

    @staticmethod
    def gene_name_list(gene_list):
        """
        This function takes in a list of genes and returns a list of gene names

        Parameter gene_list: a list of genes to get the names from.
        Precondtion gene_list is a list of Genes.

        Returns: a list of Gene names[str]
        """
        assert type(gene_list) == list or type(gene_list) == set
        new_list = []
        for gene in gene_list:
            new_list.append(gene.get_gene_name())
        return new_list

    @staticmethod
    def gene_id_list(gene_list):
        """
        This function takes in a list of genes and returns a list of gene IDs

        Parameter gene_list: a list of genes to get the IDs from.
        Precondtion gene_list is a list of Genes.

        Returns: a list of Gene IDs[str]
        """
        assert type(gene_list) == list or type(gene_list) == set
        new_list = []
        for gene in gene_list:
            new_list.append(gene.get_gene_id())
        return new_list

    @staticmethod
    def gene_list_w_full_disease_info(gene_list):
        """
        This function looks through a list of genes and finds the genes that
        have disease info. This includes having information on chromosome,
        allelic change, and protein change. It then saves those genes into a new
        list.

        Parameter gene_list: a list of genes
        Preconditon: gene_list must be a list of genes

        Returns: a list of genes with full disease info.
        """
        assert type(gene_list) == list
        assert len(gene_list) > 0
        assert type(gene_list[0]) == Gene
        aset = set()
        print("Creating list of genes with full disease info.")
        for gene in gene_list:
            disease_list = gene.disease_list()
            for disease in disease_list:
                if disease.has_full_info():
                    aset.add(gene)
        alist = list(aset)
        return alist

    @staticmethod
    def gene_list_w_full_dis_info_tab(data, file="zf_genes_w_disease.tsv", print_=True):
        """
        This function searches a tab - seperated file and creates a list of genes
        (with all of the information about proteins and genes stored in each
        gene).

        Parameter file: a file to search for information about the genes.
        Precondtion: file is the name of a file[str]

        Returns: a list of genes
        """
        assert type(file) == str, "File must be a string."
        filename = helper.edit_filename(file, "sibling", "output")
        if not os.path.exists(filename):
            raise helper.DoesNotExistError(
                "Gene info tab with disease info does not exist.")
        gene_list = []
        gene_set = set()
        if print_:
            print("Scanning " + file + " and building set of genes \nwith full " +
                  "mutation info.")
        with open(filename, "r", encoding='utf-8') as f:
            next(f)
            reader = csv.reader(f, dialect="excel", delimiter='\t')
            for row in reader:
                gene_name = row[0]
                gene_id = row[1]
                gene = Gene(gene_name, gene_id)
                disease_name = row[6]
                source = row[8]
                disease = MutatedDisease(disease_name, gene_name, source)
                chromosome = row[2]
                position = literal_eval(row[3])
                nucleotide_change = row[4]
                protein_change = row[5]
                dbSNP = row[7]
                disease.set_chromosome(chromosome)
                disease.set_position(position[0], position[1])
                disease.set_allele_change(nucleotide_change)
                disease.set_amino_change(protein_change)
                disease.set_SNP(dbSNP)
                if gene not in gene_set:
                    gene.set_disease(disease, gene_name, source, True)
                    gene_list.append(gene)
                    gene_set.add(gene)
                else:
                    gene = gene_list[-1]
                    gene.set_disease(disease, gene_name, source, True)
        if print_:
            print("Set of genes with full mutation info created successfully.")
        return gene_list

    @staticmethod
    def known_zf_genes():
        """
        Returns: a set of all known Zinc finger genes that have a score threshold
        below 17.7
        """
        return {'TRPS1'}

    @staticmethod
    def all_gene_stats(fasta="Homo_sapiens.GRCh38.pep.all.fa"):
        """
        This function creates a list of all genes in the human body. It then
        searches multiple databases for disease information about those genes, and
        then prints how many genes are associated with disease. It can also
        output files about those genes.
        """
        import output
        fasta = helper.edit_filename(fasta, "sibling", "hmmer")
        with open(fasta, "r") as f:
            data = f.read()
        full_dict_proteins = Protein.full_protein_dict()
        protein_set = Protein.create_protein_set(data)
        start_time = time.time()
        gene_list = Gene.create_gene_list(protein_set, full_dict_proteins)
        gene_name_list = Gene.gene_name_list(gene_list)
        database.search(gene_list, gene_name_list)
        genes_w_dis = count.genes_w_disease(gene_list)
        len_gene_list = len(gene_list)
        print("")
        print("Stats:")
        print("  " + str(genes_w_dis) + " genes have diseases. ")
        print("  There are " + str(len_gene_list) + " total genes.")
        print("  " + str(round(genes_w_dis / len_gene_list * 100, 2)) +
              "% of genes are associated with disease.")
        output.info(gene_list, full_dict_proteins, create_txt=False,
                    create_tab=False, create_chart=False)
        return gene_list


class Disease(object):
    """
    A class representing a disease.

    This class maintains all information about diseases, including where the
    disease was extracted from.

    INSTANCE ATTRIBUTES:
        _description:       the disease name[str]
        _source:            where the information comes from [str arr]
        _gene_name:         the name of the gene this disease affects
    """

    def __init__(self, description, gene_name, source):
        """
        Initializes a protein.

        Parameter description: a description of the disease
        Preconditon: description is a string

        Parameter domains: a list of protein domains
        Precondition: domains is a list
        """
        assert type(description) == str
        self._description = description
        self._source = source
        self._gene_name = gene_name

    def __eq__(self, other):
        """
        Two diseases are equal if they have the same name.

        Returns: True if this gene is equal to another gene. False otherwise
        """
        if len(self._description) > 10:
            return self._description[0:10].lower() == other._description[0:10].lower()
        else:
            return self._description.lower() == other._description.lower()

    def __hash__(self):
        """
        Ensures two equivalent diseases hash to the same value.
        """
        if len(self._description) > 10:
            return hash(self._description[0:10].lower())
        else:
            return hash(self._description.lower())

    def get_disease_name(self):
        """
        Returns: the disease's name
        """
        return self._description

    def get_gene_name(self):
        """
        Returns: the gene's name that this disease affects
        """
        return self._gene_name

    def get_source(self):
        """
        Returns: the list of sources[arr].
        """
        return self._source

    def has_full_info(self):
        """
        Returns: True if this disease has chromosome, position, and nucleotide
        change information. False otherwise
        """
        return False

    @staticmethod
    def source_list(disease_list):
        """
        Returns: the largest list of sources in a list of disease
        """
        assert type(disease_list) == list or type(disease_list) == set
        arr = []
        for disease in disease_list:
            source = disease.get_source()
            if source not in arr:
                arr.append(source)
        return arr

    @staticmethod
    def disease_list_no_rep(disease_list):
        """
        Returns: the disease list with minimum replicates
        """
        assert type(disease_list) == list or type(disease_list) == set
        arr = []
        for disease in disease_list:
            if disease not in arr:
                arr.append(disease)
        return arr

    @staticmethod
    def delete_repeated_info(gene_list):
        """
        This procedure searches through all of the genes in a gene list and
        removes disease that are present more than once.

        Parameter gene_list: a list of Genes.
        Preconditon: gene_list is a list of genes.
        """
        for gene in gene_list:
            diseases = gene.disease_list()
            dis_set = set(diseases)
            new_dis_list = list(dis_set)
            gene.set_full_disease_list(new_dis_list)

    @staticmethod
    def create_GRCh38_file(gene_list, filename='GRCh38_disease.txt'):
        """
        Creates a text file to output all of the disease nucleotide information.

        This procedure outputs a text file of all of the disease nucleotide info.
        This is done to later search the disease info on population frequency
        databases such as gnomAD or ExAC.

        Parameter gene_list: a list of Genes.
        Preconditon: gene_list is a list of Genes

        Parameter filename: a string to name the output file.
        Preconditon: filename is a non-empty string
        """
        filename = helper.edit_filename(filename, 'sibling', 'output')
        with open(filename, 'w+') as f:
            for gene in gene_list:
                diseases = gene.disease_list()
                for disease in diseases:
                    chromo = disease.get_chromosome()
                    full_pos = disease.get_full_position()
                    pos = full_pos[0]
                    loc = full_pos[1]
                    if loc == 'GRCh38':
                        f.write('chr%s:%s-%s\n' % (chromo, pos, pos))

    @staticmethod
    def convert_to_GRCh37_from_file(gene_list, filename='hglft_genome_2f9c3_f08c20.bed'):
        """
        Converts the disease information in a list of genes to GRCh37.

        This procedure converts all GRCh38 positions to GRCH37. The file has to be
        manually created using liftover. After creating a file from liftover,
        the filename has to be changes to match the outputted filename.

        Website: https://genome.ucsc.edu/cgi-bin/hgLiftOver

        Parameter gene_list: a list of Genes.
        Preconditon: gene_list is a list of Genes

        Parameter filename: a string to name the output file.
        Preconditon: filename is a non-empty string
        """
        filename = helper.edit_filename(filename, 'sibling', 'databases')
        with open(filename, 'r') as f:
            for gene in gene_list:
                diseases = gene.disease_list()
                for disease in diseases:
                    chromo = disease.get_chromosome()
                    full_pos = disease.get_full_position()
                    loc = full_pos[1]
                    if loc == 'GRCh38':
                        line = f.readline()
                        match = re.search("chr(\w+):(\d+)-(\d+)", line)
                        try:
                            line_chromo = int(match.group(1))
                        except ValueError:
                            line_chromo = match.group(1)
                        line_pos = match.group(2)
                        if chromo == line_chromo:
                            disease.set_position(line_pos, 'GRCh37')
                        else:
                            print(disease.get_disease_name())
                            print(chromo)
                            print(match)
                            raise helper.IncorrectError("The output file from liftover" +
                                                        " in the same order as the list of Genes.")


class MutatedDisease(Disease):
    """
    A class representing a disease.

    This class maintains all information about diseases, including where the
    disease was extracted from.

    INSTANCE ATTRIBUTES:
        _description:       the disease name[str]
        _source:            where the information comes from [str arr]
        _dbSNP:             the SNP number to get information about the disease[str]
        _allele_change:     the allele change associated with this disease[str]
        _amino_change:      the amino acid change associated with this disease[str]
        _chromosome:        the chromosome this gene is in [int or str]
        _position:          the posiiton of this gene[int]
        _hit_domain:        whether or not this disease hit a domain after
                            running mutated_gene_list[bool]
        _gene_name:         the name of the gene this disease affects
    """

    def __init__(self, description, gene_name, source):
        """
        Initializes a protein.

        Parameter description: a description of the disease
        Preconditon: description is a string

        Parameter domains: a list of protein domains
        Precondition: domains is a list
        """
        super().__init__(description, gene_name, source)
        self._dbSNP = '-'
        self._allele_change = ''
        self._amino_change = ''
        self._chromosome = 0
        self._position = 0
        self._hit_domain = False

    def __eq__(self, other):
        """
        Two diseases are equal if they have the same name.

        Returns: True if this gene is equal to another gene. False otherwise
        """
        if type(other) == MutatedDisease:
            return self._amino_change == other._amino_change \
                and self._chromosome == other._chromosome

    def __hash__(self):
        """
        Ensures two equivalent diseases hash to the same value.
        """
        return hash((self._amino_change, str(self._chromosome), self._description[0:4].lower()))

    def get_disease_name(self):
        """
        Returns: the disease's name
        """
        return self._description

    def get_gene_name(self):
        """
        Returns: the gene's name that this disease affects
        """
        return self._gene_name

    def get_source(self):
        """
        Returns: the list of sources[arr].
        """
        return self._source

    def get_dbSNP(self):
        """
        Returns: the dbSNP for this disease
        """
        return self._dbSNP

    def get_allele_change(self):
        """
        Returns: the allele change for this disease
        """
        return self._allele_change

    def get_amino_change(self):
        """
        Returns: the amino acid change for this disease
        """
        return self._amino_change

    def get_chromosome(self):
        """
        Returns: the chromosome for this gene
        """
        try:
            return int(self._chromosome)
        except ValueError:
            return self._chromosome

    def get_position(self):
        """
        Returns: the posiiton for this gene
        """
        return int(self._position[0])

    def get_string_full_position(self):
        """
        Returns: the full posiiton for this gene (so its position and position type)
        as a string.
        """
        return str(self._position)

    def get_full_position(self):
        """
        Returns: the full posiiton for this gene (so its position and position type)
        as a tuple
        """
        return self._position

    def has_full_info(self):
        """
        Returns: True if this disease has chromosome, position, and nucleotide
        change information. False otherwise
        """
        if self._amino_change != '' and self._allele_change != '' and \
                self._chromosome != 0:
            try:
                int(self._amino_change[-1])
                return False
            except ValueError:
                return True
        return False

    def get_SNP_info(self):
        """
        This procedure searches NLM for info about SNPs associated with
        this disease.

        Return: True if there was SNP info to get. False otherwise
        """
        if self._dbSNP != '-':
            url1 = "https://www.ncbi.nlm.nih.gov/snp/?term="
            ua = UserAgent()
            headers = {"User-Agent": ua.random}
            request = Request(url1 + self._dbSNP, headers=headers)
            response = urlopen(request)
            respData = response.read()
            response.close()
            page_soup = soup(respData, "html.parser")
            try:
                info = page_soup.findAll(
                    "dl", {"class": "snpsum_dl_left_align"})[0].text
                allele_change = re.findall(r"Alleles:(\w{1,2}>\S+)", info)[0]
                chromo, pos = re.findall(
                    r"Chromosome: (\d{1,2}|X|Y):(\d+)", info)[0]
                self._chromosome = chromo
                self._position = (int(pos), 'GRCh38')
                self._allele_change = allele_change
            except:
                print(self._dbSNP + " was not recognized.")

    def get_hit_domain(self):
        """
        Returns: True if this disease affected a domain. False otherwise.
        """
        return self._hit_domain

    def set_SNP(self, dbSNP):
        """
        This procedure adds a dbSNP to this disease.
        """
        self._dbSNP = dbSNP

    def set_amino_change(self, AA):
        """
        This procedure adds the amino acid change associated with this disease.
        """
        self._amino_change = AA

    def set_allele_change(self, DNA):
        """
        This procedure adds the amino acid change associated with this disease.
        """
        self._allele_change = DNA

    def set_chromosome(self, chromo):
        """
        This procedure adds the chromosome of the gene this disease is in.
        """
        self._chromosome = chromo

    def set_position(self, pos, type_pos):
        """
        This procedure adds the position of this mutation to the disease.

        Parameter pos: the position to set.
        Preconditon: pos is an int

        Parameter type_pos: whether the position is relative, chromosomal, or unknown
        Preconditon: type_pos is relative, chromosomal, or unknown
        """
        assert type(
            type_pos) == str, "Second element in tuple must be a string."
        pos = int(pos)
        self._position = (pos, type_pos)

    def set_hit_domain(self, abool=True):
        """
        This procedure sets hit domain to be equal to abool. A domain is hit
        when the mutation is located within one of the domains within a protein.

        Parameter abool: a boolean to set hit domain to.
        Preconditon: abool is a bool.
        """
        assert type(abool) == bool, "abool must be a bool."
        self._hit_domain = abool
