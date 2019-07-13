from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
import xml.etree.ElementTree as ET


protein_sequence = input("Please enter the FASTA of the protein you would like "
                         "to BLAST: ")
print("Parsing....", flush=True)

""" use blastp to find similar proteins to the query and return their sequence id
and RefSeq ID
"""
result_handle = NCBIWWW.qblast("tblastn", "nr", protein_sequence)
with open("./data/my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
result_handle = open("./data/my_blast.xml")
blast_record = NCBIXML.parse(result_handle)
print("***********Results************", flush=True)

ids = []
tree = ET.parse('./data/my_blast.xml')
root = tree.getroot()
for child in root.iter():
    if child.tag == "Hit_id":
        bars = []
        for i in range(len(child.text)):
            if child.text[i] == '|':
                bars.append(i)
        _id = child.text[bars[2]+1:bars[3]]
        ids.append(_id)
print(ids)
# blastp the query sequence to find refseq ids of similar proteins
result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
with open("./data/my_blastp.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
result_handle = open("./data/my_blastp.xml")
blast_record = NCBIXML.parse(result_handle)
f = open("./data/my_blastp.xml", "r")
refseqs = []
tree = ET.parse("./data/my_blastp.xml")
root = tree.getroot()
for child in root.iter():
    if child.tag == "Hit_id" and "ref" in child.text:
        refseqs.append(child.text.strip("ref").strip("|"))
print(refseqs)

# Get the genbank files of all the similar proteins
Entrez.email = "cindyfang70@gmail.com"
with open("my_genbank.gb", "w") as out_handle_2:
    for _id in ids:
        result_handle_2 = Entrez.efetch(db="nucleotide", id=_id, rettype="gb",
                                      retmode="text")
        out_handle_2.write(result_handle_2.read())
        result_handle_2.close()
        print(result_handle_2, flush=True)

<<<<<<< HEAD
# refseqs =  ['WP_028357638.1', 'WP_032962436.1', 'WP_126980089.1', 'WP_126708169.1', 'WP_081462204.1', 'WP_126705769.1', 'WP_106448375.1', 'WP_013516228.1', 'WP_068832127.1', 'WP_049360089.1', 'WP_068832973.1', 'WP_014743597.1', 'WP_051629968.1', 'WP_044587254.1', 'XP_001301229.1', 'WP_036972373.1']
print("*****Matches:*****")
results = []
=======
all_genomes = []
wanted_loci = {}
>>>>>>> 9ce45d87ce8f2c49101948a593d8e5cb16e49b8e
with open("my_genbank.gb", "rU") as handle:
    for record in SeqIO.parse(handle, "genbank"):
        # genome = Genome(record, refseqs)
        # all_genomes.append(genome)
        for feature in record.features:
<<<<<<< HEAD
            if feature.type == "source":
                if "organism" in feature.qualifiers:
                    temp_organism_name = feature.qualifiers["organism"][0]
            if feature.type == "CDS":
                if "locus_tag" in feature.qualifiers:
                    temp_locus = feature.qualifiers["locus_tag"][0]
                if "inference" in feature.qualifiers:
                    for refseq in refseqs:
                        if refseq in feature.qualifiers["inference"][0]:
                            result = (temp_organism_name, temp_locus,
                                      refseq)
                            results.append(result)
                            print(result)

=======
            if feature.type =="CDS":
                for line in feature.type:
                    ßprint(line)


for genome in all_genomes:
    genome.get_organism_name()
    print("************LOCUS TAG**************")
    print(genome.find_locus_tag())
>>>>>>> 9ce45d87ce8f2c49101948a593d8e5cb16e49b8e

def convert_to_int(locus_tag: str) -> tuple:
    numbers_in_id = ""
    chars_before_numbers = ""
    for i in range(len(locus_tag)):
        chars_before_numbers = chars_before_numbers + locus_tag[i]
        if locus_tag[i] == "_":
            underscore_index = i
            break
    for i in range(underscore_index + 1, len(locus_tag)):
        if locus_tag[i].isnumeric() and i > underscore_index:
            numbers_in_id = numbers_in_id + locus_tag[i]
    return int(numbers_in_id), chars_before_numbers


<<<<<<< HEAD
def find_upstream(locus_tup: tuple) -> list:
    i = 0
    new_locus_numbers = []
    locus_numbers = locus_tup[0]
    while i < 10:
        i += 1
        locus_numbers += 5
        new_locus = str(locus_numbers)
        new_locus = locus_tup[1] + new_locus
        new_locus_numbers.append(new_locus)
    return new_locus_numbers


def find_downstream(locus_tup: tuple) -> list:
    i = 0
    new_locus_numbers = []
    locus_numbers = locus_tup[0]
    while i < 10:
        i += 1
        locus_numbers -= 5
        new_locus = str(locus_numbers)
        new_locus = locus_tup[1] + new_locus
        new_locus_numbers.append(new_locus)
    return new_locus_numbers


def find_protein_products(locus_tags: list) -> list:
    protein_products = []
    with open("my_genbank.gb", "rU") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    for locus in locus_tags:
                        if feature.qualifiers["locus_tag"][0] == locus:
                            protein_products.append((locus,
                                                     feature.qualifiers
                                                     ["product"][0]))

    return protein_products


for result in results:
    print("upstream locus tags:")
    print(find_protein_products(find_upstream(convert_to_int(result[1]))))
    print("downstream locus tags")
    print(find_protein_products(find_downstream(convert_to_int(result[1]))))







"""
how to identify if the protein is in a morphogenetic region?

find locus tags of refseq id and +/-15 without blasting (given genome)
read product of locus tags 

flowchart

lambda phage 
hk97 phage
mu phage 
p22 phage 
p2 phage 
spp1 phage 

read genomic arrangement for phages
=======
"""
locus tags don't match up with refseq ids (refseq ids do not appear in the
genomes of the organisms)

difference between tblastn and blastp 

how are the results of tblastn and blastp related 

get only locus tags related to the refseq ids 

>>>>>>> 9ce45d87ce8f2c49101948a593d8e5cb16e49b8e

"""
