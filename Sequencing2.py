from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from Genome import Genome
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
print("****** ids ********")
with open("my_genbank.gb", "w") as out_handle_2:
    for _id in ids:
        print(_id, flush=True)
        result_handle_2 = Entrez.efetch(db="nucleotide", id=_id, rettype="gb",
                                      retmode="text")
        out_handle_2.write(result_handle_2.read())
        result_handle_2.close()
        print(result_handle_2, flush=True)

all_genomes = []
with open("my_genbank.gb", "rU") as handle:
    for record in SeqIO.parse(handle, "genbank"):
        genome = Genome(record, refseqs)
        all_genomes.append(genome)

for genome in all_genomes:
    genome.get_organism_name()
    print("************LOCUS TAG!!!!!**************")
    print(genome.find_locus_tag())



