from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
from typing import List
class Genome:
    """
    A class that represents a singular genome parsed from the Genbank file.

    organism_name: name of the organism
    genbank_genome: the entire genome from the genbank file

    """
    def __init__(self, organism_name: str, genbank_genome: list) -> None:
        self.organism_name = organism_name
        self.genbank_genome = genbank_genome

    def print_name_and_genome(self):
        print(f"*****{self.organism_name} Genome*****")
        print(self.genbank_genome)

    def find_locus_tag(self) -> list:
        for item in refseqs:
            for i in range(len(self.genbank_genome)):
                if item in line:
                    return self._search_for_locus_tag()

    def _search_for_locus_tag(self) -> list:
        loci = []
        j = 0
        while j <= 5:
            j += 1
            if "old_locus_tag" in self.genbank_genome[i - 1]:
                quotes = []
                for j, char1 in enumerate(self.genbank_genome[i - 1]):
                    if char1 == '"':
                        quotes.append(j)
                    locus = self.genbank_genome[i - 1][quotes[0] + 1:quotes[1]]
                    loci.append(locus)
        return loci




protein_sequence = input("Please enter the FASTA of the protein you would like "
                         "to BLAST: ")
print("Parsing....", flush=True)

""" use blastp to find similar proteins to the query and return their sequence id
and RefSeq ID
"""
result_handle = NCBIWWW.qblast("tblastn", "nr", protein_sequence)
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
result_handle = open("my_blast.xml")
blast_record = NCBIXML.parse(result_handle)
print("***********Results************")
f = open('my_blast.xml', 'r')
ids = []
for line in f:
    if "Hit_id>gi" in line:
        bars = []
        for i, char in enumerate(line):
            if char == '|':
                bars.append(i)
        _id = line[bars[2]+1:bars[3]]
        ids.append(_id)

# blastp the query sequence to find refseq ids of similar proteins
result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence)
with open("my_blastp.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
result_handle = open("my_blastp.xml")
blast_record = NCBIXML.parse(result_handle)
f = open("my_blastp.xml", "r")
refseqs = []
print("*******RefSeqs********")
for line in f:
    if "Hit_id" in line and "ref" in line:
        bars = []
        for i, char in enumerate(line):
            if char == "|":
                bars.append(i)
        refseq = line[bars[0]+1:bars[1]]
        print(refseq)
        refseqs.append(refseq)

# Get the genbank files of all the similar proteins


Entrez.email = "cindyfang70@gmail.com"
with open("my_genbank.xml", "w") as out_handle_2:
    for _id in ids:
        print(_id)
        result_handle_2 = Entrez.efetch(db="nucleotide", id=_id, rettype="gb",
                                      retmode="text")
        out_handle_2.write(result_handle_2.read())
        result_handle_2.close()
        print(result_handle_2)

print("*******organism names + loci*******")

with open("my_genbank.xml", "r") as genbank:
    lines = genbank.read().splitlines()
organism_names = []
temp_locus = ""
genomes = []
genome = []
for i in range(len(lines)):
    if "/locus" in lines[i]:
        temp_locus = lines[i][11:-1]

    if "/organism" in lines[i]:
        line = lines[i].strip()[11:-1]
        if line not in organism_names:
            organism_names.append(line)
        if "ACCESSION" in lines[i]:
            genome.append(lines[i-1:i+1])
            i += 1
        while i < len(lines) and "ACCESSION" not in lines[i]:
            genome.append(lines[i])
            i += 1
        genomes.append(genome)
        genome.append(lines[i])

result_genomes = []
for i in range(len(genomes)):
    new_genome = Genome(organism_names[i], genomes[i])
    result_genomes.append(new_genome)

print(result_genomes)
for genome in result_genomes:
    genome.print_name_and_genome()
    print(genome.find_locus_tag())



""" 
search for refseq id instead of sequence string 
use blastp on query to get refseq ids 
id #s increment by 5 at a time 

what if the blastp has a gb id instead of refseq? 
how to find corresponding locus tag of the refseq id in the genbank? 


-get the temp locus, see if the refseqid it refers to matches any of the ones
from blastp 
if it does, then return the locus 
use the returned loci to find the proteins upstream and downstream
fetch sequences 15 upstream and 15 downstream
eutils ncbi package to fetch sequences 
aligh locus tag w organism
fetch accession number of genome from genbank
query the genbank with the refseq id, then look up to 5 lines above for the 
locus tag
take old locus tag instead of locus tag
use locus tag to find genome using eutilities 


"""


