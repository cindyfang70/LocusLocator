from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

class Genome:
    """
    A class that represents a singular genome parsed from the Genbank file.

    organism_name: name of the organism
    genbank_genome: the entire genome from the genbank file

    """
    def __init__(self, genbank_genome: list) -> None:
        self.organism_name = ""
        self.genbank_genome = genbank_genome

    def get_organism_name(self):
        for line in self.genbank_genome:
            if "ORGANISM" in line:
                self.organism_name = line.replace("  ORGANISM", "")

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
with open("./data/my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
result_handle = open("./data/my_blast.xml")
blast_record = NCBIXML.parse(result_handle)
print("***********Results************", flush=True)
f = open('./data/my_blast.xml', 'r')
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
with open("./data/my_blastp.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
result_handle = open("./data/my_blastp.xml")
blast_record = NCBIXML.parse(result_handle)
f = open("./data/my_blastp.xml", "r")
refseqs = []
print("*******RefSeqs********", flush=True)
for line in f:
    if "Hit_id" in line and "ref" in line:
        bars = []
        for i, char in enumerate(line):
            if char == "|":
                bars.append(i)
        refseq = line[bars[0]+1:bars[1]]
        print(refseq, flush=True)
        refseqs.append(refseq)

# Get the genbank files of all the similar proteins


Entrez.email = "cindyfang70@gmail.com"
with open("my_genbank.xml", "w") as out_handle_2:
    for _id in ids:
        print(_id, flush=True)
        result_handle_2 = Entrez.efetch(db="nucleotide", id=_id, rettype="gb",
                                      retmode="text")
        out_handle_2.write(result_handle_2.read())
        result_handle_2.close()
        print(result_handle_2, flush=True)

with open("my_genbank.xml", "r") as genbank:
    lines = genbank.read().splitlines()
organism_names = []
temp_locus = ""
genomes = []

i = 0
while i < len(lines):
    genome = []
    if "/locus" in lines[i]:
        temp_locus = lines[i][11:-1]
    i += 1

j = 0
while j < len(lines):
    if "DEFINITION" in lines[j] and "complete genome" in lines[j]:
        print("this should be the definition of the genome", flush=True)
        print(lines[j-1:j+2], flush=True)
        genome.append(lines[j-1:j+2])
        j += 2
        while j < len(lines) and "ACCESSION" not in lines[j]:
            genome.append(lines[j])
            j += 1
        new_genome = Genome(genome)
        new_genome.get_organism_name()
        genomes.append(new_genome)

for genome in genomes:
    genome.print_name_and_genome()
    print(genome.find_locus_tag(), flush=True)



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


