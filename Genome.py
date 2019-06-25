class Genome:

    """
    A class that represents a singular genome parsed from the Genbank file.

    organism_name: name of the organism
    genbank_genome: the entire genome from the genbank file

    """
    def __init__(self, genbank_genome: list, refseq_ids: list) -> None:
        self.organism_name = ""
        self.genbank_genome = genbank_genome
        self.refseq_ids = refseq_ids

    def get_organism_name(self):
        for line in self.genbank_genome:
            if "ORGANISM" in line:
                self.organism_name = line.replace("  ORGANISM", "")

    def print_name_and_genome(self):
        print(f"*****{self.organism_name} Genome*****")
        print(self.genbank_genome)

    def find_locus_tag(self) -> list:
        for item in self.refseq_ids:
            for i in range(len(self.genbank_genome)):
                if item in self.genbank_genome[i]:
                    return self._search_for_locus_tag()

    def _search_for_locus_tag(self) -> list:
        loci = []
        j = 0
        while j <= 5:
            j += 1
            if "old_locus_tag" in self.genbank_genome[j - 1]:
                quotes = []
                for i, char1 in enumerate(self.genbank_genome[i - 1]):
                    if char1 == '"':
                        quotes.append(j)
                    locus = self.genbank_genome[i - 1][quotes[0] + 1:quotes[1]]
                    loci.append(locus)
        return loci

