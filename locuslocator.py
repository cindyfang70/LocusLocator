from Sequencing2 import *

print("Hello and welcome to the Locus Locator!")
print("Would you like to (1) upload your own genomes or (2) have "
      "me fetch the genomes for you?")
user_choice = input("Please enter 1 for option (1) and 2 for option (2)")

if user_choice == "2":
    filename = "my_genbank.gb"
    protein_sequence = input("Please enter the FASTA sequence of the protein "
                             "you would like to locate.")
    print("Parsing...", flush=True)
    ids = qblast(protein_sequence)
    refseqs = blastp(protein_sequence)
    get_genomes(ids)
    results = find_refseq_matches_in_genome(refseqs, filename)
    print(results)
    for result in results:
        print("upstream locus tags for", result[0], ":")
        print(find_protein_products(find_upstream(convert_to_int(result[1]))))
        print("downstream locus tags for", result[0], ":")
        print(find_protein_products(find_downstream(convert_to_int(result[1]))))
elif user_choice == "1":
    user_refseq = input("Please enter the refseq you would like to locate:")
    print("Before we begin, please ensure that the genome files you would like "
          "me to search are in the same folder as this program!")
    file = input("Please enter the file name of the genome in which you "
                  "would like to look for the refseq:")
    all_files = []
    all_files.append(file)
    while file != "f":
        file = input("Please enter the name of the next file you would like me "
                     "to search. If there are no more files to search through, "
                     "please enter f.")
        if file != "f":
            all_files.append(file)
    refseqs = [user_refseq]
    results = []
    for file in all_files:
        result = find_refseq_matches_in_genome(refseqs, file)
        results.append(result)
    print(results)
    for result in results:
        print(result)
        print(result[0][1])
        print("upstream locus tags for", result[0][0], ":")
        print(find_protein_products(find_upstream(convert_to_int(result[0][1]))))
        print("downstream locus tags for", result[0][0], ":")
        print(find_protein_products(find_downstream(convert_to_int(result[0][1]))))
