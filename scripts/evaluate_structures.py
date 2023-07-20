import argparse
from rna_folding.evaluate import f1_score


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotypes", required=True, type=str, help="List of genotypes")
    parser.add_argument("-p", "--phenotypes", required=True, type=str, help="List of phenotypes and genotype IDs that map to them")
    parser.add_argument("-r", "--reference", required=True, type=str, help="Reference RNA secondary structure in dot-bracket notation")
    parser.add_argument("-a", "--abstract", required=False, type=int, help="Level of abstraction on which to compare structures")

    args = parser.parse_args()
    
    seqid_to_db = {}
    if args.file:
        # Parse the phenotype output to be a sequence ID to RNA SS dict
        with open(args.phenotypes, "r") as query_file:
            for line in query_file:
                l = line.split()
                db = l[0]
                for i in l[1:]:
                    id = int(i)
                    if id in seqid_to_db:
                        seqid_to_db[id].append(db)
                    else:
                        seqid_to_db[id] = [db]
        
        # Parse the genotype input file (line number (-1) = sequence ID
        query_seq_list = []
        with open(args.genotypes, "r") as ref_file:    
            for id, line in enumerate(ref_file):
                query_seq_list.append(line.strip())
        print(query_seq_list)

        
