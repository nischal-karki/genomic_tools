"""This script extracts genes from a genome based on a search term in a GFF file and outputs them in GenBank format.
    Usage:
        python select_genes.py <genome_file> <gff_file> <search_term> [offset=0] [output=standard output] [-fa] [-h] [--help]
        Arguments:
            genome_file: The genome file in fasta format.
            gff_file: gff file corresponding to the genome file.
            search_term: The search term to use.
            offset: Amount extra around the sequence to use when extracting the sequences. Default is 0.
            output: The output file to write the gb text to. Default is to print the output.
        Flags:
            -fa, --fasta: Output the fasta sequences after the gb output.
            -h, --help: Print this help message.
"""

import os
import sys
if os.path.abspath("../Useful Python/") not in sys.path:
    sys.path.insert(0, os.path.abspath("../Useful Python/"))
from multiprocessing import Pool
from extract_seq_from_genome_based_on_gff_feature import *
from format_genbank import makegb
from gen_seq_from_motifs import *

def select_genes(genome_file: str, gff_file: str, search_term: str, offset: int = 0, output: str | typing.TextIO = sys.stdout):
    found = extract_seq_from_genome_based_on_gff_feature(
        genome_file, gff_file, search_term, offset
        )
    fasta_seq = ""
    for gene, (seq, start, features) in found.items():
        features = reformat_features_for_gb(seq, start, features)
        
        if isinstance(output,str):
            real_out = open(output + gene + ".gb", "w")
        else:
            real_out = output
            print(gene,file=output)
        print(
            makegb(
                gene= gene,
                seq= seq,
                features= features,
            ),
            file=real_out,
            flush=True
        )
        if isinstance(output,str):
            real_out.close()
        fasta_seq += f">{gene}\n{seq}\n"
    return fasta_seq

if __name__ == "__main__":
    import sys
    fasta = False
    if "-fa" in sys.argv:
        sys.argv.remove("-fa")
        fasta = True
    if "--fasta" in sys.argv:
        sys.argv.remove("--fasta")
        fasta = True
    if "-h" in sys.argv or "--help" in sys.argv or len(sys.argv) < 4:
        print( __doc__ )
        sys.exit()
    args = {i:j for i,j in enumerate(sys.argv)}
    genome_file = args.get(1)
    gff_file = args.get(2)
    search_term = args.get(3)
    offset = int(args.get(4,0))
    output = args.get(5,sys.stdout)
    fasta_seq = select_genes(genome_file, gff_file, search_term, offset, output)
    if fasta:
        print(fasta_seq)
