"""
Convert a GFF file and a genome file into an annotation FASTA file.
"""
from sys import argv

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.

    Args:
        seq (str): The DNA sequence.

    Returns:
        str: The reverse complement of the input sequence.
    """
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return "".join([complement[i.upper()] if i.upper() == i else complement[i.upper()].lower() for i in seq][::-1])

def gff_genome_to_annfa(gff_file, genome_file, output_file):
    """
    Convert a GFF file and a genome file into an annotation FASTA file.

    Args:
        gff_file (str): The path to the GFF file.
        genome_file (str): The path to the genome file.
        output_file (str): The path to the output annotation FASTA file.

    Returns:
        None
    """
    with open(gff_file, 'r', encoding='utf-8') as gff:
        with open(genome_file, encoding='utf-8') as genome:
            seq = "".join(genome.readlines()).strip().split('>')[1:]
        seq = {i.split(None,1)[0]:"".join(i.split()[1:]) for i in seq}
        with open(output_file, 'w', encoding='utf-8') as output:
            for line in gff:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip().split('\t')
                    print( f">{''.join(line[8:])}_{line[0]}_{line[2]}", file=output )
                    if line[6] == '+':
                        print( seq[line[0]][int(line[3])-1:int(line[4])], file=output )
                    else:
                        s = seq[line[0]][int(line[3])-1:int(line[4])]
                        print( reverse_complement(s)
                              , file=output )

    return

if __name__ == '__main__':
    print(argv)
    if len(argv) == 4:
        gff_genome_to_annfa(argv[1], argv[2], argv[3])
    elif len(argv) == 3:
        gff_genome_to_annfa(argv[1], argv[2], 'output.fa')
    else:
        print('Usage:\npython gff_genome_to_annfa.py <gff_file> <genome_file> [output_file]')
