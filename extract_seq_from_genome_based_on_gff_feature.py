import sys
import typing
import os
import io
from format_genbank import makegb
from smith_waterman import smith_waterman
def translate(seq: str) -> str:
    return "".join(codon_table[seq[i:i+3]] for i in range(0,len(seq),3))
def align(seq1: str, seq2: str) -> str:
    return smith_waterman(seq1, seq2)[0]

def align_translated_product_to_genomic_DNA(protein: str, dna: str, **kwargs) -> str:
    """
    Aligns a protein to a DNA sequence by translating the DNA sequence in all six frames.
    The return value is the 3 frame translation in the strand that the protein aligns to.
    Args:
        protein: The protein sequence to align
        dna: The DNA sequence to align
        **kwargs: Is used internally by the function and will be ignored.
    Returns:
        Strand (+ or -)
        List of the 3 frame alignments
    """
    
    alignments = [align(protein, translate(dna[i:])) for i in range(3)]
    
    formatted_alns = []
    for aln in alignments:
        formatted_alns.append( "" )
        tar = aln.split("-")
        for j in tar:
            if j == "":
                continue
            loc = protein.index( j )
            prefix = "-" * ( loc - len(formatted_alns[-1]) )
            formatted_alns[-1] += prefix + j
        formatted_alns[-1] += "-" * ( len(protein) - len(formatted_alns[-1]) )

    recreated_protein = ""
    for i in range( len(protein) ):
        for j in formatted_alns:
            if j[i] == "-":
                continue
            recreated_protein += j[i]
            break
    test_frames = recreated_protein == protein

    if test_frames or kwargs.get("return_anyways", False):
        return ("+", alignments)
    elif kwargs.get("tried_once", True):
        return ("-", align_translated_product_to_genomic_DNA(protein, reverse_complement(dna), tried_once=False)[1])
    else:
        raise ValueError("Could not find the protein in the DNA sequence.")
        
    
def identifyCDS(target_protein, query_DNA):
    """
    Identifies the CDS of the query DNA based on the target protein and returns as features.
    Args:
        target_protein: The protein sequence to search
        query_DNA: The DNA sequence to search
    Returns:
        A list of dictionaries of the CDS features.
        {
            "type":"CDS",
            "location":(start, end),
            "strand":strand,
            "qualifiers": {
                "translation":translation_product
            }
        }
    """
    try:
        strand, alignments = align_translated_product_to_genomic_DNA(target_protein, query_DNA)
    except ValueError:
        raise ValueError(f"Could not find the protein {target_protein} in the DNA sequence.")
    
    if strand == "-":
        query_DNA = reverse_complement(query_DNA)
    
    frags = []
    for k, aln in enumerate(alignments):
        n_dna = query_DNA[k:]
        trans_prot = translate( n_dna )
        for frag in aln.split("-"):
            if frag == "" or len(frag) < 5:
                continue
            try:
                loc = trans_prot.index(frag) * 3 + k
            except ValueError:
                continue
            frags.append( (loc, frag) )
        

    cds_features = []
    for loc, frag in frags:
        if strand == "+":
            start =  loc + 1
            end = len(frag) * 3 + start
            location = f"{start}..{end}"
        elif strand == "-":
            end = len(query_DNA) - loc
            start = end - len(frag) * 3 + 1
            location = f"complement({start}..{end})"
        cds_features.append(
            {
                "type":"CDS",
                "location": location,
                "strand":strand,
                "translation":frag
            }
        )
    print( "Note that the annotations might not be correct, since splicing is not considered." )
    return cds_features


# Codon table
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

codon_table = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
codon_table.update( {x.replace("U","T"):y for x,y in codon_table.items()} )
codon_table.update( {x.lower():y.lower() for x,y in codon_table.items()} )


# SOMEONE CHECK IF THIS IS CORRECT, I TRUSTED GITHUB COPILOT FOR THIS.
complement = {
    "A":"T",
    "C":"G", 
    "N":"N", 
    "W":"W", 
    "S":"S", 
    "K":"M", 
    "R":"Y", 
    "B":"V", 
    "D":"H",
}
complement.update({y:x for x,y in complement.items()})
complement.update({x.lower():y.lower() for x,y in complement.items()})

def find_features_between(
        gff_file: str|os.PathLike|typing.IO,
        chromosome: str,
        start: int,
        end: int
        ) -> list[str]:
    """
    Extracts the features from a gff file between the start and end.
    Args:
        gff_file: The gff file to search
        chromosome: The chromosome to search
        start: The start of the search
        end: The end of the search
    Returns:
        A list of the features found.
    """
    features = []
    if isinstance(gff_file,list):
        f = gff_file
    else:
        f = open(gff_file, "r") if isinstance(gff_file,str) else gff_file
        if "b" in f.mode:
            f = io.TextIOWrapper(f)

    for line in f:
        vals = line.split("\t")
        if vals[0] != chromosome:
            continue
        vals = [*map(int,vals[3:5])]
        compare = [start,end]
        if max([*compare,*vals]) in compare:
            if min([*compare,*vals]) in compare:
                features.append(line)
    if isinstance(f, typing.IO):
        f.close()
    return features

def find_gene_loc(
        gff_file: str|os.PathLike|typing.IO|list[str], 
        search_terms: str|list[str]
        ) -> dict[str,list[str]]:
    """
    Extracts the gene from a gff file based on the search terms
    Args:
        gff_file: The gff file to search
        search_terms: The search terms to use
    Returns:
        A dictionary of the search terms and the features found.
    """
    
    if isinstance(search_terms,str):
        search_terms = [search_terms]
    
    if isinstance(gff_file,list):
        gff = gff_file
    else:
        f = open(gff_file, "r") if isinstance(gff_file,str) else gff_file
        if "b" in f.mode:
            f = io.TextIOWrapper(f)
        gff = f.readlines()
        f.close()
    
    f_locs = {}
    for line in gff:
        if not( "gene" in line.split("\t") or "mRNA" in line.split("\t") ):
            continue
        for feature in search_terms:
            if feature in line:
                f_locs[feature] = f_locs.get(feature,[])
                f_locs[feature].append( line )
    
    for feature in search_terms:
        if feature not in f_locs:
            print("!"*20)
            print(f"{feature} not found!")
            print("!"*20)
            continue

        if len(f_locs[feature]) == 1:
            f_locs[feature] = f_locs[feature][0]
        else:
            found = f_locs.pop(feature)
            start_locs = set( int(x.split("\t")[3]) for x in found )
            end_locs = set( int(x.split("\t")[4]) for x in found )
            if len(start_locs) == 1 and len(end_locs) == 1:
                f_locs[feature] = found[0]
            else:
                for i, val in enumerate(found):
                    f_locs[f"{feature}.found_{i+1}"] = val
    return f_locs

def get_gene_fa(
        genome_file: str|os.PathLike|typing.IO, 
        locations: dict[str,list[str]],
        offset: int = 0
        ) -> dict[str,tuple[str,str,int,int]]:
    """
    Extracts the gene sequences from the genome file based on the locations
    Args:
        genome_file: The genome file to extract the sequences from
        locations: The locations of the genes
    Returns:
        A dictionary of the genes and their sequence and chromosome/start/stop position of the gene.
    """
    
    chromosomes = { x.split(None,1)[0]:"" for x in locations.values() }
    if chromosomes == {}:
        raise ValueError("No chromosomes found in the locations.")

    f = open(genome_file, "r") if isinstance(genome_file,str) else genome_file
    if "b" in f.mode:
        f = io.TextIOWrapper(f)

    keep = False
    for line in f:
        if line.startswith(">"):
            chromosome = line[1:].strip()
            if chromosome in chromosomes:
                keep = True
            else:
                keep = False
        elif keep:
            chromosomes[chromosome] += line.strip()
    f.close()

    not_found = [x for x,y in chromosomes.items() if y == ""]
    if not_found != []:
        raise ValueError(f"Could not find {', '.join(not_found)} found in the genome file.")

    genes = {}
    for gene, locs in locations.items():
        first = int(1e9)
        last = 0
        chromosome = locs.split("\t")[0]
        vals = [*map(int,locs.split("\t")[3:5])]

        first = min([first,*vals])
        last = max([last,*vals])
        
        chr_range = first-offset,last+offset
        chr_range = [0 if x < 0 else x for x in chr_range]
        genes[gene] = ( chromosomes[chromosome][chr_range[0]:chr_range[1]], chromosome, *chr_range )

    return genes

def extract_seq_from_genome_based_on_gff_feature(
        genome_file: str|os.PathLike|typing.IO, 
        gff_file: str|os.PathLike|typing.IO, 
        search_terms: str|list[str], 
        offset: int = 0,
        ) -> dict[str,tuple[str,int,list[str]]]:
    """
    Extracts the gene sequences from the genome file based on the gff file and search terms.
    Args:
        genome_file: The genome file to extract the sequences from
        gff_file: The gff file to search
        search_terms: The search terms to use
        offset: The offset to use when extracting the sequences
    Returns:
        A dictionary of the genes and their sequence and features.
    """
    locations = find_gene_loc(gff_file, search_terms)
    genes = get_gene_fa(genome_file, locations, offset)
    features = {}
    result = {}
    for x, (seq, chromosome, start, end) in genes.items():
        features[x] = find_features_between(gff_file, chromosome, start, end)
        result[x] = (seq, min(start,end), features[x])
    return result

def reformat_features_for_gb(seq, start, features: list[str]) -> list[dict[str,str]]:
    """
    Reformats the features for a makegb call.
    Args:
        features: The features to reformat
    Returns:
        A list of dictionaries of the features
    """
    f = []
    for line in features:
        vals = line.split("\t")
        s = int(vals[3]) - start
        e = int(vals[4]) - start
        s,e = min(s,e),max(s,e)
        s = 0 if s < 0 else s
        e = 0 if e > len(seq) else e
        f.append({
            "type":vals[2],
            "location":f"{s}..{e}" if vals[6] == "+" else f"complement({s}..{e})",
            "strand":vals[6],
            "notes": vals[8],
        })
    return f

def find_sub_seq(seq: str, sub_seq: str,check_rev_complement=False) -> list[(int,int)]:
    """
    Finds the sub sequence in the sequence
    Args:
        seq: The sequence to search
        sub_seq: The sub sequence to find
    Returns:
        A list of the indexes of the sub sequence in the sequence
    """
    
    seq = seq.lower()
    sub_seq = sub_seq.lower()
    i = seq.find(sub_seq,0)
    found = []
    while i != -1:
        found.append((i,i+len(sub_seq)))
        i = seq.find(sub_seq,i+1)
        
    if check_rev_complement:
        return found + find_sub_seq(reverse_complement(seq),sub_seq,False)
    return found

def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complement of a sequence
    Args:
        seq: The sequence to reverse complement
    Returns:
        The reverse complement of the sequence
    """
    return "".join([complement[x] for x in seq[::-1]])

def main(genome_file: str, gff_file: str, search_term: str, offset: int = 0, output: str | typing.TextIO = sys.stdout):
    locations = find_gene_loc(gff_file, search_term)
    if locations == {}:
        print(f"No {search_term} found in the gff file.")
        return
    genes = get_gene_fa(genome_file, locations, offset)
    fasta_seq = ""
    for gene, (seq, chromosome, start, end) in genes.items():
        features = find_features_between(gff_file, chromosome, start, end)
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

def help():
    print(
        """
    Usage:
        python extract_seq_from_genome_based_on_gff_feature.py <genome_file> <gff_file> <search_term> [offset=0] [output=standard output]
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
    )

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
        help()
        sys.exit()

    args = {i:j for i,j in enumerate(sys.argv)}
    genome_file = args.get(1)
    gff_file = args.get(2)
    search_term = args.get(3)
    offset = int(args.get(4,0))
    output = args.get(5,sys.stdout)
    fasta_seq = main(genome_file, gff_file, search_term, offset, output)
    if fasta:
        print(fasta_seq)
