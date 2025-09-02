"""
This script is used for generating the oligos for any CAS system incorporating Hammerhead to be cloned into a vector (by default it is designed for pNOC-CAS9 and derivatives).
Further, it can also be used to find the guide RNAs for a given target sequence using various CRISPR/Cas systems.

Global Options:
    --mode=find_guide, -m find_guide: Run the script in find_guide mode.
    --homology_5: The 5' homology sequence. Default is tccctccatccacagaatcg.
    --homology_3: The 3' homology sequence. Default is gtaccatgggaaagaaagga.
    -h, --help: Print this help message.

Oligo Generation:
        python guide_helper.py <target_sequence> <case_formatting>
    Arguments:
        target_sequence: The target sequence of the CAS9 site
        case_formatting: True if you want the oligos to swap between upper and lower case, False otherwise

Finding Guide RNA:
        python guide_helper.py <target_sequence> [genome] [--cas_type=0 --start_index=0 --end_index=-1 --output_file=None --guide_blast_result_location=None]
    
    You can also use space separated options instead of =.
    Arguments:
        target_sequence: The target sequence of the CAS9 site. Fasta file or sequence string is accepted.
        genome: The name of blast database to blast against for off-target sites. Optional, blast search will be skipped if not provided.
    Options:
        --cas_type:
            The type of CRISPR/Cas system to use. Can be either an integer or a string. Default is 0.
            Supported CAS types are:
                SpCas9, SpCas9_VRER, SpCas9_VQR, xCas9, SpCas9_NG, SaCas9, AsCpf1, AsCpf1_RR, AsCpf1_RVR, LbCpf1, LbCpf1_RR, FnCas12a, CjCas9, NmCas9, StCas9, TdCas9
        --start_index: The start index of the target sequence. Optional, default is 0.
        --end_index: The end index of the target sequence. Optional, default is -1.
        --output_file: Print out the guide sequences to a file instead of stdout. This includes the detailed ordering information if --detailed_output is provided.
        --guide_blast_result_location: The location to save the blast results for the guide sequences.
        --gb_file: The location to save the genbank file for the target sequence.
        --detailed_output: Print out the detailed ordering information for each guide RNA.
    Special Options:
        --protein_fa:
            The location of the protein fasta file.
            If provided, the script will generate the genbank file for the target sequence and the protein sequence.
            The start and end index will be ignored.
            CDS will be identified based on the protein sequence.
            The protein sequence must be encoded by the target sequence, splice information is not needed.
            If gb_file is not provided, the genbank file will be generated as the fasta filename but with gb extension.
        --gene_start: The start index of the gene in the target sequence. Default is 0.
        --gene_end: The end index of the gene in the target sequence. Default is -1.
        --gene_length: The length of the gene in the target sequence. Ignores gene_end if provided.
        --search_offset: The offset to search for the guide RNA from the start and end of the CDS. Default is 25.

"""

def help():
    print(__doc__)

import os
import sys

from extract_seq_from_genome_based_on_gff_feature import reverse_complement, identifyCDS
from format_genbank import makegb

from re import finditer, IGNORECASE
from shutil import which
from subprocess import Popen, PIPE

from warnings import warn
import math


def gene_frag_to_order(target:str, homology_5: str = "tccctccatccacagaatcgat", homology_3: str = "ggtaccatgggaaagaaagga", case_formatting=False) -> str:
    '''
    Given a target sequence, return the oligos to order
    '''
    sg_RNA_scaffold = "gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttggccggcatggtcccagcctcctcgctggcgccggctgggcaacatgcttcggcatggcgaatgggac"

    non_cutters = {
        "avrII_cutsite": "CCTAGG",
        "asiSI_cutsite": "GCGATCGC",
        "pacI_cutsite" : "TTAATTAA",
        "sgrAI_cutsite": "CACCGGTG",
        "bbvci_cutsite": "CCTCAGC",
        "fseI_cutsite" : "GGCCGGCC",
    }
    
    cutsite_front, cutsite_back = "", ""
    for i, (enzyme, recog_site) in enumerate(non_cutters.items()):
        if recog_site not in sg_RNA_scaffold and recog_site not in target:
            if cutsite_front == "":
                cutsite_front = enzyme, recog_site
            elif cutsite_back == "":
                cutsite_back = enzyme, recog_site
            else:
                break

    guide_with_HH = "$reverse_compliment6$ctgatgagtccgtgaggacgaaacgagtaagctcgtc$target$"

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complement.update({i.lower(): j.lower() for i, j in complement.items()})

    if case_formatting:
        target = target.upper()

    oligo = guide_with_HH.replace("$target$", target)
    target_complement = ''.join(complement[base] for base in target)
    oligo = oligo.replace("$reverse_compliment6$", target_complement[5::-1])
    common_cutters = { "BamHI": "GGATCC", "EcoRI": "GAATTC", "HindIII": "AAGCTT", "XhoI": "CTCGAG" }
    multiple_cloning_site = ''.join([common_cutters[i] for i in common_cutters])

    front_frag = homology_5.upper() + cutsite_front[1].lower() + oligo + sg_RNA_scaffold.lower() + cutsite_front[1].upper() + multiple_cloning_site
    back_frag = multiple_cloning_site + cutsite_back[1].lower() + oligo + sg_RNA_scaffold.lower() + cutsite_back[1].upper() + homology_3.lower()
    direct = homology_5.upper() + cutsite_front[1].lower() + oligo + sg_RNA_scaffold.lower() + cutsite_front[1].upper() + homology_3.lower()
    return front_frag, back_frag, cutsite_front[0], cutsite_back[0], direct

def ordering_info( target:str, homology_5: str = "tccctccatccacagaatcgat", homology_3: str = "ggtaccatgggaaagaaagga", case_formatting=False ) -> str:
    front_insert, back_insert, cutsite_f, cutsite_b, direct_insert = gene_frag_to_order(target, homology_5, homology_3, case_formatting)
    front_primers = pick_primers(front_insert)
    back_primers = pick_primers(back_insert)
    direct_primers = pick_primers(direct_insert)
    import random
    gc_content = sum( [1 for insert in [ front_insert, back_insert, direct_insert ] for nt in insert if nt in "GC"] )
    total = sum( [len(insert) for insert in [ front_insert, back_insert, direct_insert ]] )
    additional_gc = int( 50 * (0.5 - gc_content / total) )
    additional_gc = max(0, additional_gc) // 2
    random_padding = list("AT" * (25-additional_gc) + "GC" * additional_gc)

    random.shuffle(random_padding)
    random_padding = "".join(random_padding)

    return f"""guide RNA: {target}
To insert into {homology_5}...{homology_3} using infusion or similar method (default pNOC-Cas9):
1) Order gene fragments and primers below.
    *Gene fragments may need to be padded to be synthesized by Twist Bioscience.
    * Random padding sequence: {random_padding}
2) PCR amplify the gene fragments using the primers below.
    *Use front and back if inserting pair of guides otherwise use direct insert.
3) Linearize the vector. [pNOC-Cas9 derivatives: Digest with ClaI and KpnI]
4) Use infusion cloning (or similar method) and mix all the fragments with digested vector.
5) Transform the mixture into competent cells and select for the correct clones.
Multiple cloning site is included for front and back inserts, so you can use one of BamHI, EcoRI, HindIII, or XhoI to cut the synthesized fragments to join a pair of guides.
Note that these enzymes cannot be used to cut pNOC-Cas9 and derivatives, use ClaI and KpnI instead.

Front and back inserts can be removed by digesting with {cutsite_f} and {cutsite_b} respectively.
Direct insert can be removed by digesting with {cutsite_f}.

Ordering information GeneFrag:
    Front Insert ({len(front_insert)}): {front_insert}

    Back Insert ({len(back_insert)}): {back_insert}

    Direct Insert ({len(direct_insert)}): {direct_insert}

Oligos to order:
    Front Insert:
        fwd: 5' {front_primers[0]} 3'
        rev: 5' {front_primers[1]} 3'
    Back Insert:
        fwd: 5' {back_primers[0]} 3'
        rev: 5' {back_primers[1]} 3'
    Direct Insert:
        fwd: 5' {direct_primers[0]} 3'
        rev: 5' {direct_primers[1]} 3'
    """

class CAS:
    """
    Class that allows for the identification of guides for a given target sequence using various CRISPR/Cas systems.
    If the supported CAS system does not fit your needs, you can edit the name, pam, length, and location attributes (See Attributes).
    Initialization:
            var = CAS(cas_type)
        Parameters:
            cas_type:
                The type of CRISPR/Cas system to use. Can be either an integer or a string.
                Print CAS.supported_cas for the dictionary of supported types. 
                The integer value is interpreted as the index of the supported_cas dictionary.
    Attributes:
        name: The name of the CRISPR/Cas system
        pam: The list of PAM sequence for the CRISPR/Cas system. Supports degenerate bases.
        length: The length of the guide RNA for the CRISPR/Cas system
        location: The location of the guide RNA in relation to the PAM sequence (5 or 3 representing 5' or 3' respectively)
    Methods:
        find_guides: 
            Given a target sequence, return a list of guide RNAs for the CRISPR/Cas system.
            Usage:
                var = CAS()
                var.find_guides(target_sequence, seq_start_index=0, seq_end_index=-1)
            Parameters:
                target_sequence: The target sequence for the CRISPR/Cas system
                seq_start_index: The start index of the target sequence. Optional, default is 0.
                seq_end_index: The end index of the target sequence. Optional, default is -1.
            Returns:
                A list of guide RNAs for the CRISPR/Cas system. 
                Found guides are in uppercase, reverse complement guides are in lowercase.
        find_guides_and_blast_genome: Given a target sequence and a genome sequence, return a list of guide RNAs for the CRISPR/Cas system and blast the genome for potential off-target sites.
            Usage:
                var = CAS()
                var.find_guides_and_blast_genome(
                    target_sequence, genome,
                    start_index=0, end_index=-1,
                    output_file=None, guide_blast_result_location=None
                )
            Parameters:
                target_sequence, start_index, end_index: See find_guides.
                genome: 
                    The genome database to blast against for off-target sites.
                    Note both blastn and the database must be accessible from the command line.
                output_file: Print out the guide sequences to a file.
                guide_blast_result_location: The location to save the blast results for the guide sequences.
            Returns:
                A list of guide RNAs for the CRISPR/Cas system.
                Blast results are saved in the guide_blast_result_location if provided.
    """
    __SpCas9__ = {"pam":["NGG"], "length":20, "location":"3"}
    __SpCas9_VRER__ = {"pam":["NGCG"], "length":20, "location":"3"}
    __SpCas9_VQR__ = {"pam":["NGAN","NGNG"], "length":20, "location":"3"}
    __xCas9__ = {"pam":["NG","GAA","GAT"], "length":20, "location":"3"}
    __SpCas9_NG__ = {"pam":["NG"], "length":20, "location":"3"}
    __SaCas9__ = {"pam":["NNGRRT"], "length":20, "location":"3"}

    __AsCpf1__ = {"pam":["TTTV"], "length":23, "location":"5"}
    __AsCpf1_RR__ = {"pam":["TYCV"], "length":23, "location":"5"}
    __AsCpf1_RVR__ = {"pam":["TATV"], "length":23, "location":"5"}

    __LbCpf1__ = {"pam":["TTTV"], "length":23, "location":"5"}
    __LbCpf1_RR__ = {"pam":["TYCV"], "length":23, "location":"5"}

    __FnCas12a__ = {"pam":["TTTN"], "length":23, "location":"5"}

    __CjCas9__ = {"pam":["NNNNRYAC"], "length":24, "location":"3"}

    __NmCas9__ = {"pam":["NNNNGATT"], "length":24, "location":"3"}

    __StCas9__ = {"pam":["NNAGAAW"], "length":24, "location":"3"}

    __TdCas9__ = {"pam":["NAAAAC"], "length":24, "location":"3"}

    supported_cas = {
        "SpCas9":__SpCas9__,
        "SpCas9_VRER":__SpCas9_VRER__,
        "SpCas9_VQR":__SpCas9_VQR__,
        "xCas9":__xCas9__,
        "SpCas9_NG":__SpCas9_NG__,
        "SaCas9":__SaCas9__,
        "AsCpf1":__AsCpf1__,
        "AsCpf1_RR":__AsCpf1_RR__,
        "AsCpf1_RVR":__AsCpf1_RVR__,
        "LbCpf1":__LbCpf1__,
        "LbCpf1_RR":__LbCpf1_RR__,
        "FnCas12a":__FnCas12a__,
        "CjCas9":__CjCas9__,
        "NmCas9":__NmCas9__,
        "StCas9":__StCas9__,
        "TdCas9":__TdCas9__
        }

    def __init__(self, cas_type: int|str=0):
        if type(cas_type) == str:
            cas_type = [i.lower() for i in self.supported_cas.keys()].index(cas_type.lower())
        if cas_type >= len(self.supported_cas):
            raise Exception("Unsupported CAS type")
        self.name = list(self.supported_cas.keys())[cas_type]
        self.pam = self.supported_cas[self.name]["pam"]
        self.length = self.supported_cas[self.name]["length"]
        self.location = self.supported_cas[self.name]["location"]
    
    @property
    def search_formats(self):
        pam = [i.replace("N","[ACTG]") for i in self.pam]
        pam = [i.replace("R","[AG]") for i in pam]
        pam = [i.replace("Y","[CT]") for i in pam]
        pam = [i.replace("W","[AT]") for i in pam]
        pam = [i.replace("S","[GC]") for i in pam]
        pam = [i.replace("K","[GT]") for i in pam]
        pam = [i.replace("M","[AC]") for i in pam]
        pam = [i.replace("B","[CGT]") for i in pam]
        pam = [i.replace("D","[AGT]") for i in pam]
        pam = [i.replace("H","[ACT]") for i in pam]
        pam = [i.replace("V","[ACG]") for i in pam]
        if self.location == "5":
            pam = [f"(?=({i})({'.'*self.length}))" for i in pam]
        else:
            pam = [f"(?=({'.'*self.length})({i}))" for i in pam]
        return pam

    def find_guides(self, seq: str, start:int=0, end:int=-1, include_pams=False) -> list[str]:
        """
        For the defined CAS system, find the guide RNAs for the target sequence.
        """
        if set(seq) - set("ATCG") != set():
            raise Exception("Invalid sequence")
        seq = seq[start:end]
        rev_seq = reverse_complement( seq )
        seq, rev_seq = str(seq), str(rev_seq)
        for pam in self.search_formats:
            fguides = finditer(pam, seq, flags=IGNORECASE)
            rguides = finditer(pam, rev_seq, flags=IGNORECASE)
            pam_group, guide_group = (0, 1) if self.location == "5" else (1, 0)
            guides = [ [j.upper() for j in i.groups()] for i in fguides ] 
            guides += [ [j.lower() for j in i.groups()] for i in rguides ]
            pams = [ guide[pam_group] for guide in guides ]
            guides = [ guide[guide_group] for guide in guides ]
            if include_pams:
                guides = [ 
                    pam+guide if self.location == "5" else guide+pam
                    for pam, guide in zip(pams, guides) ]
        return guides

    def find_guides_and_blast_genome(self, 
        target:str, genome: str,
        start:int=0, end:int=-1,
        guide_blast_result_location: str|os.PathLike=None,
        blast_cutoff:float=1.0
    ):
        """
        For the defined CAS system, find the guide RNAs for the target sequence blast the genome for potential off-target sites.
        """
        # Find guides using for the target sequence.
        guides = [ (i,guide) for i,guide in enumerate(self.find_guides(target, start, end, True))]
        actual_guides = [ i for i in self.find_guides(target, start, end) ]

        # Look for the blastn executable.
        blastn = which("blastn")
        if not blastn:
            warn("Blastn not found. Skipping BLAST search.")
            return [ actual_guides[i] for i,_ in guides ]
        
        # Check if the genome database exists.
        if guide_blast_result_location and ( isinstance(guide_blast_result_location, str) or isinstance(guide_blast_result_location, os.PathLike) ):
            guide_blast_result_location = os.path.abspath(guide_blast_result_location)
            os.makedirs(guide_blast_result_location, exist_ok=True)
        else:
            guide_blast_result_location = "."
        guide_index = 1
        fail_index = 1
        selected_guides = []
        for guide_iteration, guide in guides:
            temp_guide_loc = os.path.join(guide_blast_result_location, "temp_guide.fa")
            with( open( temp_guide_loc, "w") ) as f:
                f.write(f">guide_{guide_iteration}\n{guide}")
            
            process = Popen( [blastn, "-db", genome,  "-query", temp_guide_loc, "-outfmt", "6 qacc sacc evalue qseq sseq qstart qend", "-ungapped", "-task", "blastn-short"], stdout=PIPE, stderr=PIPE)
            result, error = process.communicate()
            result = result.decode("utf-8")
            result = [ i for i in result.split("\n") if i != "" if int(i.split("\t")[-1]) == len(guide) ]
            result = [ i for i in result if float(i.split("\t")[2]) < blast_cutoff ]
            if result:
                selected_guides.append(guide_iteration)
                res_loc = os.path.join(guide_blast_result_location, f"guide_{guide_index}.blast.txt")
                guide_index += 1
            else:
                res_loc = os.path.join(guide_blast_result_location, f"failed_guide_{fail_index}.blast.txt")
                fail_index += 1
            
            if guide_blast_result_location != ".":
                res_file = open(res_loc, "w")
                process = Popen( [blastn, "-db", genome,  "-query", temp_guide_loc, "-ungapped", "-task", "blastn-short"], stdout=res_file, stderr=PIPE)
                result, error = process.communicate( timeout=60)
                res_file.close()
                
        if os.path.exists(temp_guide_loc):
            os.remove(temp_guide_loc)
        
        return [ actual_guides[i] for i in selected_guides ]
    
def calculate_tm(sequence, salt_concentration=50e-3, dna_concentration=0.5e-6):

    sequence = sequence.upper()
    # Thermodynamic parameters for nearest-neighbor pairs (ΔH in kJ/mol and ΔS in J/(mol·K))
    nn_params = {
        "AA": (-33.1, -92.9), "TT": (-33.1, -92.9),
        "AT": (-30.1, -85.4), "TA": (-30.1, -89.1),
        "CA": (-35.6, -95.0), "TG": (-35.6, -95.0),
        "GT": (-35.1, -93.7), "AC": (-35.1, -93.7),
        "CT": (-32.6, -87.9), "AG": (-32.6, -87.9),
        "GA": (-34.3, -92.9), "TC": (-34.3, -92.9),
        "CG": (-44.4, -113.8), "GC": (-41.0, -102.1),
        "GG": (-33.5, -83.3), "CC": (-33.5, -83.3)
    }
    # Terminal base pair adjustments (ΔH in kJ/mol and ΔS in J/(mol·K))
    terminal_params = {
        "A": (9.6, 17.2), "T": (9.6, 17.2),
        "G": (0.4, -11.7), "C": (0.4, -11.7)
    }

    dH = 0
    dS = 0
    for i,j in zip( sequence, sequence[1:]):
        dH += nn_params[i+j][0]
        dS += nn_params[i+j][1]

    dH += terminal_params[sequence[0]][0]
    dS += terminal_params[sequence[0]][1]
    dH += terminal_params[sequence[-1]][0]
    dS += terminal_params[sequence[-1]][1]

    # Melting temperature calculation
    tm = (dH * 1000) / ( dS + 8.314 * math.log( dna_concentration ) ) - 273.15
    tm += 16.6 * math.log10(salt_concentration)
    return tm
    
    
def pick_primers(sequence, min_tm=50):
    rev_seq = reverse_complement(sequence)
    fwd = sequence[:15]
    rev = rev_seq[:15]
    while True:
        fwd_tm = calculate_tm(fwd)
        rev_tm = calculate_tm(rev)
        if fwd_tm < min_tm or fwd[-1].upper() in "AT":
            fwd += sequence[len(fwd)]
        elif rev_tm - fwd_tm > 3:
            fwd += sequence[len(fwd)]

        if rev_tm < min_tm or rev[-1].upper() in "AT":
            rev += rev_seq[len(rev)]
        elif fwd_tm - rev_tm > 3:
            rev += rev_seq[len(rev)]
        if abs(fwd_tm - rev_tm) < 5 and fwd_tm > min_tm and rev_tm > min_tm and fwd[-1] not in "AT" and rev[-1].upper() not in "AT":
            break
    return fwd, rev

def make_guide_gb(dna_seq, guides, seq_name, homology_5, homology_3,  other_annotations={}, output=sys.stdout):
    """
    Generates a GenBank file with guide RNA annotations DNA sequence.
    Args:
        dna_seq: The DNA sequence to be annotated.
        guides: A list of guide RNA sequences.
        seq_name: The name of the DNA sequence.
        other_annotations: A dictionary of other annotations to be included in the GenBank file.
        output: The file to write the GenBank file to. Default is sys.stdout.
    Returns:
        The text of the GenBank file.
    """

    for i, guide in enumerate(guides):
        complement = False
        guide_loc = dna_seq.lower().find(guide.lower())
        if guide_loc == -1:
            guide_loc = dna_seq.lower().find(reverse_complement(guide.lower()))
            if guide_loc == -1:
                raise ValueError(f"Guide {guide} not found in sequence\n{dna_seq}")
            complement = True
        
        start = guide_loc + 1
        end = start + len(guide) - 1

        other_annotations.append(
            {
                "location": f"{start}..{end}" if not complement else f"complement({start}..{end})",
                "label": f"guide RNA {i+1}",
                "type": "misc_RNA",
                "note": f'"{ordering_info(guide, homology_5, homology_3, True)}"',
            }
        )
        other_annotations[-1]["note"] = other_annotations[-1]["note"].replace("\n", "\n" + " " * 2)
    gb =  makegb(
            gene=seq_name,
            organism="Nannochloropsis",
            seq= dna_seq,
            features=other_annotations,
        )
    if type(output) == str:
        f = open(output, "w")
    else:
        f = output
    
    print( gb, file=f )

    if type(output) == str:
        f.close()
    
    return output

def parse_options(argv):
    """
    Parse the options from the command line.
    """
    parsed_options = {}
    remove = []
    for i in argv:
        if i.startswith("-") and "=" in i:
            parsed_options[i.split("=")[0]] = i.split("=")[1]
            remove.append(i)
    for i, j in zip(argv, argv[1:]):
        if i.startswith("-") and not "=" in i:
            remove += [i, j]
            i = "--" + i.replace("-","")
            parsed_options[i] = j
            

    for i in remove:
        argv.remove(i) if i in argv else None

    for i in argv:
        if i.startswith("-") and i not in remove:
            i = "--" + i.replace("-","")
            parsed_options[i] = True
            remove.append(i)
    argv = [i for i in argv if not i.startswith("-")]

    parsed_options["target"] = argv[1]
    parsed_options["find_guide"] = parsed_options.get("--mode", "") == "find_guide" or parsed_options.get("--m", "") == "find_guide"

    genome = argv[2] if len(argv) > 2 else ""
    if genome == "":
        return parsed_options

    db_locations = os.environ.get("BLASTDB",".") + ":."
    db_locations = [ f for db_location in db_locations.split(":") for f in os.listdir(db_location) if db_location != "." ]
    required_files = [genome+".nhr", genome+".nin", genome+".nsq"]
    if not all( [f in db_locations for f in required_files] ):
        warn(f"Genome database {genome} not found in BLASTDB environment variable or current directory. Skipping BLAST search.")
    else:
        parsed_options["genome"] = genome
    return parsed_options

def main(argv):

    parsed_options = parse_options(argv)
    if parsed_options.get("--help", False) or parsed_options.get("-h", False):
        help()
        return

    target = parsed_options["target"]
    if os.path.exists(target):
        target = "".join( open(target).readlines() )
        target = [i for i in target.split(">") if i != ""][0]
        target = "".join( target.split("\n")[1:] )
    elif set(target.upper()) - set("ATCG") != set():
        raise Exception("Invalid sequence or file for target sequence.")
    homology_5 = parsed_options.get("--homology_5", "tccctccatccacagaatcgat")
    homology_3 = parsed_options.get("--homology_3", "ggtaccatgggaaagaaagga")
    case_formatting = parsed_options.get("--case_formatting", False)
    if not parsed_options.get("find_guide", False):
        if len(target) > 30:
            print("Target sequence is too long. I don't think this is the function you are trying to run.")
            import time
            time.sleep(1)
            print("Here you go anyway:")
            time.sleep(2)
        print( ordering_info(target, homology_5, homology_3, case_formatting) )
        return
    
    target_name = ".".join(parsed_options["target"].split(".")[:-1])
    genome = parsed_options.get("genome",None)
    cas = CAS(parsed_options.get("--cas_type",0))
    start = int(parsed_options.get("--start_index",0))
    end = int(parsed_options.get("--end_index",-1))
    output_file = parsed_options.get("--output_file",sys.stdout)
    detailed_output = parsed_options.get("--detailed_output",False)
    guide_blast_result_location = parsed_options.get("--guide_blast_result_location",None)
    homology_5 = parsed_options.get("--homology_5","tccctccatccacagaatcgat")
    homology_3 = parsed_options.get("--homology_3","ggtaccatgggaaagaaagga")
    protein_fa = parsed_options.get("--protein_fa",None)
    gb_file = parsed_options.get("--gb_file",None)
    search_offset = int(parsed_options.get("--search_offset",25))
    homology_5 = parsed_options.get("--homology_5","tccctccatccacagaatcgat")
    homology_3 = parsed_options.get("--homology_3","ggtaccatgggaaagaaagga")

    def get_guides(s, e):
        print(f"Finding guides for {cas.name} system in target sequence from {s} to {e}.")
        if genome is not None and genome != "":
            guides = cas.find_guides_and_blast_genome(
                target, genome,
                s, e,
                guide_blast_result_location
            )
            for i in os.listdir(guide_blast_result_location):
                if i.endswith(".blast.txt") and (i.startswith("guide_") or i.startswith("failed_guide_")):
                    os.rename( 
                        os.path.join(guide_blast_result_location, i), 
                        os.path.join(guide_blast_result_location, f"{target_name}.{i}") )
        else:
            guides = cas.find_guides(target, s, e)
        return guides

    
    guides = []
    if protein_fa is not None:
        protein_seq = "".join( open(protein_fa).readlines() )
        protein_seq = [i for i in protein_seq.split(">") if i != ""][0]
        protein_seq = "".join( protein_seq.split("\n")[1:] )

        gene_start = int(parsed_options.get("--gene_start",0))
        gene_end = int(parsed_options.get("--gene_end",-1))
        if "gene_length" in parsed_options:
            gene_end = gene_start + int(parsed_options["gene_length"])
        try:
            cds_loc = identifyCDS(protein_seq, target[gene_start:gene_end])
        except ValueError:
            print(f"Protein sequence({len(protein_seq)}) does not match target sequence({len(target)}).")
            print("Check the protein sequence and template sequence.")
            print("Also try providing the location of the gene in the target sequence using --gene_start=X and --gene_end=Y.")
            sys.exit(1)
        
        # Add the offset for each CDS
        for cds in cds_loc:
            loc = cds["location"].replace("complement(","").replace(")","")
            start, end = map(int, loc.split(".."))
            start, end = start+gene_start, end+gene_start
            if cds["strand"] == "-":
                loc = f"complement({start}..{end})"
            else:
                loc = f"{start}..{end}"
            cds["location"] = loc

            if cds["translation"].startswith( protein_seq[:5]):
                start_codon = cds["location"].split("..")
            if cds["translation"].endswith( protein_seq[-5:]):
                stop_codon = cds["location"].split("..")
        
        if "complement" in start_codon[0]:
            start = int( start_codon[1].split( ")" )[0] )
            end = int( stop_codon[0].split( "(" )[1] )
        else:
            start = int(start_codon[0])
            end = int(stop_codon[1])
        
        guides = get_guides(start-search_offset, start+search_offset)
        start, end = end-search_offset, end+search_offset
        if gb_file is None:
            gb_file = os.path.basename(protein_fa).replace(".fa",".gb")
            gb_file = os.path.join( os.path.dirname(output_file), gb_file )
    else:
        cds_loc = {}
    
    guides += get_guides(start, end)

    if gb_file is not None:
        make_guide_gb(target, guides, target_name, homology_5, homology_3, cds_loc, gb_file)
    
    if output_file is not None and output_file != sys.stdout and output_file != "":
        f = open(output_file, "w")
    else:
        f = sys.stdout
    for i, guide in enumerate(guides):
        print(f"Guide {i+1}: {guide}", file=f)
        if detailed_output:
            print(ordering_info(guide, homology_5, homology_3, True), file=f)
    if output_file is not None:
        f.close()

if __name__ == "__main__":
    if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
        help()
        sys.exit(1)

    argv = sys.argv.copy()
    main(argv)

