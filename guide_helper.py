"""
This script is used for generating the oligos for the CAS9 system incorporating Hammerhead to be cloned in pNOC-Cas9 vector.
Further, it can also be used to find the guide RNAs for a given target sequence using various CRISPR/Cas systems.

Oligo Generation:
        python generate_CAS9_primers.py <target_sequence> <case_formatting>
    Arguments:
        target_sequence: The target sequence of the CAS9 site
        case_formatting: True if you want the oligos to swap between upper and lower case, False otherwise

Finding Guide RNA:
        python generate_CAS9_primers.py <target_sequence> [genome] [--cas_type=0 --start_index=0 --end_index=-1 --output_file=None --guide_blast_result_location=None]
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
        --output_file: Print out the guide sequences to a file.
        --guide_blast_result_location: The location to save the blast results for the guide sequences.
        --gb_file: The location to save the genbank file for the target sequence.
    Special Options:
        --protein_fa:
            The location of the protein fasta file.
            If provided, the script will generate the genbank file for the target sequence and the protein sequence.
            The start and end index will be ignored.
            CDS will be identified based on the protein sequence.
            The protein sequence must be encoded by the target sequence, splice information is not needed.
            If gb_file is not provided, the genbank file will be generated as the fasta filename but with gb extension.
        --search_offset: The offset to search for the guide RNA from the start and end of the CDS. Default is 25.

Global Options:
    --mode=find_guide, -m find_guide: Run the script in find_guide mode.
    -h, --help: Print this help message.
"""

def help():
    print(__doc__)

import os
import sys

from extract_seq_from_genome_based_on_gff_feature import reverse_complement, identifyCDS
from format_genbank import makegb

from re import findall, IGNORECASE
from shutil import which
from subprocess import Popen, PIPE

from warnings import warn
from typing import TextIO


def oligos_to_order(target:str, case_formatting=False) -> str:
    '''
    Given a target sequence, return the oligos to order
    '''
    five_prime = "cga$reverse_compliment6$ctgatgagtccgtgaggacgaaacgagtaagctcgtc$target$g"
    three_prime = "$target6_reverse$gactactcaggcactcctgctttgctcattcgagcag$target_compliment$caaa"

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complement.update({i.lower(): j.lower() for i, j in complement.items()})
    if case_formatting:
        target = target.upper()
    oligo5 = five_prime.replace("$target$", target)
    target_complement = ''.join(complement[base] for base in target)
    oligo5 = oligo5.replace("$reverse_compliment6$", target_complement[5::-1])
    
    oligo3 = three_prime.replace("$target_compliment$", target_complement)
    oligo3 = oligo3.replace("$target6_reverse$", target[5::-1])

    return oligo5, oligo3

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
        pam = [f"{'.'*self.length}{i}" for i in pam]
        return pam

    def find_guides(self, seq: str, start:int=0, end:int=-1) -> list[str]:
        """
        For the defined CAS system, find the guide RNAs for the target sequence.
        """
        if set(seq) - set("ATCG") != set():
            raise Exception("Invalid sequence")
        rev_seq = reverse_complement( seq )
        seq, rev_seq = str(seq), str(rev_seq)
        for pam in self.search_formats:
            fguides = findall(pam, seq, flags=IGNORECASE)
            rguides = findall(pam, rev_seq, flags=IGNORECASE)
            guides = [ i.upper() for i in fguides ] + [ i.lower() for i in rguides ]
        return guides
    

    def find_guides_and_blast_genome(self, 
        target:str, genome: str,
        start:int=0, end:int=-1,
        output_file: str|os.PathLike|TextIO|None=None,
        guide_blast_result_location: str|os.PathLike=None
    ):
        """
        For the defined CAS system, find the guide RNAs for the target sequence blast the genome for potential off-target sites.
        """
        # Find guides using for the target sequence.
        guides = self.find_guides(target, start, end)

        # Look for the blastn executable.
        blastn = which("blastn")
        if not blastn:
            warn("Blastn not found. Skipping BLAST search.")
            return guides
        
        # If the output file is specified, create a new file to write the guides to the file.
        if output_file and ( isinstance(output_file, str) or isinstance(output_file, os.PathLike) ):
            output_file = open(output_file, "w") 
        if guide_blast_result_location and ( isinstance(guide_blast_result_location, str) or isinstance(guide_blast_result_location, os.PathLike) ):
            guide_blast_result_location = os.path.abspath(guide_blast_result_location)
            os.makedirs(guide_blast_result_location, exist_ok=True)
        else:
            guide_blast_result_location = "."
        guide_index = 1
        for guide in guides:
            with( open( os.path.join(guide_blast_result_location,f"guide_{guide_index}.fa"), "w") ) as f:
                f.write(f">guide_{guide_index}\n{guide}")
            
            process = Popen( [blastn, "-db", genome,  "-query", "temp_guide.fa", "-outfmt", "6 qacc sacc evalue qseq sseq qstart qend", "-ungapped", "-task", "blastn-short"], stdout=PIPE, stderr=PIPE)
            result, error = process.communicate( timeout=60)
            result = result.decode("utf-8").split("\n")
            result = [ i for i in result if i != "" if int(i.split("\t")[-1]) == len(guide) if float(i.split("\t")[2]) < 1 ]
            result = len(result) > 1
            if result:
                os.remove( os.path.join(guide_blast_result_location,f"guide_{guide_index}.fa") )
                continue
            
            if guide_blast_result_location == ".":
                os.remove( os.path.join(guide_blast_result_location,f"guide_{guide_index}.fa") )
            else:
                res_loc = os.path.join(guide_blast_result_location, f"guide_{guide_index}.blast.txt")
                res_file = open(res_loc, "w")
                process = Popen( [blastn, "-db", genome,  "-query", "temp_guide.fa", "-ungapped", "-task", "blastn-short"], stdout=res_file, stderr=PIPE)
                result, error = process.communicate( timeout=60)
                result = result.decode("utf-8").split("\n")
                res_file.close()
                guide_index += 1
        return guides

def make_guide_gb(dna_seq, guides, seq_name, other_annotations={}, output=sys.stdout):
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

    for guide in guides:
        complement = False
        guide_loc = dna_seq.find(guide)
        if guide_loc == -1:
            guide_loc = dna_seq.find(reverse_complement(guide))
            if guide_loc == -1:
                raise ValueError(f"Guide {guide} not found in sequence")
            complement = True
        
        start = guide_loc + 1
        end = start + len(guide) - 1

        fwd, rev = oligos_to_order(guide,True)

        other_annotations.append(
            {
                "location": f"{start}..{end}" if not complement else f"complement({start}..{end})",
                "label": "guide RNA",
                "type": "misc_RNA",
                "note": f"""guide RNA: {guide}
                Oligos to order:
                    {fwd=}
                    {rev=}"""
            }
        )
    gb =  makegb(
            gene=seq_name,
            organism="Nannochloropsis",
            seq= dna_seq,
            features=other_annotations,
        )
    if type(output) == str:
        f = open(output, "w")    
    
    print( gb, file=output )

    if type(output) == str:
        f.close()
    
    return output

if __name__ == "__main__":
    if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
        help()
        sys.exit(1)
    find_guide = False
    if "--mode=find_guide" in sys.argv:
        find_guide = True
        sys.argv.remove("--mode=find_guide")
    elif "-m\tfind_guide" in "\t".join(sys.argv):
        find_guide = True
        sys.args.remove("-m")
        sys.args.remove("find_guide")
    if find_guide:
        if len(sys.argv) < 2:
            help()
            sys.exit(1)
        if len(sys.argv) < 3:
            sys.argv.append("")
        target, genome = sys.argv[1], sys.argv[2]
        parsed_options = {
            i.split("=")[0]:i.split("=")[1] for i in sys.argv
            if i.startswith("--") and "=" in i
        }
        if os.path.exists(target):
            target = "".join( open(target).readlines() )
            target = [i for i in target.split(">") if i != ""][1]
            target = "".join( target.split("\n")[1:] )
        elif set(target) - set("ATCG") != set():
            raise Exception("Invalid sequence or file for target sequence.")
        db_locations = os.environ.get("BLASTDB",".") + ":."
        db_locations = [ f for db_location in db_locations.split(":") for f in os.listdir(db_location) if db_location != "." ]
        required_files = [genome+".nhr", genome+".nin", genome+".nsq"]
        if not all( i in db_locations for i in required_files ):
            warn("Genome database not found. Skipping BLAST search.")
            genome = None

        cas = CAS(parsed_options.get("--cas_type",0))
        start = int(parsed_options.get("--start_index",0))
        end = int(parsed_options.get("--end_index",-1))
        output_file = parsed_options.get("--output_file",sys.stdout)
        guide_blast_result_location = parsed_options.get("--guide_blast_result_location",None)
        
        def get_guides(start, end):
            if genome is not None or genome != "":
                guides = cas.find_guides_and_blast_genome(
                    target, genome,
                    start, end,
                    output_file, guide_blast_result_location
                )
            else:
                guides = cas.find_guides(target, start, end)
            return guides

        protein_fa = parsed_options.get("--protein_fa",None)
        gb_file = parsed_options.get("--gb_file",None)
        search_offset = int(parsed_options.get("--search_offset",25))
        
        guides = []
        
        if protein_fa is not None:
            protein_seq = "".join( open(protein_fa).readlines() )
            protein_seq = [i for i in protein_seq.split(">") if i != ""][1]
            protein_seq = "".join( protein_seq.split("\n")[1:] )
            
            cds_loc = identifyCDS(protein_seq, target)
            start_codon = cds_loc[0]["location"].split("..")
            stop_codon = cds_loc[-1]["location"].split("..")
            
            
            if "complement" in start_codon[0]:
                start = int( start_codon[1].split( ")" )[0] )
                end = int( stop_codon[0].split( "(" )[1] )
            else:
                start = int(start_codon[0])
                end = int(stop_codon[1])
            guides = get_guides(start-search_offset, start+search_offset)
            start , end = end-search_offset, end+search_offset
            if gb_file is None:
                gb_file = os.path.basename(protein_fa).replace(".fa",".gb")
                gb_file = os.path.join( os.path.dirname(output_file), gb_file )
        else:
            cds_loc = {}
        
        guides += get_guides(start, end)

        if gb_file is not None:
            make_guide_gb(target, guides, os.path.basename(gb_file), cds_loc, gb_file)
        

        for i, guide in enumerate(guides):
            print(f"Guide {i+1}: {guide}", file=output_file)
        sys.exit(0) 

    target = sys.argv[1]
    case_formatting = False if len(sys.argv) < 3 else sys.argv[2].lower() == "true"
    five, three = oligos_to_order(target, case_formatting)
    print(f"5' Oligo: {five}")
    print(f"3' Oligo: {three}")

