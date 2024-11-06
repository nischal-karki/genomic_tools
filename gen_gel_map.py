"""
This script takes a genbank, fasta, or snapgene file and a list of restriction enzymes as input and returns a gel map of the restriction sites.
"""
import sys
import warnings
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis
warnings.filterwarnings("ignore")

class gen_frag_sizes:
    def __init__(self, file_path, restriction_enzymes):
        file_format = {"gb": "genbank", "fasta": "fasta", "dna": "snapgene"}[file_path.split(".")[-1]]
        self.seq = SeqIO.read(file_path, file_format)
        self.rb = RestrictionBatch(restriction_enzymes)
        self.frag_sites = []
        self.frag_sizes = []
    
    def gen_frag_sizes(self):
        loc = min(
            [ 
                (site,j) 
                for site, i in self.rb.search(self.seq.seq, linear=False).items() 
                for j in i
            ], 
            key=lambda x: x[1]
        )
        self.rebuild_seq(loc[1])
        frag_sizes, sites = [], []
        new_locs = [(loc[0],0)] +  [ 
                (site,j) 
                for site, i in self.rb.search(self.seq.seq, linear=False).items() 
                for j in i
            ]
        new_locs.sort(key=lambda x: x[1])

        for count, (site,new_loc) in enumerate(new_locs[1:]):
            frag_sizes.append(new_loc - new_locs[count][1])
            sites.append((site,new_locs[count][0]))

        self.frag_sizes = frag_sizes
        self.frag_sites = sites

    def __str__(self):
        if not self.frag_sizes:
            self.gen_frag_sizes()
        out = ""
        for size, sites in zip(self.frag_sizes, self.frag_sites):
            out += f"{sites[0]} to {sites[1]}: {size}\n"
        return out

    def rebuild_seq(self, loc):
        self.seq = self.seq[loc:] + self.seq[:loc]

if __name__ == "__main__":
    fp = sys.argv[1]
    re = sys.argv[2:]
    print(gen_frag_sizes(fp, re))
