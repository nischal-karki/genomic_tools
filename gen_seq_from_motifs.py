"""
"""
from itertools import product
from math import prod

def gen_seq_from_motifs(motif_files: str|list[str]) -> dict[str,float]:
    """
    Generate sequences from motifs
    """
    import random

    # Create a list of motifs
    if isinstance(motif_files, str):
        motif_files = [motif_files]

    for motif_file in motif_files:
        with open(motif_file, "r") as f:
            motif_probs = {"A":[],"C":[],"G":[],"T":[]}
            for line in f.readlines():
                if line[0] == ">":
                    continue
                vals = line.split()
                try:
                    cur = {i:float(j) for i,j in zip("ACGT",vals)}
                    cur = {i:j if j > 0.001 else 0 for i,j in cur.items()}
                    total = sum(cur.values())
                    cur = {i:j/total for i,j in cur.items()}
                    for i in cur:
                        motif_probs[i].append(cur[i])
                except ValueError:
                    continue

        motif_len = len(motif_probs["A"])
        motifs = {
            i:prod( motif_probs[k][l] for l,k in enumerate(i) ) 
                for i in product("ACGT",repeat=motif_len)
            }
        motifs = {"".join(k):v*100 for k,v in motifs.items() if v > 1e-5}
    return motifs