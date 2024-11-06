import sys
import os
if os.path.abspath("../Useful Python/") not in sys.path:
    sys.path.insert(0, os.path.abspath("../Useful Python/"))
from multiprocessing import Pool
from extract_seq_from_genome_based_on_gff_feature import *
from format_genbank import makegb
from gen_seq_from_motifs import *
"""
FAD d9  -> 593563
Aureo4  -> 622551 (paper says FAD d6)
FAD d5  -> 593992
FAD d6  -> 606752
FAD d12 -> 591536 (paper says -> 591534 but 591534 is MGDG synthase)
FAD w3  -> 639463
"""

features = {
    "FAD d9": "proteinId=593563",
    "Aureo4": "proteinId=622551",
    "FAD d5": "proteinId=593992",
    "FAD d6": "proteinId=606752",
    "FAD d12": "proteinId=591536",
    "FAD w3": "proteinId=639463"
}
genome_file = "Nanoce1779_2_AssemblyScaffolds.fasta"
gff_file = "Nanoce1779_2_GeneCatalog_20180119.gff3"
offset = 3000
fa = open("Selected_genes/selected.fa","w")

for feature, search_term in features.items():
    genes = extract_seq_from_genome_based_on_gff_feature(genome_file, gff_file, search_term, offset)
    for gene, (seq, start, features) in genes.items():
        print(f"Gene: {gene}",flush=True)
        features = reformat_features_for_gb(seq, start, features)
        motifs = gen_seq_from_motifs("CTACGTCAGC.motif")
        print(f"Identifying {len(motifs)} motifs in {gene}",flush=True)
        def find_and_list_features(motif):
            motif,probability = motif
            if probability < 0.01:
                return []
            motif_r = reverse_complement(motif)
            locs = find_sub_seq(seq, motif)
            rlocs = find_sub_seq(seq, motif)
            locs += [(j,i) for i,j in rlocs]
            return [{
                "type": "misc_feature",
                "label": f"CTACGTCAGC like motif: {motif}",
                "location": f"{loc[0]}..{loc[1]}" if loc[0] < loc[1] else f"complement({loc[1]}..{loc[0]})",
                "note": f"{motif} probability: {probability}"
            } for loc in locs]
        motif_features = map(find_and_list_features,motifs.items())
        features += [motif_feature for motif in motif_features for motif_feature in motif]
        print(f"Writing to file with {len(features)} features.",flush=True)
        with open(f"Selected_genes/{feature}.gb", "w") as f:
            f.write(makegb(locus=gene, seq=seq, features=features))
        fa.write(f">{gene}\n{seq}\n")
        print(f"Saved {feature}.gb")
fa.close()