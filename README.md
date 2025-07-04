# Genomic Tools Guide  

This guide provides instructions on how to use each Python file in this repository, both as standalone scripts and by importing their functions into other Python programs.  

---

## 1. **`smith_waterman.py`**  
### Standalone Usage:  
```bash  
python smith_waterman.py <sequence1> <sequence2> [options]  
```  
#### Options:  
- `--match_score=<number>`: Set the match score (default: 3).  
- `--mismatch_score=<number>`: Set the mismatch penalty (default: -3).  
- `--gap_penalty=<number>`: Set the gap penalty (default: -2).  
- `--characters_per_line=<number>`: Set the number of characters per line in the alignment output (default: 60).  

### Importing Functions:  
```python  
from smith_waterman import smith_waterman  

seq1 = "ACGT"  
seq2 = "AGCT"  
aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2)  
print(aligned_seq1, aligned_seq2, score)  
```  

---

## 2. **`select_genes.py`**  
### Standalone Usage:  
```bash  
python select_genes.py  
```  
This script extracts specific genes from a genome file based on GFF features.  

### Importing Functions:  
```python  
from select_genes import extract_seq_from_genome_based_on_gff_feature  

genes = extract_seq_from_genome_based_on_gff_feature(genome_file, gff_file, search_term, offset=3000)  
print(genes)  
```  

---

## 3. **`plot_gel.py`**  
### Standalone Usage:  
```bash  
python plot_gel.py <file1> <file2> ... <restriction_enzymes> [options]  
```  
#### Options:  
- `--label`: Add a label to the gel plot.  
- `--glob`: Use glob patterns to select files.  

### Importing Functions:  
```python  
from plot_gel import plot_gel  

fig, ax = plt.subplots()  
plot_gel(ax, "example.gb", ["EcoRI", "BamHI"])  
plt.show()  
```  

---

## 4. **`guide_helper.py`**  
### Standalone Usage:  
#### Oligo Generation:  
```bash  
python guide_helper.py <target_sequence> <case_formatting>  
```  
#### Finding Guide RNA:  
```bash  
python guide_helper.py <target_sequence> [genome] [options]  
```  
#### Options:  
- `--cas_type=<type>`: Specify the CRISPR/Cas system.  
- `--output_file=<path>`: Save guide sequences to a file.  

### Importing Functions:  
```python  
from guide_helper import gene_frag_to_order  

oligos = gene_frag_to_order(target_sequence, case_formatting=True)  
print(oligos)  
```  

---

## 5. **`gff_genome_to_annfa.py`**  
### Standalone Usage:  
```bash  
python gff_genome_to_annfa.py <gff_file> <genome_file> [output_file]  
```  

### Importing Functions:  
```python  
from gff_genome_to_annfa import gff_genome_to_annfa  

gff_genome_to_annfa("example.gff", "genome.fasta", "output.fa")  
```  

---

## 6. **`generate_CAS9_primers.py`**  
### Standalone Usage:  
```bash  
python generate_CAS9_primers.py <target_sequence> <case_formatting>  
```  

### Importing Functions:  
```python  
from generate_CAS9_primers import gene_frag_to_order  

primers = gene_frag_to_order(target_sequence, case_formatting=True)  
print(primers)  
```  

---

## 7. **`gen_seq_from_motifs.py`**  
### Standalone Usage:  
This script does not have standalone usage.  

### Importing Functions:  
```python  
from gen_seq_from_motifs import gen_seq_from_motifs  

motifs = gen_seq_from_motifs("example.motif")  
print(motifs)  
```  

---

## 8. **`gen_gel_map.py`**  
### Standalone Usage:  
```bash  
python gen_gel_map.py <file_path> <restriction_enzymes>  
```  

### Importing Functions:  
```python  
from gen_gel_map import gen_frag_sizes  

gel_map = gen_frag_sizes("example.gb", ["EcoRI", "BamHI"])  
gel_map.gen_frag_sizes()  
print(gel_map.frag_sizes)  
```  

---

## 9. **`format_genbank.py`**  
### Standalone Usage:  
This script does not have standalone usage.  

### Importing Functions:  
```python  
from format_genbank import makegb  

genbank_file = makegb(locus="Example", seq="ATGC", features=[{"type": "gene", "location": "1..4"}])  
print(genbank_file)  
```  

---

This guide provides a quick reference for using each script effectively. For more details, refer to the script comments or documentation.  