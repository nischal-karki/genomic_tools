# Code generated by CoPilot based on Smith-Waterman algorithm.

def smith_waterman(seq1, seq2, match_score=3, mismatch_score=-3, gap_penalty=-2):
    # Initialize the scoring matrix
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = None

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback to get the alignment
    aligned_seq1, aligned_seq2 = [], []
    i, j = max_pos
    while score_matrix[i][j] != 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), max_score

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
        print("""
Usage: 
    python smith_waterman.py <sequence1> <sequence2>    
    
    Arguments:
        sequence1:
            The first sequence to align or a file containing two sequences separated by a newline or two fasta formatted sequences.
            If a file is provided, the sequences will be read from the file and the first two sequences found will be aligned.
        sequence2: 
            The second sequence to align. 
            Ignored if sequence1 is a file containing at leaset two sequences separated by a newline or fasta formatted sequences.
    Options:
        -h, --help: Display this help message
        --match_score=<number>: The score for a match (default: 3)
        --mismatch_score=<number>: The score for a mismatch (default: -3)
        --gap_penalty=<number>: The penalty for a gap (default: -2)
        --characters_per_line=<number>: The number of characters to display per line (default: 60) for the alignment
""")
        sys.exit(1)

    def pretty_format(seq, characters_per_line=60, gaps=" "):
        out = ""
        for i in range(0, len(seq), characters_per_line):
            for j in range(i, min(i+characters_per_line, len(seq)), 10):
                out += seq[j:j+10] + gaps
            out += "\n"
        return out
    def parse_fasta(file):
        file_content = "".join(open(file).readlines())
        if ">" in file_content:
            seqs = ""
            for seq in file_content.split(">"):
                if seq != "":
                    seq = "".join( seq.split("\n")[1:] )
                    seqs += seq + "\n"
        else:
            seqs = file_content
        seqs = [seq for seq in seqs.split("\n") if seq != ""]
        return seqs

    from os.path import isfile

    seq2 = None
    if isfile(sys.argv[1]):
        seqs = parse_fasta(sys.argv[1])
        seq1 = seqs[0]
    else:
        seq1 = sys.argv[1]
        seqs = []
    if len(seqs) > 1:
        seq2 = seqs[1]
    else:
        if len(sys.argv) < 3:
            print("Error: Only one sequence found in file. Please provide a second sequence. Use -h or --help for more information.")
            sys.exit(1)
        seq2 = parse_fasta(sys.argv[2])[0] if isfile(sys.argv[2]) else sys.argv[2]

    options = {option.split("=")[0]: option.split("=")[1] for option in sys.argv if option.startswith("--") and "=" in option}
    match_score = int(options.get("--match_score", 3))
    mismatch_score = int(options.get("--mismatch_score", -3))
    gap_penalty = int(options.get("--gap_penalty", -2))
    characters_per_line = int(options.get("--characters_per_line", 60))
    
    aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2, match_score=match_score, mismatch_score=mismatch_score, gap_penalty=gap_penalty)
    
    print(f"Alignment score: {score}")
    print(f"Seq 1:")
    print("\t", pretty_format(seq1), sep="")
    print(f"Seq 2: {seq2}")
    print("\t", pretty_format(seq2), sep="")
    print(f"Alignment:")
    aligned_seq1 = pretty_format(aligned_seq1).split('\n')
    aligned_seq2 = pretty_format(aligned_seq2).split('\n')
    for aln_seq1, aln_seq2 in zip(aligned_seq1, aligned_seq2):
        print("\t", aln_seq1, sep="")
        print("\t", aln_seq2, sep="")


