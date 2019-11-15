from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    stop_codons = ['TAA', 'TGA', 'TAG']
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in stop_codons:
            break
        protein += table[codon]
    return protein


def read_seq(input_file):
    with open(input_file, "r") as f:
        lines = f.read().splitlines(True)
    seq = [line for line in lines if not line.startswith(">")]
    seq = "".join(seq)
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


def is_same_protein(alignment):
    # Hsp_identity / Hsp_align - len
    hsp = alignment.hsps[0]
    return round(hsp.identities / hsp.align_length, 2) == 1


if __name__ == "__main__":
    seq_string = read_seq("serum_albumin.fasta")
    sequence_data = translate(seq_string[110:])
    print(sequence_data)
    result_sequence = NCBIWWW.qblast(
        "blastp", "pdb", sequence_data, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
    )
    # Parametrizar word size, matrix_name
    print(result_sequence)
    blast_records = NCBIXML.read(result_sequence)
    exact_protein = next(alignment for alignment in blast_records.alignments if is_same_protein(alignment))
    if exact_protein:
        print("Esta es igual")

# ADN:
    # TGA
    # TAA
    # TAG
    #
    # ARN:
    # UAG
    # UAA
    # UGA

