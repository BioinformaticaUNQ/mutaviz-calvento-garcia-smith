import argparse
import json

from mutaviz.models.mutaviz import Mutaviz

# this function is for extract the sequence from the fasta file
def read_seq(input_file):
    with open(input_file, "r") as f:
        lines = f.read().splitlines(True)
    seq = [line for line in lines if not line.startswith(">")]
    seq = "".join(seq)
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


# this function is for reading the mutations from the mutations file
def read_mutations(mutations_file):
    with open(mutations_file, "r") as f:
        json_data = json.load(f)
        return {int(k): str(v) for k, v in json_data.items()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', help='Path of the Fasta file containing the problem sequence')
    parser.add_argument(
        '--gap-costs', help="BLAST gap costs, first existence and then extension. Default: '11 1'", default="11 1"
    )
    parser.add_argument('--matrix-name', help='BLAST matrix name. Default: BLOSUM62', default="BLOSUM62")
    parser.add_argument(
        '--mutations', help='Path of the mutations file, format must be a json with index and mutation i.e. {10: "A"}'
    )
    parser.add_argument('--name', help='Name of the program run')
    parser.add_argument('--open-pymol', help='If true, opens PyMOL with both PDB files. Default false', default="")
    parser.add_argument('--seq-end', help='Ending position of the given sequence. You can use GenBank info')
    parser.add_argument(
        '--seq-start', help='Starting position of the given sequence (starts at 1). You can use GenBank info'
    )
    parser.add_argument(
        '--seq-type', help='Type of the fasta sequence (DNA, RNA or PROTEIN). Default: DNA', default="DNA"
    )
    parser.add_argument('--threshold', help='BLAST threshold. Default: 10', default=10)
    parser.add_argument('--word-size', help='BLAST word size. Default: 6', default=6)
    parser.add_argument('--output-path', help='Path where the output files will be stored', default='./outputs')
    args = parser.parse_args()

    seq_string = read_seq(args.fasta)
    mutations = read_mutations(args.mutations)
    start = args.seq_start and int(args.seq_start) - 1 or 0
    end = args.seq_end and int(args.seq_end) or None
    open_pymol = args.open_pymol == 'true'

    mutaviz = Mutaviz(
        seq_string=seq_string[start:end],
        mutations=mutations,
        seq_name=args.name,
        seq_type=args.seq_type,
        output_path=args.output_path)

    mutaviz.process(
        word_size=int(args.word_size), threshold=int(args.threshold),
        matrix_name=args.matrix_name, gap_costs=args.gap_costs, open_pymol=open_pymol
    )
