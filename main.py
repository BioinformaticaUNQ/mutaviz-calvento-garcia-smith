from backend.models.mutaviz import Mutaviz


def read_seq(input_file):
    with open(input_file, "r") as f:
        lines = f.read().splitlines(True)
    seq = [line for line in lines if not line.startswith(">")]
    seq = "".join(seq)
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


if __name__ == "__main__":
    seq_string = read_seq("backend/serum_albumin_dna.fasta")
    muta = Mutaviz(seq_string[110:1871], {10: "G"}, "testing")
    result = muta.process()
    print(result)