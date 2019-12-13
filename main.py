import os

from backend.models.mutaviz import Mutaviz
from pymol2 import PyMOL


def read_seq(input_file):
    with open(input_file, "r") as f:
        lines = f.read().splitlines(True)
    seq = [line for line in lines if not line.startswith(">")]
    seq = "".join(seq)
    seq = seq.replace("\n", "")
    seq = seq.replace("\r", "")
    return seq


if __name__ == "__main__":
    os.environ["LM_LICENSE_FILE"] = "pymol-edu-license.lic"
    seq_string = read_seq("backend/serum_albumin_dna.fasta")
    muta = Mutaviz(seq_string[0:], {10: "A"}, "testing")
    result = muta.process()
    print(result)
    pymol = PyMOL()
    pymol.start()
    pymol.cmd.load(result[0], "el1")
    pymol.cmd.load(result[1], "el2")
    pymol.cmd.png("coso.png")
