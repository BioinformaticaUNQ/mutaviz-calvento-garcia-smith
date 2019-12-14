# mutaviz

## Setup

### Create conda env
`conda env create -f environment.yml`

### Update conda env
`conda env update -n mutaviz --file environment.yml  --prune`

### Activate conda env
`conda activate mutaviz`

### Deactivate conda env
`conda deactivate`

### Set modeller key in
`<where is conda installed>/envs/mutaviz/lib/modeller-9.23/modlib/modeller/config.py`

## Usage

### Bash command line
For options run `python main.py --help`, you will see usage help

```
usage: main.py [-h] [--fasta FASTA] [--gap-costs GAP_COSTS]
               [--matrix-name MATRIX_NAME] [--mutations MUTATIONS]
               [--name NAME] [--open-pymol OPEN_PYMOL] [--seq-end SEQ_END]
               [--seq-start SEQ_START] [--seq-type SEQ_TYPE]
               [--threshold THRESHOLD] [--word-size WORD_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  --fasta FASTA         Path of the Fasta file containing the problem sequence
  --gap-costs GAP_COSTS
                        BLAST gap costs, first existence and then extension.
                        Default: '11 1'
  --matrix-name MATRIX_NAME
                        BLAST matrix name. Default: BLOSUM62
  --mutations MUTATIONS
                        Path of the mutations file, format must be a json with
                        index and mutation i.e. {10: "A"}
  --name NAME           Name of the program run
  --open-pymol OPEN_PYMOL
                        If true, opens PyMOL with both PDB files. Default
                        false
  --seq-end SEQ_END     Ending position of the given sequence. You can use
                        GenBank info
  --seq-start SEQ_START
                        Starting position of the given sequence (starts at 1).
                        You can use GenBank info
  --seq-type SEQ_TYPE   Type of the fasta sequence (DNA, RNA or PROTEIN).
                        Default: DNA
  --threshold THRESHOLD
                        BLAST threshold. Default: 10
  --word-size WORD_SIZE
                        BLAST word size. Default: 6
```

#### Usage example
`python main.py --fasta=mutaviz/human_serum_albumin_dna.fasta --seq-start=111 --seq-end=1871 --mutations=mutations.json --name='testing' --gap-costs='11 1' --seq-type=DNA --open-pymol=true`

### Importing it as a Python module
```python
from mutaviz.models.mutaviz import Mutaviz

seq = ... # A String containing dna, arn or protein chain
mutations_dict = ... # A dict containing mutations for mutating the original sequence. i.e. {50: 'C'},
                     # where 50 is the position in the sequence and 'C' is the atom you want to use.
mutaviz = Mutaviz(seq_string=seq, mutations=mutations_dict, seq_name="testing", seq_type="DNA")
mutaviz.process(word_size=6, threshold=10, matrix_name="BLOSUM62", gap_costs="11 1", open_pymol=False)
```