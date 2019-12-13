from functools import reduce
from io import StringIO

from Bio.Blast import NCBIWWW, NCBIXML

from backend.models.aligner import Aligner, AlignmentFormatter
from backend.models.file_name_generator import FileNameGenerator
from backend.models.modeller import Modeller
from backend.models.synth import Synthesizer
from modeller import *


class Mutaviz:
    def __init__(self, seq_string, mutations, sequence_name):
        self.__mutations = mutations
        self.__seq_string = seq_string
        self.__sequence_name = sequence_name
        self.__protein_chain = None
        self.__most_similar_structure = None
        self.__mutated_sequence = None
        self.__original_pdb_filename = None

    def process(self):
        print("Processing chain")
        self.__protein_chain = self.synthesize(self.__seq_string)
        print("Chain " + self.__protein_chain)
        self.blast()
        if self.is_same_protein():
            self.mutate()
            mutated_protein = self.synthesize(self.__mutated_sequence)
            alignment_file = self.align(mutated_protein)
            print(alignment_file)
            model_filename = self.model_structure(alignment_file)
            return [model_filename, self.__original_pdb_filename]
        else:
            alignment_file = self.align(self.protein_chain)
            original_model_filename = self.model_structure(alignment_file)
            self.mutate()
            mutated_protein = self.synthesize(self.__mutated_sequence)
            alignment_file = self.align(mutated_protein)
            print(alignment_file)
            model_filename = self.model_structure(alignment_file)
            return [model_filename, original_model_filename]

    def model_structure(self, alignment_file):
        return 'backend/modeller/' + Modeller().execute(
            alignment_file=alignment_file, pdb_id=self.__pdb_key(), sequence=self.__sequence_name
        )['name']

    def synthesize(self, sequence):
        return Synthesizer.accepting(Synthesizer.ADN, sequence[0:]).run()

    def mutate(self):
        self.__mutated_sequence = self.__seq_string

        for (position, mutation) in self.__mutations.items():
            changing_seq = list(self.__mutated_sequence)
            changing_seq[position] = mutation
            self.__mutated_sequence = ''.join(changing_seq)

        return self.__mutated_sequence

    def blast(self):
        print("Blast")
        # scan_result = NCBIWWW.qblast(
        #     "blastp", "pdb", self.__protein_chain, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
        # )
        with open("backend/Z7ZU2J3C016-Alignment.xml", "r") as f:
            file = f.read()
        scan_result = StringIO(file)
        print("Blast query done")
        blast_records = NCBIXML.read(scan_result)
        print("Finding best match")
        self.__most_similar_structure = reduce(lambda result, output: self.__most_similar_between(result, output), blast_records.alignments)
        # self.__most_similar_structure = blast_records.alignments[1]

    def __pdb_key(self):
        return self.__most_similar_structure.accession.split("_")[0]

    def __matching_sequence(self):
        return self.__most_similar_structure.hsps[0].sbjct

    def align(self, protein):

        aligner = Aligner(path="backend/alignments", sequence_name=self.__sequence_name, sequence_1=protein,
                          pdb_key=self.__pdb_key(), sequence_2=self.__matching_sequence())
        align_file_path = aligner.file_align()
        print(align_file_path)
        self.__original_pdb_filename = aligner.pdb_file_path

        pir_file_path = FileNameGenerator().random(extension='pir', path='backend/alignments')
        result_pir_file_path = FileNameGenerator().random(extension='pir', path='backend/alignments')
        AlignmentFormatter(align_file_path, pir_file_path, aligner.structure_info(), "backend/alignments").to_pir()

        env = environ()
        aln = alignment(env)
        mdl = model(env, file=self.__original_pdb_filename)
        aln.append_model(mdl, align_codes=self.__pdb_key(), atom_files=self.__original_pdb_filename)
        aln.append(file=pir_file_path, align_codes=self.__sequence_name)
        aln.align2d()
        aln.write(file=result_pir_file_path, alignment_format='PIR')
        print(result_pir_file_path)
        return result_pir_file_path

    @property
    def protein_chain(self):
        return self.__protein_chain

    @property
    def original_sequence(self):
        return self.__seq_string

    def __most_similar_between(self, alignment, another_alignment):
        if self.identity_percentage(alignment) > self.identity_percentage(another_alignment):
            return alignment
        return another_alignment

    def identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 2)

    def is_same_protein(self):
        # Hsp_identity / Hsp_align - len
        print(self.identity_percentage(self.__most_similar_structure))
        return self.identity_percentage(self.__most_similar_structure) == 1
