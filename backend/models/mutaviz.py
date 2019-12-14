from functools import reduce
from io import StringIO
from os import mkdir
from os.path import isdir
from pathlib import Path
from shutil import copyfile, rmtree

from Bio.Blast import NCBIWWW, NCBIXML
from modeller import *
from pymol2 import PyMOL

from backend.models import logs
from backend.models.aligner import Aligner, AlignmentFormatter
from backend.models.file_name_generator import FileNameGenerator
from backend.models.modeller import Modeller
from backend.models.synth import Synthesizer


class Mutaviz:
    def __init__(self, seq_string, mutations, sequence_name):
        self.__mutations = mutations
        self.__seq_string = seq_string
        self.__sequence_name = sequence_name
        self.__protein_chain = None
        self.__most_similar_structure = None
        self.__mutated_sequence = None
        self.__original_pdb_filename = None
        self.__logger = logs.logger()

    def process(self):
        self.__add_required_dirs()
        self.__debug("Processing sequence " + self.__seq_string)
        self.__protein_chain = self.synthesize(self.__seq_string)
        self.__debug("Resulting chain " + self.__protein_chain)
        self.__blast()
        results = self.__process_and_model_blast_result()
        self.__print_results(results)
        self.__clean_files()

    def __print_results(self, results):
        pymol = PyMOL()
        pymol.start()
        pymol.cmd.load(results[0], "result_1")
        pymol.cmd.load(results[1], "result_2")
        pymol.cmd.png(self.__outputs_path() + "model_alignment.png")

    def __process_and_model_blast_result(self):
        if self.is_same_protein():
            self.__debug("Exact protein found!")
            self.mutate()
            mutated_protein = self.synthesize(self.__mutated_sequence)
            self.__debug("Mutated protein: " + mutated_protein)
            alignment_file = self.align(mutated_protein)
            output_pdb_file = self.__move_to_outputs(self.__original_pdb_filename, self.__pdb_key() + ".pdb")
            self.__debug("Alignment file: " + alignment_file)
            model_filename = self.model_structure(alignment_file)
            model_pdb_file = self.__move_to_outputs(model_filename,
                                                    self.__sequence_name + "_" + self.__pdb_key() + ".pdb")

            return [output_pdb_file, model_pdb_file]
        else:
            alignment_file = self.align(self.protein_chain)
            problem_sequence_structure_file = self.model_structure(alignment_file)
            output_pdb_file = self.__move_to_outputs(problem_sequence_structure_file, self.__pdb_key() + ".pdb")
            self.mutate()
            mutated_protein = self.synthesize(self.__mutated_sequence)
            alignment_file = self.align(mutated_protein)
            print(alignment_file)
            model_filename = self.model_structure(alignment_file)
            model_pdb_file = self.__move_to_outputs(model_filename,
                                                    self.__sequence_name + "_" + self.__pdb_key() + ".pdb")
            return [output_pdb_file, model_pdb_file]

    def __move_to_outputs(self, src, filename):
        return copyfile(src, self.__outputs_path() + filename)

    def __outputs_path(self):
        return str(Path(__file__).parent.parent) + '/outputs/'

    def __debug(self, message):
        return self.__logger.info(message)

    def model_structure(self, alignment_file):
        return 'backend/modeller/' + Modeller().execute(
            alignment_file=alignment_file, pdb_id=self.__pdb_key(), sequence=self.__sequence_name
        )['name']

    def synthesize(self, sequence):
        return Synthesizer.accepting(Synthesizer.ADN, sequence[0:]).run()

    def mutate(self):
        self.__mutated_sequence = self.__seq_string

        self.__debug("Mutating sequence")
        for (position, mutation) in self.__mutations.items():
            changing_seq = list(self.__mutated_sequence)
            changing_seq[position] = mutation
            self.__mutated_sequence = ''.join(changing_seq)

        self.__debug("Mutated sequence: " + self.__mutated_sequence)
        return self.__mutated_sequence

    def __blast(self):
        self.__debug("Performing blast")
        # scan_result = NCBIWWW.qblast(
        #     "blastp", "pdb", self.__protein_chain, word_size=2, threshold=200000, matrix_name="BLOSUM62", gapcosts="11 1"
        # )
        with open("backend/serum_albumin_result.xml", "r") as f:
            file = f.read()
        scan_result = StringIO(file)
        self.__debug("Blast query done")
        blast_records = NCBIXML.read(scan_result)
        self.__debug("Finding best match")
        self.__most_similar_structure = reduce(lambda result, output: self.__most_similar_between(result, output), blast_records.alignments)

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
        if self.identity_percentage(alignment) >= self.identity_percentage(another_alignment):
            return alignment
        return another_alignment

    def identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 1)

    def is_same_protein(self):
        # Hsp_identity / Hsp_align - len
        identity_percentage = self.identity_percentage(self.__most_similar_structure)
        self.__debug("Matching structure identity percentage " + str(identity_percentage))

        return identity_percentage == 1

    def __add_required_dirs(self):
        parent_path = str(Path(__file__).parent.parent)
        paths = ["/atom_files", "/modeller", "/alignments", "/outputs"]
        for path in paths:
            if not isdir(parent_path + path):
                mkdir(parent_path + path)

    def __clean_files(self):
        parent_path = str(Path(__file__).parent.parent)
        paths = ["/atom_files", "/modeller", "/alignments"]
        for path in paths:
            rmtree(parent_path + path)
