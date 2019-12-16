import os
from functools import reduce
from os import mkdir
from os.path import isdir
from pathlib import Path
from shutil import copyfile, rmtree

from Bio.Blast import NCBIXML, NCBIWWW
from Bio.PDB import PDBList
from pymol2 import PyMOL

from mutaviz.models import logs
from mutaviz.models.aligner import Aligner
from mutaviz.models.modeller import Modeller
from mutaviz.models.synth import Synthesizer


class Mutaviz:
    def __init__(self, seq, mutations, seq_name, seq_type, output_path='outputs'):
        """
        :param seq: Sequence you want to synthesize and model
        :type seq: str
        :param mutations: Mutations dict containing with array index as key and mutation as value
        :type mutations: dict
        :param seq_name: Name of the whole test, used for theoretical sequence alignment
        :type seq_name: str
        :param seq_type: Type of given sequence (DNA, RNA or PROTEIN)
        :type seq_type: str
        :param output_path: Outputs path for output files, default to 'outputs'
        :type output_path: str
        """
        self.__sequence_type = seq_type
        self.__mutations = mutations
        self.__sequence = seq
        self.__sequence_name = seq_name
        self.__protein_chain = None
        self.__most_similar_structure = None
        self.__mutated_sequence = None
        self.__logger = logs.logger()
        self.__outputs_path = str(Path(output_path).absolute()) + "/"

    def process(self, word_size=6, threshold=10, matrix_name="BLOSUM62", open_pymol=False, gap_costs="11 1"):
        """
        Process sequence string finding the matching structure, mutating it and then aligning both
        :param word_size: BLAST word size param, defaults to 6
        :type word_size: int
        :param threshold: BLAST threshold param, defaults to 10
        :type threshold: int
        :param matrix_name: BLAST matrix name, defaults to BLOSUM62
        :type matrix_name: str
        :param open_pymol: If True, opens PyMOL displaying both structures aligned. If False, prints a PNG with the
         structures aligned (you can find it with the output files), defaults to False
        :type open_pymol: bool
        :param gap_costs: BLAST gap costs param, defaults to 11 1 (existence and the extension)
        :type gap_costs: str
        """
        self.__add_required_dirs()
        self.__debug("Processing sequence " + self.__sequence)
        self.__protein_chain = self.synthesize(self.__sequence)
        self.__debug("Resulting chain " + self.__protein_chain)
        self.__blast(word_size, threshold, matrix_name, gap_costs)
        results = self.__process_and_model_blast_result()
        self.__print_results(open_pymol, results)
        self.__clean_files()

    def __print_results(self, open_pymol, results):
        if open_pymol:
            script_file = str(Path(__file__).parent) + "/pymol_script.py"
            cmd = f"pymol {script_file} {results[0]} {results[1]}"
            os.system(cmd)
        else:
            detached_pymol = PyMOL()
            detached_pymol.start()
            detached_pymol.cmd.load(results[0])
            detached_pymol.cmd.load(results[1])
            filename = results[0].split("/")[-1].split(".")[0]
            second_filename = results[1].split("/")[-1].split(".")[0]
            detached_pymol.cmd.png(f"{self.__outputs_path + filename}_{second_filename}_alignment.png")

    def __process_and_model_blast_result(self):
        if self.__is_same_protein():
            self.__debug("Exact protein found!")

            self.mutate()
            mutated_protein = self.synthesize(self.__mutated_sequence)
            self.__debug("Mutated protein: " + mutated_protein)

            pdb_file_path = self.__fetch_pdb()
            output_pdb_file = self.__move_to_outputs(
                pdb_file_path, f"{self.__sequence_name}_{self.__pdb_key}_selected_model.pdb"
            )

            alignment_file = self.__align(mutated_protein, self.__pdb_key, pdb_file_path)
            alignment_file = self.__move_to_outputs(
                alignment_file, f"{self.__sequence_name}_mutated_{self.__pdb_key}_alignment.pir"
            )
            self.__debug("Alignment file: " + alignment_file)

            model_pdb_file = self.__model_structure(
                alignment_file, self.__pdb_key, self.__sequence_name + "_mutation_theoretical_model"
            )

            self.__debug(f"Results: [{output_pdb_file}, {model_pdb_file}]")
            return [output_pdb_file, model_pdb_file]
        else:
            pdb_file_path = self.__fetch_pdb()
            self.__move_to_outputs(pdb_file_path, f"{self.__sequence_name}_{self.__pdb_key}_selected_model.pdb")

            alignment_file = self.__align(self.__protein_chain, self.__pdb_key, pdb_file_path)
            alignment_file = self.__move_to_outputs(
                alignment_file, f"{self.__sequence_name}_model_{self.__pdb_key}_alignment.pir"
            )

            modeled_pdb_name = self.__sequence_name + "_theoretical_model"
            modeled_pdb_file = self.__model_structure(alignment_file, self.__pdb_key, modeled_pdb_name)

            mutation_modeled_pdb_file = self.__mutate_and_model(modeled_pdb_name, modeled_pdb_file)
            self.__debug(f"Results: [{modeled_pdb_file}, {mutation_modeled_pdb_file}]")
            return [modeled_pdb_file, mutation_modeled_pdb_file]

    def __mutate_and_model(self, model_pdb_id, modeled_pdb_file):
        self.mutate()
        mutated_protein = self.synthesize(self.__mutated_sequence)
        self.__debug("Mutated protein: " + mutated_protein)

        alignment_file = self.__align(mutated_protein, model_pdb_id, modeled_pdb_file)
        alignment_file = self.__move_to_outputs(alignment_file, f"{self.__sequence_name}_mutated_alignment.pir")
        print(alignment_file)

        return self.__model_structure(
            alignment_file, model_pdb_id, self.__sequence_name + "_mutation_theoretical_model"
        )

    def __move_to_outputs(self, src, filename):
        return copyfile(src, self.__outputs_path + filename)

    def __debug(self, message):
        return self.__logger.info(message)

    def __model_structure(self, alignment_file, pdb_id, output_filename):
        model_filename = 'mutaviz/modeller/' + Modeller().execute(
            alignment_file=alignment_file, pdb_id=pdb_id, sequence=self.__sequence_name
        )['name']

        return self.__move_to_outputs(model_filename, output_filename + ".pdb")

    def synthesize(self, sequence):
        return Synthesizer.accepting(self.__sequence_type, sequence).run()

    def mutate(self):
        self.__mutated_sequence = self.__sequence

        self.__debug("Mutating sequence")
        for (position, mutation) in self.__mutations.items():
            changing_seq = list(self.__mutated_sequence)
            changing_seq[position] = mutation
            self.__mutated_sequence = ''.join(changing_seq)

        self.__debug("Mutated sequence: " + self.__mutated_sequence)
        return self.__mutated_sequence

    def __blast(self, word_size, threshold, matrix_name, gap_costs):
        # If you have a local blast result replace scan_result variable with these lines
        # with open("path/to/blast_query_result.xml", "r") as f:
        #     file = f.read()
        # scan_result = StringIO(file)

        self.__debug("Performing blast")
        scan_result = NCBIWWW.qblast(
            "blastp", "pdb", self.__protein_chain, word_size=word_size,
            threshold=threshold, matrix_name=matrix_name, gapcosts=gap_costs
        )
        self.__debug("Blast query done")
        blast_records = NCBIXML.read(scan_result)
        self.__debug("Finding best match")
        self.__most_similar_structure = reduce(
            lambda result, output: self.__most_similar_between(result, output), blast_records.alignments
        )

    @property
    def __pdb_key(self):
        return self.__most_similar_structure.accession.split("_")[0]

    def __matching_sequence(self):
        return self.__most_similar_structure.hsps[0].sbjct

    def __align(self, protein, pdb_key, pdb_file_path):
        aligner = Aligner(path="mutaviz/alignments", sequence_name=self.__sequence_name, sequence_1=protein,
                          pdb_key=pdb_key, sequence_2=self.__matching_sequence(), pdb_file_path=pdb_file_path)
        return aligner.execute_alignment()

    def __most_similar_between(self, alignment, another_alignment):
        if self.__identity_percentage(alignment) >= self.__identity_percentage(another_alignment) or\
                self.__evalue(alignment) <= self.__evalue(another_alignment):
            return alignment
        return another_alignment

    def __evalue(self, alignment):
        hsp = alignment.hsps[0]
        return hsp.evalue

    def __identity_percentage(self, alignment):
        hsp = alignment.hsps[0]
        return round(hsp.identities / hsp.align_length, 2)

    def __is_same_protein(self):
        # Hsp_identity / Hsp_align - len
        identity_percentage = self.__identity_percentage(self.__most_similar_structure)
        self.__debug("Matching structure identity percentage " + str(identity_percentage))

        return identity_percentage == 1

    def __add_required_dirs(self):
        parent_path = str(Path(__file__).parent.parent)
        paths = ["/atom_files", "/modeller", "/alignments"]
        for path in paths:
            if not isdir(parent_path + path):
                mkdir(parent_path + path)

    def __clean_files(self):
        self.__in_program_execution_dirs(lambda parent_path, path: rmtree(parent_path + path))

    def __in_program_execution_dirs(self, block):
        parent_path = str(Path(__file__).parent.parent)
        paths = ["/atom_files", "/modeller", "/alignments"]
        for path in paths:
            block(parent_path, path)

    def __fetch_pdb(self):
        pdb_file_path = PDBList().retrieve_pdb_file(self.__pdb_key, pdir='mutaviz/atom_files', file_format="pdb")
        new_file_name = 'mutaviz/atom_files/%s.pdb' % self.__pdb_key
        os.rename(pdb_file_path, new_file_name)
        return new_file_name
