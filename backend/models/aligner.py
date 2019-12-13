import os

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.PDB import PDBList, PDBParser

from backend.models.file_name_generator import FileNameGenerator

from modeller import *


class Aligner:
    def __init__(self, path=None, sequence_name=None, sequence_1=None, pdb_key=None, sequence_2=None):
        self.__path = path
        self.__sequence_name = sequence_name
        self.__pdb_key = pdb_key
        self.__sequence_1 = sequence_1
        self.__sequence_2 = sequence_2

    def file_align(self):
        self.__fetch_structure_info()
        fasta_path = self.__in_file(lambda file: self.write_fasta(file))

        return fasta_path

    def __in_file(self, write_block):
        file_path = self.random_file_name()
        file = open(file_path, 'w')
        write_block(file)
        file.close()
        return file_path

    def write_fasta(self, file):
        file.write(">%s" % self.__sequence_name)
        file.write("\n")
        file.write(self.__sequence_1)

    def random_file_name(self):
        return FileNameGenerator().random(extension='fasta', path=self.__path)

    def __fetch_structure_info(self):
        new_file_name = self.__fetch_pdb()

        structure = PDBParser().get_structure(self.__pdb_key, new_file_name)
        chain = list(structure.get_chains())[0]
        residues = list(chain.get_residues())
        atoms = list(filter(lambda residue: not residue.get_id()[0].strip(), residues))
        first_atom_residue = atoms[0].get_id()[1]
        last_atom_residue = atoms[-1].get_id()[1]

        self.__structure_info = {'chain': chain.get_id(), 'start': first_atom_residue, 'end': last_atom_residue}
        return self.__structure_info

    def __fetch_pdb(self):
        pdb_file_path = PDBList().retrieve_pdb_file(self.__pdb_key, pdir='backend/atom_files', file_format="pdb")
        new_file_name = 'backend/atom_files/%s.pdb' % self.__pdb_key
        # renombro para que tenga el nombre que espera modeller (xxxx.pdb)
        os.rename(pdb_file_path, new_file_name)
        self.__pdb_file_path = new_file_name
        return self.__pdb_file_path

    @property
    def pdb_file_path(self):
        return self.__pdb_file_path

    def structure_info(self):
        structure_info = {"seq_name": self.__sequence_name, "pdb_key": self.__pdb_key}
        structure_info.update(self.__structure_info)
        return structure_info


class AlignmentFormatter:
    def __init__(self, in_file=None, out_file=None, structure_info=None, path=None):
        self.__path = path
        self.__structure_info = structure_info
        self.__in_file = in_file
        self.__out_file = out_file

    def to_pir(self):
        e = environ()
        a = alignment(e, file=self.__in_file, align_codes='all', alignment_format='FASTA')
        a.write(file=self.__out_file, alignment_format='PIR')
        return self.__out_file
        # with open(self.__out_file, "r") as f:
        #     read = f.read()
        #     file_lines = read.strip().split(">")
        # sequences = list(map(lambda line: self.__get_sequence(line.split("\n")), file_lines))
        # sequences = list(filter(lambda seq: seq, sequences))
        # return self.__generate_pir_file(sequences)

    def __generate_pir_file(self, sequences):
        file_path = self.random_file_name()
        file = open(file_path, 'w')
        self.__write_pir(file, sequences)
        file.close()
        return file_path

    def __get_sequence(self, line):
        return "\n".join(list(filter(lambda line: not (line.startswith("sequence") or line.startswith("P1") or not line.strip()),
                           line)))

    def __write_pir(self, file, sequences):
        structure_info = self.__structure_info
        sequence_name = structure_info["seq_name"]
        file.write(">P1;%s" % sequence_name)
        file.write("\n")
        file.write("%s::::::::-1.00:-1.00" % (sequence_name))
        file.write("\n")
        file.write(sequences[0])
        file.write("\n")
        file.write("\n")
        pdb_key = structure_info["pdb_key"]
        file.write(">P1;%s" % pdb_key)
        file.write("\n")
        file.write(
            f"structureX:{pdb_key}:{structure_info['start']}:{structure_info['chain']}:{structure_info['end']}"
            f":{structure_info['chain']}:::-1.00:-1.00"
        )
        file.write("\n")
        file.write(sequences[1])

    def random_file_name(self):
        return FileNameGenerator().random(extension='pir', path=self.__path)

