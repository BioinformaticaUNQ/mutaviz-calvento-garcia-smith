from modeller import *

from backend.models.file_name_generator import FileNameGenerator


class Aligner:
    def __init__(self, path, sequence_name, sequence_1, pdb_key, sequence_2, pdb_file_path):
        self.__pdb_file_path = pdb_file_path
        self.__path = path
        self.__sequence_name = sequence_name
        self.__pdb_key = pdb_key
        self.__sequence_1 = sequence_1
        self.__sequence_2 = sequence_2

    def execute_alignment(self):
        align_file_path = self.__fasta_file()
        print("Align file: " + align_file_path)

        result_pir_file_path = self.__align(align_file_path)
        print(result_pir_file_path)

        return result_pir_file_path

    def __align(self, align_file_path):
        alignments_path = 'backend/alignments'
        pir_file_path = FileNameGenerator().random(extension='pir', path=alignments_path)
        self.__fasta_to_pir(align_file_path, pir_file_path)
        env = environ()
        aln = alignment(env)
        mdl = model(env, file=self.__pdb_file_path)
        aln.append_model(mdl, align_codes=self.__pdb_key, atom_files=self.__pdb_file_path)
        aln.append(file=pir_file_path, align_codes=self.__sequence_name)
        aln.align2d()
        result_pir_file_path = FileNameGenerator().random(extension='pir', path=alignments_path)
        aln.write(file=result_pir_file_path, alignment_format='PIR')

        return result_pir_file_path

    def __fasta_to_pir(self, align_file_path, pir_file_path):
        e = environ()
        a = alignment(e, file=align_file_path, align_codes='all', alignment_format='FASTA')
        a.write(file=pir_file_path, alignment_format='PIR')

    def __fasta_file(self):
        return self.__in_file(lambda file: self.__write_fasta(file))

    def __in_file(self, write_block):
        file_path = self.__random_file_name()
        file = open(file_path, 'w')
        write_block(file)
        file.close()
        return file_path

    def __write_fasta(self, file):
        file.write(">%s" % self.__sequence_name)
        file.write("\n")
        file.write(self.__sequence_1)

    def __random_file_name(self):
        return FileNameGenerator().random(extension='fasta', path=self.__path)
