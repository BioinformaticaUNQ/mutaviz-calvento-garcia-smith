import os
from pathlib import Path

from modeller import *
from modeller.automodel import *
from functools import reduce


class Modeller:
    def execute(self, alignment_file='aln_rat_3v03.pir', pdb_id='3V03', sequence='NM_134326'):
        generated_model = self.__generate_model(alignment_file, pdb_id, sequence)

        return reduce(lambda result, output: self.__best_between(result, output), generated_model.chains[0].seq.outputs)

    def __generate_model(self, alignment_file, pdb_id, sequence):
        os.chdir(str(Path(__file__).parent.parent) + '/modeller')
        generated_model = automodel(self.__environment(), alnfile=alignment_file, knowns=pdb_id, sequence=sequence,
                                    assess_methods=assess.DOPE)
        # code of the target
        generated_model.starting_model = 1  # index of the first model
        generated_model.ending_model = 5  # index of the last model
        # (determines how many models to calculate)
        generated_model.make()  # do the actual homology modeling
        os.chdir(str(Path(__file__).parent.parent.parent))
        return generated_model

    def __environment(self):
        env = environ()  # create a new MODELLER environment to build this model in
        # directories for input atom files
        env.io.atom_files_directory = str(Path(__file__).parent.parent.parent)
        env.io.hetatm = True

        return env

    def __best_between(self, x, y):
        if x['DOPE score'] > y['DOPE score']:
            return y
        return x
