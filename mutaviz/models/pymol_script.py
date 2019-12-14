import sys
import pymol

args = sys.argv
pymol.cmd.load(args[2])
pymol.cmd.load(args[3])
