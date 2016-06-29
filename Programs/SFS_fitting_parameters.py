#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Parameter file for the SFS_fitting.py program

"""

# Integer, Number of bins to ignore at the beginning of the SFS (to keep all bin values, clear=0, to ignore singletons, clear=1, and so on)
clear = 1

# Integer, Number of replicates for the sfs simulations by the simSFS program (for clean results, rep=10000000, make it smaller for faster and noisier results)
rep=1000

# Float, Starting point for t optimization
range_start=0.8

# Float, Ending value for t optimization
range_end=3.0

# Float, Step value for t optimization
step=0.01

# String, path to the folder where output files will be
path="../Output_files/"

# String; path to the file containing the observed SFS (must be folded, see the example for format)
sfs_file="../Data/FoldSFS_YRI.txt"   

# Integer, sample size of the observed SFS (number of sequences, 2n for diploids)
sample_size=216

# String, path to the folder where the executable program simSFS is
simSFS_path="simulSFS/"
