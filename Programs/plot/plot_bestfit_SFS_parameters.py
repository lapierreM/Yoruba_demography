#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Parameter file for the plot_bestfit_SFS.py program

"""

# Integer, Number of bins to ignore at the beginning of the SFS (to keep all bin values, clear=0, to ignore singletons, clear=1, and so on)
clear = 1

# String; path to the file containing the observed SFS (must be folded, see the example for format)
sfs_file="../Data/FoldSFS_YRI.txt"

# String, path to the file where the output file .bestfit is
bestfit_file="../Output_files/example_FoldSFS_YRI.txt.bestfit"

# String, path to the file where the stairway plot SFS is
swp_file="../Output_files/example_swp_SFS_YRI"

# Integer, sample size of the observed SFS (number of sequences, 2n for diploids)
sample_size=216

# String, name of the saved figure
figure_name="Fig_bestfit_SFS.eps"
