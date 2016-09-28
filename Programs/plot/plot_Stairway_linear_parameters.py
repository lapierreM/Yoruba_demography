#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Parameter file for the plot_Stairway.py program

"""

# String, path to the folder containing the summary files of the Stairway plot
path="../Data/"

# Float, Stairway inference "factor" (if genome length was divided by 10 for the stairway inference to reduce computation time, factor=0.1)
factor=0.1

# String, type of plot needed ("year_Ne" for Ne depending on years, or "mut_theta" for theta per site depending on mutation per site)
graph="year_Ne"