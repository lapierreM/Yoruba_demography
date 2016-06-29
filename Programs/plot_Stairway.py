#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
	Copyright (C) 2016  Lapierre/Marguerite

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

	for more information, please contact: marguerite.lapierre@mnhn.fr
"""

import re
import sys
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def plot_Stairwayplot(file, factor, graph):
	mutation_per_site=[]
	theta=[]
	theta_per_site_median=[]
	theta_per_site_2_5=[]
	theta_per_site_97_5=[]
	year=[]
	Ne_median=[]
	Ne_2_5=[]
	Ne_97_5=[]

	with open(file,"r") as file1:
		read=0
		for ligne in file1:
			liste=ligne.split('\t')
			if liste[0]=="mutation_per_site":
				read=1
			
			elif len(liste)>1 and read==1:
				mutation_per_site.append(float(liste[0]))
				theta.append(int(liste[1]))
				theta_per_site_median.append(float(liste[2]))
				theta_per_site_2_5.append(float(liste[3]))
				theta_per_site_97_5.append(float(liste[4]))
				year.append(float(liste[5]))
				Ne_median.append(float(liste[6]))
				Ne_2_5.append(float(liste[7]))
				Ne_97_5.append(float(liste[8]))
	
	mutation_per_site_rescaled=[x*factor for x in mutation_per_site]
	theta_per_site_rescaled=[x*factor for x in theta_per_site_median]
	theta_per_site_2_5_rescaled=[x*factor for x in theta_per_site_2_5]
	theta_per_site_97_5_rescaled=[x*factor for x in theta_per_site_97_5]
	year_rescaled=[x*factor for x in year]
	Ne_rescaled=[y*factor for y in Ne_median]
	Ne_2_5_rescaled=[y*factor for y in Ne_2_5]
	Ne__97_5_rescaled=[y*factor for y in Ne_97_5]
		
	#
	# plot parameters and colors
	#
	minorLocator = MultipleLocator(10000)
	gris=(127/255.0,127/255.0,127/255.0)

	plt.figure(1)
	
	if graph=="year_Ne":
		plt.grid(which='both',linestyle='-',linewidth=0.5,color='0.75')
		
		plt.plot(year_rescaled,Ne_rescaled,color=gris,linewidth=2,label=file.split('/')[-1])
		plt.plot(year_rescaled,Ne_2_5_rescaled,'k')
		plt.plot(year_rescaled,Ne__97_5_rescaled,'k')

	
		#plt.xlim(0,700000) # if x axis scale needs to be precised
	
		plt.xlabel("year")
		plt.ylabel("Ne")
		
		ax = plt.gca()

		ax.tick_params(width=1, length=7,which='major')

		ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
		ax.get_yaxis().tick_left()

		
		#plt.legend()
		plt.savefig('Fig_stairway.eps',format='eps')

	elif graph=="mut_theta":

		plt.plot(mutation_per_site_rescaled,theta_per_site_rescaled,'k',linewidth=2,label=file.split('/')[-1])
		plt.plot(mutation_per_site_rescaled,theta_per_site_2_5_rescaled,'k')
		plt.plot(mutation_per_site_rescaled,theta_per_site_97_5_rescaled,'k')
	
		ax = plt.gca()
		ax.xaxis.set_minor_locator(minorLocator)

		ax.tick_params(width=1, length=7,which='major')
		ax.tick_params(width=1, length=4,which='minor')

		ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
		ax.get_yaxis().tick_left()

		plt.grid(which='both',linestyle='-',linewidth=0.5,color='0.75')
		#plt.legend()
	
#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

#
#Read parameter values in file
#
with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
			
		if list[0][:13]=='stairway_file':
			stairway_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
		
		elif list[0][:8]=='factor':
			factor=float(list[1])
		
		elif list[0][:11]=='graph':
			graph=re.sub(r'"','',(list[1][:-1].strip()))

	
plot_Stairwayplot(stairway_file,factor,graph)

plt.show()

