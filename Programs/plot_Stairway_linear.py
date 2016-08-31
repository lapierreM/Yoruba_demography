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
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties

rcParams.update({'font.size': 18})

gris=(127/255.0,127/255.0,127/255.0)
violet=(110/255.0,16/255.0,144/255.0)
vert=(110/255.0,192/255.0,56/255.0)
jaune=(226/255.0,184/255.0,1/255.0)
bleu=(54/255.0,125/255.0,162/255.0)
rouge=(207/255.0,35/255.0,43/255.0)

def plot_Stairwayplot(file, facteur, graphique, linear, couleur):
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
	
	mutation_per_site_rescaled=[x*facteur for x in mutation_per_site]
	theta_per_site_rescaled=[x*facteur for x in theta_per_site_median]
	theta_per_site_2_5_rescaled=[x*facteur for x in theta_per_site_2_5]
	theta_per_site_97_5_rescaled=[x*facteur for x in theta_per_site_97_5]
	year_rescaled=[x*facteur for x in year]
	Ne_rescaled=[y*facteur for y in Ne_median]
	Ne_2_5_rescaled=[y*facteur for y in Ne_2_5]
	Ne__97_5_rescaled=[y*facteur for y in Ne_97_5]
		
	plt.figure(1)
	FontProperties().set_family('sans-serif')

	
	if graphique=="year_Ne":
		plt.grid(which='both',linestyle='-',linewidth=0.5,color='0.75')
		if linear==1:
			plt.plot([1,(2.48*2*Ne_rescaled[0]*24)],[Ne_rescaled[0],0],'--',color=vert, linewidth=2)
		plt.plot(year_rescaled,Ne_rescaled,color=couleur,linewidth=2,label=file.split('/')[-2])
	
		#plt.xlabel("year")
		#plt.ylabel("Ne")
		
################################################

#
#Read parameter values in file
#
with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
			
		if list[0][:4]=='path':
			path=re.sub(r'"','',(str(list[1][:-1]).strip()))
		
		elif list[0][:6]=='factor':
			factor=float(list[1])
		
		elif list[0][:5]=='graph':
			graph=re.sub(r'"','',(list[1][:-1].strip()))
			
			
summary_file=path+"Linear_1000_200files_summary"
plot_Stairwayplot(summary_file,factor, graph,0,jaune)
summary_file=path+"Linear_10000_200files_summary"
plot_Stairwayplot(summary_file,factor, graph,0,bleu)
summary_file=path+"Linear_100000_200files_summary"
plot_Stairwayplot(summary_file,factor, graph,1,violet)

#plt.xscale('log')

ax = plt.gca()

ax.tick_params(width=1, length=7,which='major')

ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()

xticks = ax.xaxis.get_major_ticks()
for i in range(1,9,2):
	xticks[i].label1.set_visible(False)

yticks = ax.yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)

plt.savefig('Fig_stairway_linear.eps',format='eps')
		
plt.show()

