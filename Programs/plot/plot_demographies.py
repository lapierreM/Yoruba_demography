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
import matplotlib.pyplot as plt
import math
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

#
#Read parameter values in file
#
with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
		
		if list[0][:11]=='models_file':
			models_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:2]=='BD':
			BD_trajectory_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
		
		elif list[0][:4]=='Cond':
			Cond_trajectory_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
				
		elif list[0][:11]=='figure_name':
			figure_name=re.sub(r'"','',(list[1][:-1].strip()))
			
t=[]
BD=[]
lK=[]
cK=[]
bnK=[]
eK=[]

#
# plot colors
#
gris=(127/255.0,127/255.0,127/255.0)
violet=(110/255.0,16/255.0,144/255.0)
vert=(110/255.0,192/255.0,56/255.0)
jaune=(226/255.0,184/255.0,1/255.0)
bleu=(54/255.0,125/255.0,162/255.0)
rouge=(207/255.0,35/255.0,43/255.0)

#
# Read the best fit t parameter value
#
with open(models_file,"r") as file_models:
	i=0
	for ligne in file_models:
		if i==1:
			liste=re.split(' |\t|\n',ligne)
			t_BD=liste[3][3:-1]
			t_lK=liste[5][3:-1]
			t_cK=liste[7][3:-1]
			t_bnK=liste[9][3:-1]
			t_eK=liste[11][3:-1]
		i+=1

#
# Read the implicit demography models (Birth-Death and Conditioned) trajectories
#
x_BD=[]
avg_BD=[]
with open(BD_trajectory_file,"r") as file_BD:
	for line in file_BD:
		liste=re.split('\t|\n',line)
		x_BD.append(float(liste[0]))
		avg_BD.append(float(liste[1]))

x_Cond=[]
avg_Cond=[]	
with open(Cond_trajectory_file,"r") as file_Cond:
	for line in file_Cond:
		liste=re.split('\t|\n',line)
		x_Cond.append(float(liste[0]))
		avg_Cond.append(float(liste[1]))
		
#
# Define the explicit demography models (Linear, Sudden and Exponential) trajectories
#
linear_x=[0,t_lK]
linear_y=[1,0]

sudden_x=[0,t_bnK,t_bnK]
sudden_y=[1,1,0]

expo_x=[x/100.0 for x in xrange(300)]
expo_y=[math.exp(-x*float(t_eK)) for x in expo_x]

#
# plot parameters
#
BD_x=[(float(t_BD)-x)*2 for x in x_BD]

Cond_x=[(float(t_cK)-x) for x in x_Cond]

minorLocator = MultipleLocator(0.1)
rcParams.update({'font.size': 18})
rcParams['axes.linewidth'] = 2

plt.figure(figsize=(20,8))

plt.plot(sudden_x,sudden_y,color=violet,label="Bottleneck Kingman",linewidth=3)
plt.plot(expo_x,expo_y,color=jaune,label="Exponential Kingman",linewidth=3)
plt.plot(linear_x,linear_y,color=vert,label="Linear Kingman",linewidth=3)

plt.plot(Cond_x,avg_Cond,'--',color=bleu,label="Conditionned Kingman",linewidth=3)
plt.plot(BD_x,avg_BD,'--',color=rouge,label="Birth-death",linewidth=3)

#plt.xlim([0,3]) if x axis limits needs to be defined
plt.ylim([0,1.03])

ax = plt.gca()
ax.xaxis.set_minor_locator(minorLocator)

ax.tick_params(width=1, length=15,which='major',direction='out')
ax.tick_params(width=1, length=12,which='minor',direction='out')

ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()

ax.xaxis.set_ticklabels([])
ax.yaxis.set_ticklabels([])

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.savefig(figure_name,format='eps')

plt.show()
