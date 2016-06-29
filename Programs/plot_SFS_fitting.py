#!/anaconda/bin/python
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
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
from matplotlib.font_manager import FontProperties

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

#
#Read parameter values in file
#
with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
			
		if list[0][:5]=='clear':
			clear=int(list[1])
		
		elif list[0][:8]=='fit_file':
			fit_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
		
		elif list[0][:11]=='figure_name':
			figure_name=re.sub(r'"','',(list[1][:-1].strip()))

t=[]
t_BD=[]
BD=[]
lK=[]
cK=[]
bnK=[]
eK=[]

#
# plot parameters and colors
#
minorLocator = MultipleLocator(0.1)

gris=(127/255.0,127/255.0,127/255.0)
violet=(110/255.0,16/255.0,144/255.0)
vert=(110/255.0,192/255.0,56/255.0)
jaune=(226/255.0,184/255.0,1/255.0)
bleu=(54/255.0,125/255.0,162/255.0)
rouge=(207/255.0,35/255.0,43/255.0)

#
# Read the .fit file
#
with open(fit_file,"r") as file_fit:
	for ligne in file_fit:
		liste=re.split(' |\n',ligne)
		if ligne[0]!='#' and len(liste)==7:
			t.append(liste[0])
			BD.append(liste[1])
			lK.append(liste[2])
			cK.append(liste[3])
			bnK.append(liste[4])
			eK.append(float(liste[5]))
			

t_ek=[1/float(x) for x in t] 
t_BD=[float(x)*2 for x in t] #time scaling 2 to put Birth-Death on the same time scale as Kingman models

rcParams.update({'font.size': 18})
FontProperties().set_family('sans-serif')

plt.figure(1)

plt.plot(t_BD,BD,color=rouge,label="Birth-death",linewidth=2)
plt.plot(t,lK,color=vert,label="Linear Kingman",linewidth=2)
plt.plot(t,cK,color=bleu,label="Conditionned Kingman",linewidth=2)
plt.plot(t,bnK,color=violet,label="Bottleneck Kingman",linewidth=2)
plt.plot(t_ek,eK,color=jaune,label="Exponential Kingman",linewidth=2)

plt.yscale('log')
#plt.xlim([0.8,3.0]) #if specific bounds for the x axis are needed

ax = plt.gca()
ax.xaxis.set_minor_locator(minorLocator)

ax.tick_params(width=1, length=7,which='major')
ax.tick_params(width=1, length=4,which='minor')

ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()

plt.savefig(figure_name,format='eps')

plt.show()
