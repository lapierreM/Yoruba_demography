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
from matplotlib.font_manager import FontProperties

gris=(127/255.0,127/255.0,127/255.0)
violet=(110/255.0,16/255.0,144/255.0)
jaune=(226/255.0,184/255.0,1/255.0)
bleu=(54/255.0,125/255.0,162/255.0)

minorLocator = MultipleLocator(0.1)
rcParams.update({'font.size': 18})

#Kingman_SFS: get the expected SFS under the standard Kingman coalescent
#arguments: sample size
#return: SFS (list)
def Kingman_SFS( n ):
	
	sfs = []
	
	an=0;
	for k in xrange(1,n):
		an += 1.0/(k+0.0)

	for k in xrange(n-1):
		sfs.append( (1.0/(k+1.0))/an )
		
	return sfs

#------------------------------------------------------------------------------------------
#MAIN
#------------------------------------------------------------------------------------------

c=0
couleurs=[jaune, bleu, violet]

#
#Read parameter values in file
#
with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
			
		if list[0][:8]=='lin_path':
			lin_path=re.sub(r'"','',(str(list[1][:-1]).strip()))

		elif list[0][:8]=='swp_path':
			swp_path=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
for rep in [1000,10000,100000]:

	mean_sfs_lin=[]

	with open(lin_path+"mean_linearSFS_"+str(rep),"r") as file_lin:
		for ligne in file_lin:
			mean_sfs_lin.append(float(ligne.split('\t')[1][:-1]))
	
	n=len(mean_sfs)+1
	
	sfs_swp=[]

	with open(swp_path+"example_swp_SFS_linear_"+str(rep),"r") as file_swp:
		for ligne in file_swp:
			sfs_swp.append(float(ligne[:-1]))
			
	#Transformation Kingman
	for i in xrange(1,n):
		mean_sfs_lin[i-1]=float(mean_sfs_lin[i-1])*i
		sfs_swp[i-1]=float(sfs_swp[i-1])*i
		
	#Normalization
	mean_sfs_lin_norm=[float(y)/(sum(mean_sfs_lin)) for y in mean_sfs_lin]
	sfs_swp_norm=[float(y)/(sum(sfs_swp)) for y in sfs_swp
	
	x=[float(i)/(n) for i in range(1,n)]

	plt.figure(1)
	plt.plot(x,mean_sfs_t_n,color=couleurs[c],label="Mean SFS, R = "+str(rep),linewidth=1.5)
	plt.plot(x,sfs_swp_norm,'.',color=couleurs[c],label="SFS stairway plot, R = "+str(rep),markersize=4,markevery=2)

	c+=1
	

#Kingman
sfs_kgm=Kingman_SFS(n+1)

#Transformation Kingman
sfs_kgm_t=[]
for i in xrange(1,n+1):
	sfs_kgm_t.append(float(sfs_kgm[i-1])*i)
sfs_kgm_t_n=[float(y)/sum(sfs_kgm_t) for y in sfs_kgm_t]

plt.plot(x,sfs_kgm_t_n,'--',color=gris,label="Kingman",linewidth=1.5)

plt.tick_params(axis='both', which='major', labelsize=18)

ax = plt.gca()
ax.xaxis.set_minor_locator(minorLocator)

ax.tick_params(width=1, length=7,which='major')
ax.tick_params(width=1, length=4,which='minor')

ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()

xticks = ax.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)
	
plt.savefig('Fig_SFS_linear.eps',format='eps')
plt.show()
