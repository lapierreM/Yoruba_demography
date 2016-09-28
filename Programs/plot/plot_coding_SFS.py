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
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.font_manager import FontProperties

#Kingman_SFS: get the expected SFS under the standard Kingman coalescent
#arguments: sample size
#return: SFS (list)
def Kingman_SFS( n ):
	
	sfs = []
	
	an=0
	for k in xrange(1,n):
		an += 1.0/(k+0.0)

	for k in xrange(n-1):
		sfs.append( (1.0/(k+1.0))/an )
		
	return sfs

#plot_SFS: plot the transformed and normalized SFS
#arguments: file containing the SFS (see example for format), color, line type
#return: length of the plotted SFS (integer)
def plot_SFS(file,couleur,trait):
	with open(file,"r") as file1:
		sfs=[]
		for ligne in file1:
			sfs.append(int(ligne.split('\t')[1][:-1]))
	n=len(sfs)

	#Transform
	sfs_t=[]
	for i in xrange(1,n):
		sfs_t.append(float(sfs[i-1])*(i*(2*n-i)/(2*n+0.0)))
	sfs_t.append(float(sfs[n-1])*n)

	#Normalize
	sfs_t_n=[float(y)/sum(sfs_t) for y in sfs_t]

	x=[float(i)/(2*n) for i in range(1,n+1)]
	
	plt.figure(1)
	plt.plot(x,sfs_t_n,trait,color=couleur,linewidth=1.5,label=file.split('/')[-1])
	
	return n

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

#
#Read parameter values in file
#
with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
		
		if list[0][:8]=='CDS_file':
			CDS_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:11]=='notCDS_file':
			notCDS_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:8]=='all_file':
			all_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:11]=='figure_name':
			figure_name=re.sub(r'"','',(list[1][:-1].strip()))

#
#plot parameters and colors
#
minorLocator = MultipleLocator(0.02)
rcParams.update({'font.size': 18})

gris=(127/255.0,127/255.0,127/255.0)
violet=(110/255.0,16/255.0,144/255.0)
vert=(110/255.0,192/255.0,56/255.0)
bleu=(54/255.0,125/255.0,162/255.0)

#
# plot coding and non-coding SFS
#
n=plot_SFS(CDS_file,bleu,'-')
n=plot_SFS(notCDS_file,vert,'-')
n=plot_SFS(all_file,violet,'--')

#
# get Kingman SFS
#
sfs_kgm=Kingman_SFS(n+1)

# Transform Kingman SFS
sfs_kgm_t=[]
for i in xrange(1,n+1):
	sfs_kgm_t.append(float(sfs_kgm[i-1])*i)
sfs_kgm_t_n=[float(y)/sum(sfs_kgm_t) for y in sfs_kgm_t]

x=[float(i)/(2*n) for i in range(1,n+1)]

#
# plot Kingman SFS
#
plt.figure(1)
FontProperties().set_family('sans-serif')

plt.plot(x,sfs_kgm_t_n,'--',color=gris,label="Kingman",linewidth=2)

plt.xlim([0,0.51])

ax = plt.gca()
ax.xaxis.set_minor_locator(minorLocator)

ax.tick_params(width=1, length=7,which='major')
ax.tick_params(width=1, length=4,which='minor')

ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()

xticks = ax.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)

plt.savefig(figure_name,format='eps')

plt.show()