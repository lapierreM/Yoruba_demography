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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
		
		elif list[0][:8]=='sfs_file':
			sfs_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:12]=='bestfit_file':
			bestfit_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:8]=='swp_file':
			swp_file=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
		elif list[0][:11]=='sample_size':
			n=int(list[1])
			
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
jaune=(226/255.0,184/255.0,1/255.0)
bleu=(54/255.0,125/255.0,162/255.0)
rouge=(207/255.0,35/255.0,43/255.0)
marron=(147/255.0,73/255.0,21/255.0)

#
# get observed SFS
#
sfs=[]
with open(sfs_file,"r") as filedata:
	for ligne in filedata:
		sfs.append(int(ligne.split('\t')[1][:-1]))

sfs_norm=[float(x)/(sum(sfs)-sfs[0]) for x in sfs]

kgm=[]
BD=[]
lK=[]
cK=[]
bnK=[]
eK=[]

with open(bestfit_file,"r") as file_bestfit:
	for ligne in file_bestfit:
		liste=re.split(' |\n',ligne)
		if ligne[0]!='i' and len(liste)==9:
			kgm.append(liste[2])
			BD.append(liste[3])
			lK.append(liste[4])
			cK.append(liste[5])
			bnK.append(liste[6])
			eK.append(liste[7])

#
#Transform
#
for i in xrange((1+clear),n/2):
	
	kgm[i-1-clear]=float(kgm[i-1-clear])*((i)*(n-i)/(n+0.0))
	BD[i-1-clear]=float(BD[i-1-clear])*((i)*(n-i)/(n+0.0))
	lK[i-1-clear]=float(lK[i-1-clear])*((i)*(n-i)/(n+0.0))
	cK[i-1-clear]=float(cK[i-1-clear])*((i)*(n-i)/(n+0.0))
	bnK[i-1-clear]=float(bnK[i-1-clear])*((i)*(n-i)/(n+0.0))
	eK[i-1-clear]=float(eK[i-1-clear])*((i)*(n-i)/(n+0.0))

for i in xrange(1,n/2):
	sfs_norm[i-1]=float(sfs_norm[i-1])*((i)*(n-i)/(n+0.0))
	
sfs_norm[(n/2)-1]=float(sfs_norm[(n/2)-1])*(n/2)
	
kgm[(n/2)-1-clear]=float(kgm[(n/2)-1-clear])*(n/2)
BD[(n/2)-1-clear]=float(BD[(n/2)-1-clear])*(n/2)
lK[(n/2)-1-clear]=float(lK[(n/2)-1-clear])*(n/2)
cK[(n/2)-1-clear]=float(cK[(n/2)-1-clear])*(n/2)
bnK[(n/2)-1-clear]=float(bnK[(n/2)-1-clear])*(n/2)
eK[(n/2)-1-clear]=float(eK[(n/2)-1-clear])*(n/2)

#
#Normalize
#
sfs_norm=[float(y)/(sum(sfs_norm)-sfs_norm[0]) for y in sfs_norm]
kgm=[float(y)/sum(kgm) for y in kgm]
BD=[float(y)/sum(BD) for y in BD]
lK=[float(y)/sum(lK) for y in lK]
cK=[float(y)/sum(cK) for y in cK]
bnK=[float(y)/sum(bnK) for y in bnK]
eK=[float(y)/sum(eK) for y in eK]

#
#Get stairway plot SFS
#

swp=[]

with open(swp_file,"r") as file_swp:
	for ligne in file_swp:
		swp.append(float(ligne[:-1]))

#Fold
swp_plie=[]
for i in xrange(n-1):
	swp_plie.append(swp[i]+swp[2*n-i-2])
swp_plie.append(swp[n-1])

swp_plie=swp_plie[clear:]

#transform

for i in xrange(1+clear,n):
	swp_plie[i-1-clear]=float(swp_plie[i-1-clear])*((i)*(2*n-i)/(2*n+0.0))

swp_plie[n-1-clear]=float(swp_plie[n-1-clear])*(n)

#normalize
swp_norm=[float(y)/(sum(swp_plie)) for y in swp_plie]

plt.figure(1)

FontProperties().set_family('sans-serif')

x=[float(i)/(n) for i in range(1+clear,(n/2)+1)]
x_data=[float(i)/(n) for i in range(1,(n/2)+1)]

plt.plot(x_data,sfs_norm,'o',label="Observed",mfc='none')
plt.plot(x,kgm,'--',color=gris,label="Kingman",linewidth=2)
plt.plot(x,BD,color=rouge,label="Birth-death",linewidth=2.5)
plt.plot(x,lK,color=vert,label="Linear Kingman",linewidth=1.5)
plt.plot(x,cK,color=bleu,label="Conditionned Kingman",linewidth=1.5)
plt.plot(x,bnK,color=violet,label="Bottleneck Kingman",linewidth=1.5)
plt.plot(x,eK,color=jaune,label="Exponential Kingman",linewidth=1.5)
plt.plot(x,swp_norm,'.',color=marron,label="Stairway plot",markersize=6)

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
