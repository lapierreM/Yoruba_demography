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

import numpy,math
from matplotlib import pyplot as plt
from numpy import random
import re

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

#
#Read parameter values in file
#

with open(sys.argv[1],"r") as file_params:
	for line in file_params:
		list=line.split('=')
			
		if list[0][:8]=='nb_simul':
			nb_simul=int(list[1])
		
		elif list[0][:1]=='n':
			n=int(list[1])
			
		elif list[0][:2]=='x0':
			x0=float(list[1])
			
		elif list[0][:2]=='dt':
			dt=float(list[1])
			
		elif list[0][:3]=='tau':
			tau=float(list[1])
	
		elif list[0][:4]=='path':
			path=re.sub(r'"','',(str(list[1][:-1]).strip()))

traj_BD=[]
t_label=[]

while len(traj_BD)<nb_simul:
	
	#Initialization
	x=x0
	tabx=[]
	t=0
	
	while x>0 and t<tau:
		tabx.append(x)
		x+=numpy.random.normal(2*dt, math.sqrt(2*x*dt))
		t+=dt
	
	gamma=random.gamma(n,1.0/n)

	exp=-math.log(random.random())/n
	
	if x>0 and gamma<=x<=(gamma+exp):
		traj_BD.append(tabx)
			

#
#Mean trajectory over time
#

avg_BD=[float(sum(col))/len(col) for col in zip(*traj_BD)]

with open(path,"w") as file_BD:
	for i in xrange(len(avg_BD)):
		file_BD.write(str(i*dt)+'\t'+str(avg_BD[i])+'\n')
