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

import sys
import numpy
import math
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
		
		elif list[0][:3]=='rep':
			rep=int(list[1])
			
		elif list[0][:2]=='x0':
			x0=float(list[1])
			
		elif list[0][:2]=='dt':
			dt=float(list[1])
			
		elif list[0][:3]=='tau':
			tau=float(list[1])
	
		elif list[0][:4]=='path':
			path=re.sub(r'"','',(str(list[1][:-1]).strip()))
			
#
#Trajectory simulation
#

tps_fix_WF=[]
t_label=[]
traj_WF=[]

while len(tps_fix_WF)<nb_simul:
	x=x0
	tabx=[]
	t=0
	tps=[]
	
	while x>0 and x<1:
		tabx.append(x)
		x+=numpy.random.normal((1-x)*dt, math.sqrt(x*(1-x)*dt)) #Wright-Fisher diffusion
		t+=dt
		tps.append(t)
	
	if x>0:
		if t<=tau:
			traj_WF.append(tabx)
			tps_fix_WF.append(t)
			t_label.append(tps)


#
#Mean trajectory over time
#

for j in range(len(traj_WF)):
	traj_WF[j]=traj_WF[j]+[1]*(int(tau/dt)-len(traj_WF[j]))
	
avg_WF=[float(sum(col))/len(col) for col in zip(*traj_WF)]	

x_WF=[]

with open(path,"w") as file_WF:
	for i in xrange(len(avg_WF)):
		file_WF.write(str(i*dt)+'\t'+str(avg_WF[i])+'\n')
		x_WF.append(i*dt)
