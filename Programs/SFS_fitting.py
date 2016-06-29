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

from math import *
import sys
import os
import re

#E_Xik: compute the expectancy of the value of index k of critical Birth-Death SFS (see Delaporte et al. 2016)
#arguments: index k of the SFS bin, sample size n, parameter tau (foundation time * sample size), sum
#return: expectancy of the value of index k of critical Birth-Death SFS (float)
def E_Xik( k , n, tau, my_hsum ):

	a = ( n - 3.0*k - 1.0 ) / k 
	b = ( n - k - 1.0 ) * ( k + 1.0 ) / ( k*tau )
	
	c = pow( 1.0 + 1.0/tau , k-1.0 ) / ( tau * tau )
	
	d = 2.0*tau*tau - (n-2.0*k-1.0)*2.0*tau - (n-k-1.0)*(k+1.0)
	
	e = log( 1.0+ tau ) 
	
	f = my_hsum[k-1]

	return a+b+c*d*(e-f);

#BD_SFS: compute the critical Birth-Death SFS (see Delaporte et al. 2016)
#arguments: sample size n, parameter tau (foundation time * sample size)
#return: critical Birth-Death SFS (list)
def BD_SFS( n , tau ):

	my_hsum = []
	sfs = []

	my_hsum.append( 0 )
	for i in xrange(1,n):
		my_hsum.append( my_hsum[i-1] + 1.0/( i*pow( (1.0 + 1.0/tau), i)) )
		
	for k in xrange(1,n):
		sfs.append( E_Xik( k, n, tau, my_hsum) )
	
	sfs=[float(y)/sum(sfs) for y in sfs]

	return sfs

#read_sfs_file: get the observed SFS from the sfs_file
#arguments: path to the SFS file (see example for format)
#return: observed SFS (list)
def read_sfs_file( infile ):

	obs_sfs = []

	with open( infile, 'r') as file_sfs:
		for ligne in file_sfs:
			liste = re.split('\t|\n',ligne)
			obs_sfs.append( float(liste[1]) )
	
	return obs_sfs

#norm_clean: normalize a SFS, with possibly omitting a given number of bins
#arguments: SFS to normalize (list), integer corresponding to the number of bins to omit (0 for none; 1 to omit singletons) for the normalization
#return: normalized SFS (list)
def norm_clean( sfs , start ):
	sfs_clean=sfs[start:]
	
	sfs_norm=[0]*len(sfs_clean)
	
	mysum=sum(sfs_clean)
		
	for k in xrange(len(sfs_clean)):
		sfs_norm[k] = sfs_clean[k]/mysum

	return sfs_norm

#sqr_dev: compute a least square distance between observed SFS and simulated SFS
#arguments: observed SFS (list), simulated SFS (list)
#return: distance (float)
def sqr_dev( sfs_obs, sfs_sim):

	d = 0
	for k in xrange( 0, len(sfs_obs) ):
		f1 = sfs_obs[k]
		f2 = sfs_sim[k]
	
		d += pow(f1-f2,2)/f2
		
	return d

#fold: fold a given SFS (should be done before cleaning the first bin)
#arguments: SFS (list) unfolded
#return: folded SFS (list)
def fold( sfs ):

	f_sfs = []

	l = len(sfs)
	
	for i in xrange( l/2 ):
		f_sfs.append( sfs[i] + sfs[l-1-i] )
	
	if l%2 == 1:
		f_sfs.append( sfs[ l/2 ] )

	return f_sfs

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


#condK: simulate (via simSFS program) SFS under the Conditioned model (Kingman that conditions all lineage to have coalesce before 't')
#arguments: sample size, parameter t for the conditionning, number of replicates for the SFS simulation, path to the simSFS program
#return: SFS (list) under the Conditioned model
def condK( n, t, rep, simSFS_path ):

	sfs = []
	ll = os.popen("\\"+simSFS_path+"simSFS c "+str(n)+" "+str(t)+" "+str(rep)).readlines()
	for l in ll:
		sfs.append( float(l[:-1]) )

	return sfs


#foundedK: simulate (via simSFS program) SFS under the Sudden model (Kingman that started at time t (it was 1 individual before that))
#arguments: sample size, parameter t for the demography, number of replicates for the SFS simulation, path to the simSFS program
#return: SFS (list) under the Sudden model
def foundedK( n, t, rep, simSFS_path ):

	sfs = []
	ll = os.popen("\\"+simSFS_path+"simSFS s "+str(n)+" "+str(t)+" "+str(rep)).readlines()
	for l in ll:
		sfs.append(float(l[:-1]))

	return sfs


#linearK: simulate (via simSFS program) SFS under the Linear model (Kingman that started at time t (it was 1 individual before that) and linearly grew until now)
#arguments: sample size, parameter t for the demography, number of replicates for the SFS simulation, path to the simSFS program
#return: SFS (list) under the Linear model
def linearK( n, t, rep, simSFS_path ):

	sfs = []
	ll = os.popen("\\"+simSFS_path+"simSFS l "+str(n)+" "+str(t)+" "+str(rep)).readlines()
	for l in ll:
		sfs.append(float(l[:-1]))

	return sfs
	

#expoK: simulate (via simSFS program) SFS under the Exponential model (Kingman with exponential growth at rate 1/t)
#arguments: sample size, parameter t for the demography, number of replicates for the SFS simulation, path to the simSFS program
#return: SFS (list) under the Exponential model
def expoK( n, t, rep, simSFS_path):

	sfs = []
	ll = os.popen("\\"+simSFS_path+"simSFS e "+str(n)+" "+str(t)+" "+str(rep)).readlines()
	for l in ll:
		sfs.append(float(l[:-1]))

	return sfs


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
		
		elif list[0][:3]=='rep':
			rep=int(list[1])
			
		elif list[0][:11]=='range_start':
			range_start=float(list[1])
			
		elif list[0][:9]=='range_end':
			range_end=float(list[1])
			
		elif list[0][:4]=='step':
			step=float(list[1])
		
		elif list[0][:4]=='path':
			path=re.sub(r'"','',(str(list[1][:-1]).strip()))
		
		elif list[0][:8]=='sfs_file':
			sfs_file=re.sub(r'"','',(list[1][:-1].strip()))
			
		elif list[0][:11]=='sample_size':
			n=int(list[1])
			
		elif list[0][:11]=='simSFS_path':
			simSFS_path=re.sub(r'"','',(list[1][:-1].strip()))

run_name=sfs_file.split('/')[-1]

print "-------------SFS fitting-------------"

#
# read sfs_file (the SFS should already be folded)
#

F_sfs_obs = read_sfs_file( sfs_file );
F_sfs_obs_norm = norm_clean(F_sfs_obs,clear)

#
# Kingman aka uniform on t
#
sfs_kgm = Kingman_SFS( n )
F_sfs_kgm = fold( sfs_kgm )
F_sfs_kgm_norm = norm_clean(F_sfs_kgm,clear)


#
# SFS fitting
#

file_fit = open(path+run_name+".fit", "w");
file_fit.write("#t d_BirthDeath d_Linear d_Conditioned d_Sudden d_Exponential\n")

dmin=-1
dmin2=-1
dmin3=-1
dmin4=-1
dmin5=-1

range_start*=1.0/step
range_end*=1.0/step

print "Parameter optimization"

for x in xrange(int(range_start),int(range_end)):

	#
	# Compute t
	#
	
	t = (x+0.0)*step
	print "Parameter value :",t
	

	#
	# Get the birth_death critical SFS
	#
	sfs_bd=BD_SFS( n , t*n )
	F_sfs_bd = fold( sfs_bd )
	F_sfs_bd_norm=norm_clean(F_sfs_bd,clear)
	d = sqr_dev( F_sfs_obs_norm, F_sfs_bd_norm)

	if dmin == -1 or d<dmin:
		dmin = d
		tmin = t

	
	#
	# Get the conditionned Kingman
	#
	sfs_cK= condK( n , t, rep, simSFS_path )
	F_sfs_cK=fold(sfs_cK)
	F_sfs_cK_norm=norm_clean(F_sfs_cK,clear)
	d2 = sqr_dev( F_sfs_obs_norm, F_sfs_cK_norm)

	if dmin2 == -1 or d2<dmin2:
		dmin2 = d2
		tmin2 = t
	
	#
	# Get the sudden Kingman
	#
	sfs_fK=foundedK( n , t, rep, simSFS_path )
	F_sfs_fK =fold(sfs_fK)
	F_sfs_fK_norm=norm_clean(F_sfs_fK,clear)
	d3 = sqr_dev( F_sfs_obs_norm, F_sfs_fK_norm)

	if dmin3 == -1 or d3<dmin3:
		dmin3 = d3
		tmin3 = t


	#
	# Get the linear Kingman
	#
	sfs_lK=linearK( n , t, rep, simSFS_path )
	F_sfs_lK = fold(sfs_lK)
	F_sfs_lK_norm=norm_clean(F_sfs_lK,clear)
	d4 = sqr_dev( F_sfs_obs_norm, F_sfs_lK_norm)

	if dmin4 == -1 or d4<dmin4:
		dmin4 = d4
		tmin4 = t


	#
	# Get the exponential Kingman
	#

	sfs_eK = expoK( n , t, rep, simSFS_path)
	F_sfs_eK=fold(sfs_eK)
	F_sfs_eK_norm=norm_clean(F_sfs_eK,clear)
	d5 = sqr_dev( F_sfs_obs_norm, F_sfs_eK_norm)

	if dmin5 == -1 or d5<dmin5:
		dmin5 = d5
		tmin5 = t
	
	
	file_fit.write( "%f %f %f %f %f %f\n"%(t,d, d4, d2, d3, d5) )
	
file_fit.close()

# Report the best fit

print "-------------Best fit computation-------------"
print "Birth-Death"
sfs_bd=BD_SFS( n , tmin*n )
F_sfs_bd = fold( sfs_bd )
F_sfs_bd_norm=norm_clean(F_sfs_bd,clear)

print "Conditioned"
sfs_cK = condK( n , tmin2, rep, simSFS_path )
F_sfs_cK = fold(sfs_cK)
F_sfs_cK_norm=norm_clean(F_sfs_cK,clear)

print "Sudden"
sfs_fK = foundedK( n , tmin3, rep, simSFS_path )
F_sfs_fK =fold(sfs_fK)
F_sfs_fK_norm=norm_clean(F_sfs_fK,clear)

print "Linear"
sfs_lK = linearK( n , tmin4, rep, simSFS_path )
F_sfs_lK = fold(sfs_lK)
F_sfs_lK_norm=norm_clean(F_sfs_lK,clear)

print "Exponential"
sfs_eK = expoK( n , tmin5, rep, simSFS_path)
F_sfs_eK = fold(sfs_eK)
F_sfs_eK_norm=norm_clean(F_sfs_eK,clear)

with open(path+run_name+".models", "w") as file_model:
	file_model.write("#SFS\tKingman     \tBirthDeath(t)      \tLinear(t)      \tConditioned(t)      \tSudden(t)      \tExponential(t)\n" )
	file_model.write("%s\t%.2e\t%.2e (t=%.2f)\t%.2e (t=%.2f)\t%.2e (t=%.2f)\t%.2e (t=%.2f)\t%.2e (t=%.2f)\n"%(run_name, sqr_dev( F_sfs_obs_norm, F_sfs_kgm_norm),  sqr_dev( F_sfs_obs_norm, F_sfs_bd_norm), tmin, sqr_dev( F_sfs_obs_norm, F_sfs_lK_norm), tmin4, sqr_dev( F_sfs_obs_norm, F_sfs_cK_norm), tmin2, sqr_dev( F_sfs_obs_norm, F_sfs_fK_norm), tmin3, sqr_dev( F_sfs_obs_norm, F_sfs_eK_norm), tmin5 ) )
	file_model.write('.\t.\t.\t'+str(sum(sfs_lK))+'\t'+str(sum(sfs_cK))+'\t'+str(sum(sfs_fK))+'\t'+str(sum(sfs_eK))+'\n')

with open(path+run_name+".bestfit", "w") as file_bestfit:
	file_bestfit.write("i ObservedSFS Kingman BirthDeath Linear Conditioned Sudden Exponential\n" )
	for k in range( len(F_sfs_obs_norm) ):
		file_bestfit.write("%d %f %f %f %f %f %f %f\n"%( k+1+clear, F_sfs_obs_norm[k],  F_sfs_kgm_norm[k], F_sfs_bd_norm[k], F_sfs_lK_norm[k] , F_sfs_cK_norm[k] , F_sfs_fK_norm[k] , F_sfs_eK_norm[k]) )
