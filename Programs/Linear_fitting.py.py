#!/usr/bin/python
# -*- coding: utf-8 -*-


from math import *
import sys
import os
import re
import numpy


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

#linearK: simulate (via simSFS program) SFS under the Linear model (Kingman that started at time t (it was 1 individual before that) and linearly grew until now)
#arguments: sample size, parameter t for the demography, number of replicates for the SFS simulation, path to the simSFS program
#return: SFS (list) under the Linear model
def linearK( n, t, rep, simSFS_path ):

	sfs = []
	ll = os.popen("\\"+simSFS_path+"simSFS l "+str(n)+" "+str(t)+" "+str(rep)).readlines()
	for l in ll:
		sfs.append(float(l[:-1]))

	return sfs

#linearK_topology: simulate (via SiteFrequencySpectrum program) SFS under the Linear model, with topology simulation
#arguments: sample size, parameter t for the demography, number of loci for the SFS simulation
#return: SFS (list) under the Linear model
def linearK_topology( n, t, rep ):

	a = 1.0/t

	sfs = []
	ll = os.popen("\./SimulTrees/SiteFrequencySpectrum -l "+str(a)+" "+str(n)+" 100 "+str(rep)+" | awk '/^[0-9]/{print $2}'").readlines();
	for l in ll:
		sfs.append( float(l[:-1]) )

	return sfs

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------

#
#Read parameter values in commandline
#

loci=int(sys.argv[1])
print "Number of loci :",str(loci)

n=int(sys.argv[2])
print "Sample size :",str(n)

tau=float(sys.argv[3])
print "Foundation time :",str(tau)

#
# Create "observed" linear SFS
#

sfs_obs=linearK_topology( 2*n , tau , loci )
sfs_obs_norm=[x/sum(sfs_obs) for x in sfs_obs]

#
# Kingman aka uniform on t
#
sfs_kgm = Kingman_SFS( 2*n )
sfs_kgm_norm = [x/sum(sfs_kgm) for x in sfs_kgm]


#
# Optimization of t value (Newton Raphson method)
#

REP=int(sys.argv[4])
print "Number of replicates for SFS simulation :",str(REP)
t=1.0
eps=float(sys.argv[5])
print "Epsilon for optimization :",str(eps)

step=1.0
stop=float(sys.argv[6])
print "Stop threashhold :",str(stop),"\n"

print "Observed (simulated) SFS :"
print sfs_obs_norm

while abs(step) > stop:
	print "t :",str(t)
	print "step =",str(step),'\n'
	
	#f(t+2eps)
	sfs_lK = linearK( 2*n, (t+2*eps) , REP )
	sfs_lK_norm=[x/sum(sfs_lK) for x in sfs_lK]
	d_tp2eps=sqr_dev(sfs_obs_norm,sfs_lK_norm)

	#f(t)
	sfs_lK = linearK( 2*n, (t) , REP )
	sfs_lK_norm=[x/sum(sfs_lK) for x in sfs_lK]
	d_t=sqr_dev(sfs_obs_norm,sfs_lK_norm)

	#f(t-2eps)
	sfs_lK = linearK( 2*n, (t-2*eps) , REP )
	sfs_lK_norm=[x/sum(sfs_lK) for x in sfs_lK]
	d_tm2eps=sqr_dev(sfs_obs_norm,sfs_lK_norm)
	
	prim_d_tpeps = ( d_tp2eps - d_t ) / (2*eps)

	prim_d_tmeps = ( d_t - d_tm2eps ) / (2*eps)

	sec_d_t = ( prim_d_tpeps - prim_d_tmeps ) / (2*eps)

	prim_d_t = (prim_d_tpeps + prim_d_tmeps) / 2

	step = -prim_d_t/sec_d_t
	t+=step
	
print "t final :",str(t)
print "step final :",str(step)
