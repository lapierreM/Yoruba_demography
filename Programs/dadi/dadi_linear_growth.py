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
import numpy
import dadi

# Linear growth demographic function

def lin_growth(params, (n1), pts):
	T = params
	# Define the grid we'll use
	xx = dadi.Numerics.default_grid(pts)
	
	# phi for the equilibrium ancestral population
	phi = dadi.PhiManip.phi_1D(xx)

	# now do the population growth event
	nu=0.000000001 # initial population size
	nu_func = lambda t: nu+((1.0-nu)*t)/T # demographic function (linear growth)
	phi = dadi.Integration.one_pop(phi, xx, T, nu=nu_func)
	
	# finally, calculate the spectrum
	sfs = dadi.Spectrum.from_phi(phi, ns, (xx,))
	
	return sfs

#------------------------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------------------------
	
nu=0.000000001 # initial population size

# Load the data
data=dadi.Spectrum.from_file("YRI.fs")
ns=data.sample_sizes

# ignore singletons
ignore='yes'
if ignore=='yes':
	data.mask[1]=True

# These are the grid point settings will use for extrapolation
pts_l=[300,400,500]

# our custom model function
func = lin_growth

# Now let's optimize parameters for this model

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.
# Parameters are: (nu, T)
upper_bound = [10]
lower_bound = [0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [3]

# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)

# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.

print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=1, maxiter=3)
                                   
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finished optimization **************************************************')

print('Best-fit parameters: {0}'.format(popt))

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))


# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure(1)
dadi.Plotting.plot_1d_comp_multinom(model, data)
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
pylab.show()
