This folder contains the programs to reproduce the results.

`SFS_fitting.py`: fitting the observed Yoruba SFS with five demographic models, by minimizing a least square distance between observed and simulated SFS, depending on the parameter value t

* command line: `python SFS_fitting.py SFS_fitting_parameters.py`

* takes in argument the parameter file containing all necessary parameter values to run the program (see this file for information on the necessary parameters)

* output: 3 files (file names take the observed SFS file name (given in the parameter file) and add a specific extension)
  * `.fit`: contains the least square distance d between the observed SFS and the simulated SFS under each model, for each tested value of the parameter t  
   first column: parameter t value  
   columns 2 to 6 : d value for the models Birth-Death ; Linear ; Conditioned ; Sudden ; Exponential  
   Columns are separated with a single space, the first line gives the columns legend.  

  * `.models`: contains the minimum least square distance d found for each model, with the parameter value t for which this minimal distance was found, and the sum of the optimized SFS for Kingman derived models  
   first line: column legend  
   second line: minimal d value (parameter t value for which this minimal d was found)  
   third line: mean sum of the optimized SFS for the Kingman derived models (all but Birth-Death)  
   Columns are separated with tabulations.

  * `.bestfit`: contains the normalized observed SFS and simulated SFS for the parameter value minimizing the least square distance d (bestfit SFS under each model).  
   first line: column legend  
   following lines: normalized folded SFS value for each bin i (first column)  
   Columns are separated with a single space.

***

`Conditioned_trajectory.py`: simulating a fixation trajectory for the Conditioned model, with the optimized parameter t value

* command line: `python Conditioned_trajectory.py Conditioned_trajectory_parameters.py`

* takes in argument the parameter file containing all necessary parameter values to run the program (see this file for information on the necessary parameters)

* output: `Conditioned_trajectory.txt` containing the mean trajectory. First column gives time and second columns gives the mean frequency of the allele.

***

`BD_trajectory.py`: simulating a fixation trajectory for the critical Birth-Death model, with the optimized parameter t value

* command line: `python BD_trajectory.py BD_trajectory_parameters.py`

* takes in argument the parameter file containing all necessary parameter values to run the program (see this file for information on the necessary parameters)

* output: `BD_trajectory.txt` containing the mean trajectory. First column gives time and second columns gives the mean frequency of the allele.

***

`Linear_fitting.py`: fitting a SFS under the Linear growth model, simulated with the topology simulation method.

* command line: `python Linear_fitting.py l n t r e s`

* arguments of the command line:
  * `l` is the number of loci for the topology simulation of the SFS that will be fitted
  * `n` is the sample size
  * `t` is the foundation time of linear growth
  * `r` is the number of replicates for the SFS simulations during the fitting
  * `e` is the epsilon value for the Newton-Raphson optimization method
  * `s` is the stop threashhold for the optimization

* output: prints in the terminal the parameter values, normalized SFS simulated and step values of the optimization (see `Output_files/example_LinearSFS_fitting.txt` for an example)

