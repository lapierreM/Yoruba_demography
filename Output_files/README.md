This folder contains examples of output files for the programs.

- 3 outputs of the `SFS_fitting.py` program :
  * `example_FoldSFS_YRI.txt.fit`: values of d (least square distance) between observed Yoruba SFS and the five models depending on parameter t value

  * `example_FoldSFS_YRI.txt.models`: minimum values of d (least square distance) between observed Yoruba SFS and the five models for the optimized parameter t value (given between parenthesis). The third line gives the total inferred tree length for the Kingman derived models. 

  * `example_FoldSFS_YRI.txt.bestfit`: bestfit SFS under the five models, observed SFS and standard Kingman SFS. The first column gives the SFS bin index. SFS are normalized.

- `example_Conditioned_trajectory.txt`: trajectory of fixation under the Conditioned model, output from the `Conditioned_trajectory.py` program. First column gives time and second columns gives the mean frequency of the allele.

- `example_BD_trajectory.txt`: trajectory of fixation under the critical Birth-Death model, output from the `BD_trajectory.py` program. First column gives time and second columns gives the mean frequency of the allele.

- `example_LinearSFS_fitting.txt`: output of the `Linear_fitting.py` program. Contains the parameter values (number of loci, sample size, foundation time, number of replicates for SFS simulation, epsilon for derivative calculation and stop threashhold), followed by the normalized simulated SFS, and the step values of the optimization (for each step: current value of the parameter t, and value of the step made just before).

- `example_swp_SFS_YRI`: output of the `simSFS` program run with the option -S on the stairway plot output summary file for the Yoruba population (`Data/YRI_start2_200files_summary`). Expected SFS under the demography reconstructed by the Stairway Plot method for the Yoruba population.
