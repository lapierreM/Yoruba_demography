The simSFS program simulates SFS under the Kingman coalescent, with or without demography (see article for details on method and models).

Command line: `./simSFS [options] sample_size number_rep`

where :
- `[options]` are:
  - `-x number`: set seed to `number` (for random generator)
  - `-c float`: Conditioned model (Kingman with Tmrca smaller than float)
  - `-s float`: Sudden model (sudden growth at time float)
  - `-e float`: Exponential model (exponential growth at rate float)
  - `-l float`: Linear model (linear growth with foundation time float)
  - `-S file`: read input demography grom stairway plot output summary file

- sample_size is the sample size of the simulated SFS

- number_rep is the number of repetitions for the mean SFS computation

Output: list (one value per line), simulated SFS

Other output are possible (must be un-commented in the code)
