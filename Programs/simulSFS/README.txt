The simSFS program simulates SFS under the Kingman coalescent, with or without demography (see article for details on method and models).

Command line: ./simSFS demography sample_size time number_rep

where :
- demography can be:
	k (Kingman without demography)
	c (Conditioned model)
	s (Sudden model)
	e (Exponential model)
	l (Linear model)

- sample_size is the sample size of the simulated SFS

- time is the characteristic time for demography modeling

- number_rep is the number of repetitions for the mean SFS computation

Output: list (one value per line), simulated SFS

other output are possible (must be un-commented in the code)