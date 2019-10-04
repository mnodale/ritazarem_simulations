# ritazarem_simulations

This project contains:
1) abstract for ICTMC 2019
2) poster for ICTMC 2019
3) source code (source_simulations.R & source_results.R) used to perform the simulations and produce output figures
4) input parameters (inputs.RData) used in the simulations & derived from the trial data
5) output data for each scenario described in the poster (sims_1.RData, etc.)
6) full set of figures for each model (Mean + SE, MSE, Bias, Power, Coverage & Convergence)
7) additional materials

To replicate the simulations described in the poster run "source_simulations.R", ensuring the path to the working directory is up-to-date & inputs.RData is located in the working directory. Output data for each scenarios will be saved in the working directory.
Be patients it may take some time (hours, if not days)!

Running "source_results.R" afterwords (without resetting the global enviroment) will produce a full set of figures for each scenario. Figures are placed in a (pre-existing) folder ("/Figures") inside the working directory & are clearly labelled for ease of identification.

Marianna Nodale, October 2019