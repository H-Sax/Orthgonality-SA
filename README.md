# File Structure
We Order files by Model, Data type (discrete or continous) and then Sampling method. 
## Data Files
Each data file is stored as a .jld file. Each file has the format Sampling_EstimatorDatatype.jld. 

e.g. LatinHypercube_SD.jld is the Latin Hypercube Sampling method with the Sobol estimtor computed on discrete data. 

## Estimators Used 
S - Sobol estimator 
H - Homma estimator 
J - Jansen Estimator 
Ja - Janon Etimator 
