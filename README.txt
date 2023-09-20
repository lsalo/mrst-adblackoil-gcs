This module provides examples and advanced functionality for three-dimensional 
simulation of geologic carbon sequestration with miscibility, relative 
permeability hysteresis and molecular diffusion. All examples and testing 
have been performed using MRST's ad-blackoil module (where we assign
brine-like properties to the oil phase, so that the gas phase (CO2-rich)
can dissolve in it).
The module is organized as follows:
1. Calculation of PVT properties: thermodynamic formulation following 
   Hassanzadeh et al. (2008). Includes functions, validation and example.
2. Relative permeability hysteresis: calculation of trapped gas saturation 
   and bounding imbibition curve using Land’s model, implementation of 
   Killough's 1974 model as in ECLIPSE for scanning curves. Includes 
   the hysteresis class and an example which uses the PUNQ model and 
   compares the result with that of ECLIPSE (Juanes et al., WRR 2006)
3. Molecular diffusion (includes the diffusion class and an example).

