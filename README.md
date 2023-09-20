# mrst-adblackoil-gcs
Provides examples and advanced functionality for 3D simulation of geologic CO2 storage as described in Saló-Salgado (2023). This code will be integrated into MRST in the near future, and updates will be provided below.

## Updates
Updates will be posted here after integration into MRST.

## Citation
### Paper
```
The paper/preprint will appear here when published.
```

### Thesis
```
@phdthesis{salo2023,
  title={Numerical Modeling of Geologic Carbon Dioxide Storage in Faulted Siliciclastic Settings},
  author={Saló-Salgado, Lluís},
  year={2023},
  school={Massachusetts Institute of Technology},
  doi={doi will be added here when available},
  note={Chapter 2}
}
```

## Description
This module provides examples and advanced functionality for three-dimensional 
simulation of geologic carbon sequestration with miscibility, relative 
permeability hysteresis and molecular diffusion. All examples and testing 
have been performed using MRST's ad-blackoil module (where we assign
brine-like properties to the oil phase, so that the gas phase (CO2-rich)
can dissolve in it).
The module is organized as follows:
1. Calculation of PVT properties: thermodynamic formulation following 
   Hassanzadeh et al., IJGGC (2008). Includes functions, validation and example.
2. Relative permeability hysteresis: calculation of trapped gas saturation 
   and bounding imbibition curve using Land’s SPE J (1968) model, implementation of 
   Killough's SPE J (1974) model as in ECLIPSE for scanning curves. Includes 
   the hysteresis class and an example which uses the PUNQ model and 
   compares the result with that of ECLIPSE (Juanes et al., WRR 2006)
3. Molecular diffusion (includes the diffusion class and an example inspired by work on the [FluidFlower](https://fluidflower.w.uib.no/)).