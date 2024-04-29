# mrst-adblackoil-gcs
Provides examples and advanced functionality for 3D simulation of geologic CO2 storage as described in Saló-Salgado (2024). 

**Important note**: As of 4/25/2024, this code has been integrated into MRST (release 2024a). Therefore, the code in this repository is no longer needed if you use version 2024a or newer. Code and examples can be found in the co2lab-mit module, please see [the release notes](https://www.sintef.no/projectweb/mrst/download/release-notes-for-mrst-2024a/) for details.

## Citation
### Paper
```
The paper/preprint will appear here when published.
```

### Thesis
```
@phdthesis{salo2024,
  title={Numerical Modeling of Geologic Carbon Dioxide Storage in Faulted Siliciclastic Settings},
  author={Saló-Salgado, Lluís},
  year={2024},
  month={February},
  school={Massachusetts Institute of Technology},
  url={https://hdl.handle.net/1721.1/153718},
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