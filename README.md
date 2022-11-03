# Code modifications in CESM2.1.3 for exoplanet and paleoclimate simulations
## Authors: Greg Cooke (University of Leeds) & Dan Marsh (University of Leeds, NCAR)

## Introduction

This repository contains the code modifications necessary for simulating plaeoclimates and exoplanets with CESM2.1.3, specifically with WACCM6.

It details how to set up simulations, why the changes were made, as well as problems that were encountered.

Namelist file changes, boundary conditions, solar spectra input files, and source code modifications will be given.

The location of restarts directories will also be given and are currently on the ARC4 supercomuter at the University of Leeds.

## Varied O<sub>2</sub> simulations (paleoclimates / Earth-analogue exoplanets)

These simulations were used to calculate new ozone estimates from 0.1% the present atmospheric level (PAL) of O<sub>2</sub> to 150% PAL of O<sub>2</sub> - see [Cooke et al. 2022a](https://doi.org/10.1098/rsos.211165). They were then used to predict time-varying direct imaging of Earth-analogue exoplanets - see [Cooke et al. 2022b](https://doi.org/10.1093/mnras/stac2604).

Since then, we (Dan, Greg, and Doug Kinnison) have updated the source code to include absorption in the Schumann-Runge bands, and also added in previously missing cross sections. Whilst these updates do not make an appreciable difference for the pre-industrial atmosphere, it is now recommended to include these updates when investigating low-O<sub>2</sub> atmospheres.

Instructions on how to set up these cases are found in the [a workaround link](tree/main/O2_Earth_analogues) _O2\_Earth\_analogues_ folder.

## Tidally locked exoplanet scenarios

Instructions on how to set up these cases are found in the _Tidally\_locked\_exoplanets_ folder.

The tidally locked simulations cover several cases for known M dwarf terrestrial exoplanets. The basic modifications include changing the rotation rate of the model, changing the radius and gravitational acceleration, and fixing the solar zenith angle (in order to fix the substellar point). The substellar point is placed either in the middle of the Pacific ocean at 180&deg; longitude, or at 30&deg; longitude in Africa. Additionally, depending on the exoplanet in question, the solar file will need to be changed and scaled to the irradiance that the planet recieves. The Stellar Wind and Irradiance Model ([SWIM](https://github.com/jack-eddy-symposium/exoplanetary-impact/tree/main/SWIM)) has been developed for this purpose.
