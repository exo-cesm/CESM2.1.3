# Code modifications in CESM2.1.3 for exoplanet and paleoclimate simulations
## Authors: Greg Cooke (University of Leeds) & Dan Marsh (University of Leeds, NCAR)

## Introduction

This repository contains the code modifications necessary for simulating plaeoclimates and exoplanets with CESM2.1.3, specifically with WACCM6.

It details how to set up simulations, why the changes were made, as well as problems that were encountered.

Namelist file changes, boundary conditions, solar spectra input files, and source code modifications will be given.

The location of restarts directories will also be given and are currently on the ARC4 supercomuter at the University of Leeds.

## Varied O<sub>2</sub> simulations (paleoclimates / Earth-analogue exoplanets)

These simulations were used to calculate new ozone estimates from 0.1% the present atmospheric level (PAL) of O<sub>2</sub> to 150% PAL of O<sub>2</sub> - see [Cooke et al. 2022a](https://doi.org/10.1098/rsos.211165). They were then used to predict time-varying direct imaging of Earth-analogue exoplanets - see [Cooke et al. 2022b](https://doi.org/10.1093/mnras/stac2604).

Instructions on how to set up these cases are found in the _O2\_Earth\_analogues_ folder

## Tidally locked exoplanet scenarios

Instructions on how to set up these cases are found in Tidally\_locked\_exoplanets

