# Tidally locked case decriptions

Tidally lokced exoplanets are exoplanets that synchronously rotate around their host star, meaning their orbital period is the same as their rotational period. This is even expected to happen to potentially habiatble exoplanets around M dwarf stars and some K dwarf stars. These exoplanets are epected to have a zero (or near zero) orbital eccentrity and obliquity (axial tilt).

In each case (apart from the no TL cases, which stands for no tidal locking), the obliquity and eccentricty in the model have been set to 0, and the solar zenith angles have been fixed so that the star is fixed in the sky. 

The standard initial simulations placed the substellar point over 180&deg;  longitude,  0&deg;  latitude, which is in the Pacific ocean. Each simulation uses the BWma1850 WACCM6 configuration. These simulations assume an initial pre-industrial (PI) Earth composition with 1000 hPa surface pressure. SSPL is notation for sub-stellar point over land, where the substellar point was placed at 30&deg;  longitude,  0&deg;  latitude, which is over Africa. 0.1% PAL is 1000 times less than the present atmospheric level of oxygen. 

The following simulations have been used for the paper titled: 'Uncertainties from stellar UV in 3D simulations of TRAPPIST-1 e lead to ambiguities in the interpretation of observations', that will be shortly submitted:

<ins>**Table 1**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-2</sup>] | Spectrum source | Sub stellar point                    | 
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ |
| TP-1 e P19 PI       | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e P19 0.1% PAL | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e P19 no TL    | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e W21 PI       | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) |

## Other simulations 

The following simulations have also been performed, or are currently in progress:

<ins>**Table 2**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-2</sup>] | Spectrum source | Sub stellar point                    |
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ |
| TP-1 e P19 PI SSPL | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 30&deg;  longitude (Africa) |
| TP-1 e W21 0.1% PAL | W21 spectrum | 900   | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e W21 no TL | W21 spectrum  | 900   | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) |
| PC b PI             | Mega-MUSCLES GJ 551 | 900   | [Mega-MUSCLES Treasury Survey GJ 551](https://archive.stsci.edu/prepds/muscles/)  | 30&deg; longitude (Pacific Ocean) |

## Fixed lower boundary conditions

The standard globally averaged lower boundary conditions are given in the following table:

<ins>**Table 3**</ins>

| Chemical | Mixing ratio |
| -------- | ------------ |
| CCL<sub>4</sub> | 2.50e-14 |
| CF<sub>2</sub>CLBR | 4.45e-15 |
| CFC11eq | 3.21e-11 |
| CH<sub>2</sub>BR<sub>2</sub> | 1.20e-12 |
| CH<sub>3</sub>BR | 5.30e-12 |
| CH<sub>3</sub>CL | 4.57e-10 |
| CH<sub>4</sub> | 8.08e-07 |
| CHBR<sub>3</sub> | 1.20e-12 |
| CO<sub>2</sub> | 0.000284 |
| H<sub>2</sub> | 5.00e-07 |
| OCS | 3.37e-10 |
| N<sub>2</sub>O | 2.73e-07 |
| Total BROY | 1.13e-11 |
| Total CLOY | 4.57e-10 |
| Total FOY | 8.89e-15 |
| O<sub>2</sub> | 0.21 |

For the 0.1% PAL simulations, the O<sub>2</sub> boundary condition was 0.00021.

## Changing the stellar spectrum

We used the [SWIM](https://github.com/jack-eddy-symposium/exoplanetary-impact/tree/main/SWIM) code, or a similar method, to change the stellar spectrum and rescale it to specific exoplanets.

## Initial conditions and getting through model crashes

Taking WACCM6 from an Earth System Model to a model that simulates an M dwarf terrestrial exoplanet is a huge change in the model, and it can be unstable. Possible ways of getting through crashes are changing CLUBB parameters, decreasing the timestep, and reducing the rotation in steady jumps through several branch cases.

## Location of restart files

Restart files are currently located on ARC4, the supercomputer at the University of Leeds.
