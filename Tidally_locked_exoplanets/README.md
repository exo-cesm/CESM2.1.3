# Tidally locked case decriptions
Tidally lokced exoplanets are exoplanets that synchronously rotate around their host star, meaning their orbital period is the same as their rotational period. This is even expected to happen to potentially habiatble exoplanets around M dwarf stars and some K dwarf stars. These exoplanets are epected to have a zero (or near zero) orbital eccentrity and obliquity (axial tilt).

In each case, the obliquity and eccentricty in the model have been set to 0, and the solar zenith angles have been fixed so that the star is fixed in the sky. 

When the simulations are finished, the details will be in the table below:

<ins>**Table 1**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-3</sup>] | Spectrum source | Sub stellar point                    | 
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ |
| TP-1 e P19 PI       | P19 spectrum Peacock et al. (2019) | 900  |  | 180&deg;  longitude (Pacific ocean) |
| TP-1 e P19 0.1% PAL | P19 spectrum Peacock et al. (2019) | 900 |                 | 30&deg; longitude (Africa) |
| TP-1 e P19 no TL    | P19 spectrum Peacock et al. (2019) | 900  |                 | 180&deg;  longitude (Pacific ocean) |
| TP-1 e W21 PI       | W21 spectrum Wilson et al. (2021)  | 900   |                 | 30&deg; longitude (Africa) |
| PC b PI             | Mega-MUSCLES  | 900   |                 | 30&deg; longitude (Pacific Ocean) |

## Changing the stellar spectrum

We used the [SWIM](https://github.com/jack-eddy-symposium/exoplanetary-impact/tree/main/SWIM) code to change the stellar spectrum and rescale it to specific exoplanets.

## Initial conditions and getting through model crashes

Taking WACCM6 from an Earth System Model to a model that simulates an M dwarf terrestrial exoplanet is a huge change in the model, and it can be unstable. Possible ways of getting through crashes are changing CLUBB parameters, decreasing the timestep, and reducing the rotation in steady jumps through several branch cases.

## Location of restart files

Restart files are currently located on ARC4, the supercomputer at the University of Leeds.
