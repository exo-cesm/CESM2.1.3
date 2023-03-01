# Tidally locked case decriptions
Tidally lokced exoplanets are exoplanets that synchronously rotate around their host star, meaning their orbital period is the same as their rotational period. This is even expected to happen to potentially habiatble exoplanets around M dwarf stars and some K dwarf stars. These exoplanets are epected to have a zero (or near zero) orbital eccentrity and obliquity (axial tilt).

In each case, the obliquity and eccentricty in the model have been set to 0, and the solar zenith angles have been fixed so that the star is fixed in the sky. 

The following simulations have been used for the paper titled: 'Uncertainties from stellar UV in 3D simulations of TRAPPIST-1 e lead to ambiguities in the interpretation of observations', that will be shortly submitted:

<ins>**Table 1**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-3</sup>] | Spectrum source | Sub stellar point                    | 
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ |
| TP-1 e P19 PI       | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e P19 0.1% PAL | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e P19 no TL    | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e W21 PI       | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) |

## Other simulations 

The following simulations have also been performed, or are currently in progress:

<ins>**Table 2**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-3</sup>] | Spectrum source | Sub stellar point                    |
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ |
| TP-1 e P19 PI SSPL | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e W21 0.1% PAL | W21 spectrum | 900   | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) |
| TP-1 e W21 no TL | W21 spectrum Wilson et al. (2021)  | 900   | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) |
| PC b PI             | Mega-MUSCLES  | 900   | [Mega-MUSCLES Treasury Survey GJ 551](https://archive.stsci.edu/prepds/muscles/)  | 30&deg; longitude (Pacific Ocean) |


## Changing the stellar spectrum

We used the [SWIM](https://github.com/jack-eddy-symposium/exoplanetary-impact/tree/main/SWIM) code to change the stellar spectrum and rescale it to specific exoplanets.

## Initial conditions and getting through model crashes

Taking WACCM6 from an Earth System Model to a model that simulates an M dwarf terrestrial exoplanet is a huge change in the model, and it can be unstable. Possible ways of getting through crashes are changing CLUBB parameters, decreasing the timestep, and reducing the rotation in steady jumps through several branch cases.

## Location of restart files

Restart files are currently located on ARC4, the supercomputer at the University of Leeds.
