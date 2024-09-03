# Tidally locked case descriptions

Tidally locked exoplanets are exoplanets that synchronously rotate around their host star, meaning their orbital period is the same as their rotational period. This is expected to happen to potentially habiatble exoplanets around M dwarf stars and some K dwarf stars. These exoplanets are expected to have a zero (or near zero) orbital eccentrity and obliquity (axial tilt).

In each case (apart from the no TL cases, which stands for no tidal locking), the solar zenith angles have been fixed so that the star is unmoving in the sky. In each case the obliquity and eccentricty in the model have been set to 0. The spectrum in each case also has to be changed, depending on which exoplanet one wishes to simulate.

The standard initial simulations placed the substellar point over 180&deg;  longitude,  0&deg;  latitude, which is in the Pacific ocean. Each simulation uses the BWma1850 WACCM6 configuration. These simulations assume an initial pre-industrial (PI) Earth composition with 1000 hPa surface pressure. SPL is notation for substellar point over land, where the substellar point was placed at 30&deg;  longitude,  0&deg;  latitude, which is over Africa. 0.1% PAL is 1000 times less than the present atmospheric level of oxygen, where PAL stands for present atmospheric level.

## Publications

* Cooke _et al_. (2023), Degenerate interpretations of O<sub>3</sub> spectral features in exoplanet atmosphere observations due to stellar UV uncertainties: a 3D case study with TRAPPIST-1e, _ApJ_, [https://iopscience.iop.org/article/10.3847/1538-4357/ad0381/meta](https://iopscience.iop.org/article/10.3847/1538-4357/ad0381/meta).
* Cooke _et al_. (2024), Lethal Surface Ozone Concentrations Are Possible on Habitable Zone Exoplanets, _Planet. Sci. J._ [https://iopscience.iop.org/article/10.3847/PSJ/ad53c3](https://iopscience.iop.org/article/10.3847/PSJ/ad53c3).
 
## TRAPPIST-1 e simulations

The TRAPPIST-1 e simulations assume the following parameters which can be edited in user_nl_cam:

* **omega** (rotation rate) = 0.00001192352565 seconds<sup>-1</sup>
* **sday** (number of seconds in a sidereal day) = 526957 seconds
* **gravit** (gravitational acceleration) = 9.1454 m s<sup>-2</sup>
* **rearth** (radius of the planet) = 5797610 metres

The following simulations have been used for the paper titled: 'Uncertainties from stellar UV in 3D simulations of TRAPPIST-1e lead to ambiguities in the interpretation of observations', that has been submitted:

<ins>**Table 1 (TRAPPIST-1 e simulations)**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-2</sup>] | Spectrum source | Substellar point                    | O<sub>2</sub> mixing ratio [PAL] |
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ | ---------|
| TP-1 e P19 PI       | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) | 1.000 |
| TP-1 e P19 PI SPL       | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 30&deg;  longitude (Africa) |  1.000 |
| TP-1 e P19 PI noTL       | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | N/A | 1.000 |
| TP-1 e P19 10% PAL     | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) | 0.100 |
| TP-1 e P19 1% PAL    | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) | 0.010 |
| TP-1 e P19 0.1% PAL       | P19 spectrum | 900  | [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) | 180&deg;  longitude (Pacific ocean) | 0.001 |
| TP-1 e W21 PI       | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) | 1.000 |
| TP-1 e W21 PI SPL       | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 30&deg;  longitude (Pacific ocean) | 1.000 |
| TP-1 e W21 PI noTL  | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | N/A | 1.000 |
| TP-1 e W21 10% PAL  | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) | 0.100 |
| TP-1 e W21 1% PAL | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) | 0.010 |
| TP-1 e W21 0.1% PAL  | W21 spectrum | 900  | [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) | 180&deg;  longitude (Pacific ocean) | 0.001 |

## Proxima Centauri b simulations

The Proxima Centauri b simulations assume the following parameters:

* **omega** (rotation rate) = 0.000006502178685 seconds<sup>-1</sup>
* **sday** (number of seconds in a sidereal day) = 966320 seconds
* **gravit** (gravitational acceleration) = 12.2 m s<sup>-2</sup>
* **rearth** (radius of the planet) = 6816970 metres

<ins>**Table 2 (Proxima Centauri b simulations)**</ins>

| Simulation name | Stellar spectrum used   | Total irradiance [Wm<sup>-2</sup>] | Spectrum source | Substellar point  | O<sub>2</sub> mixing ratio [PAL] |
| --------------- | ----------------------- | ---------------------------------- | --------------- | ------------------------------------ | ---------|
| PC b PI       | GJ 551 MUSCLES spectrum | 884  | [GJ 551 MUSCLES spectrum](https://archive.stsci.edu/prepds/muscles/) | 180&deg;  longitude (Pacific ocean) | 1.000 |
| PC b PI SPL      | GJ 551 MUSCLES spectrum | 884  | [GJ 551 MUSCLES spectrum](https://archive.stsci.edu/prepds/muscles/) | 30&deg;  longitude (Africa) | 1.000 |
| PC b 10% PAL      | GJ 551 MUSCLES spectrum | 884  | [GJ 551 MUSCLES spectrum](https://archive.stsci.edu/prepds/muscles/) | 180&deg;  longitude (Pacific ocean) | 0.100 |
| PC b 1% PAL      | GJ 551 MUSCLES spectrum | 884  | [GJ 551 MUSCLES spectrum](https://archive.stsci.edu/prepds/muscles/) | 180&deg;  longitude (Pacific ocean) | 0.010 |
| PC b 0.1% PAL     | GJ 551 MUSCLES spectrum | 884  | [GJ 551 MUSCLES spectrum](https://archive.stsci.edu/prepds/muscles/) | 180&deg;  longitude (Pacific ocean) | 0.001 |

## Source file modifications

In src.cam, the following files require modification:
chemistry.F90, mo_jshort.F90, mo_photo.F90, mo_tgcm_ubc.F90, orbit.F90, radiation.F90,  upper_bc.F90

In src.cice, the following files require modification:
ice_therm_shared.F90

In src.share, the following files require modification:
shr_const_mod.F90

radiation.F90 and orbit.F90 fix the solar zenith angle. upper_bc.F90 and mo_tgcm_ubc.F90 are required to change the upper boundary conditions for several chemicals. chemistry.F90 sets the minimum values for ozone and methane. mo_photo.F90 and mo_jshort.F90 are required to include absorption from carbon dioxide and water vapour in the Schumann-Runge bands. ice_therm_shared.F90 sets a new minimum temperature for ice - in some simulations this minimum temperature may need changing if the ice gets colder than -100 degrees Celcius. shr_const_mod.F90 sets the rotation rate, gravity etc. in the ocean and other models.

## Fixed lower boundary conditions

The standard globally averaged lower boundary conditions are given in the following table:

<ins>**Table 3 (lower boundary conditions)**</ins>

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

For the lower oxygen simulations, the O<sub>2</sub> boundary condition was changed. For example, in the 0.1% PAL simulations, the O<sub>2</sub> boundary condition was 0.00021.

## Changing the stellar spectrum

We used the [SWIM](https://github.com/jack-eddy-symposium/exoplanetary-impact/tree/main/SWIM) code, or a similar method, to change the stellar spectrum and rescale it to specific exoplanets.

## Initial conditions and getting through model crashes

Taking WACCM6 from an Earth System Model to a model that simulates an M dwarf terrestrial exoplanet is a huge change in the model, and it can be unstable. Possible ways of getting through crashes are changing CLUBB parameters, decreasing the timestep, and reducing the rotation in steady jumps through several branch cases. Other problems encountered included the minimum snow temperature being set in the model at -100 degrees Celcius. In some siulations, it is necessary to decrease this to a lower value in the ice_therm_shared.F90 file.

## Location of restart files

Restart files are currently located on ARC4, the supercomputer at the University of Leeds. Please contact Greg Cooke or Dan Marsh for more details.
