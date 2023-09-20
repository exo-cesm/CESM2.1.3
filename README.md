# Code modifications in CESM2.1.3 for exoplanet and paleoclimate simulations
## Authors: Greg Cooke (gjc53@cam.ac.uk, University of Cambridge, University of Leeds) & Dan Marsh (University of Leeds)

## Introduction

This repository contains the code modifications necessary for simulating plaeoclimates and exoplanets with CESM2.1.3, specifically with WACCM6.

It details how to set up simulations, why the changes were made, as well as problems that were encountered.

Namelist file changes, boundary conditions, solar spectra input files, and source code modifications are given.

The location of restarts directories will also be given and are currently on the ARC4 supercomputer at the University of Leeds.

If anything is unclear, please get in touch with Greg at gjc53@cam.ac.uk.

## Publications

* Cooke G. J., Marsh D. R., Walsh C., Youngblood A., _Astrophsical Journal_, Accepted, 2023.
* Ji A., Kasting J. F., Cooke G. J., Marsh D. R. and Tsigaridis K. 2023, Comparison between ozone column depths and methane lifetimes computed by one- and three-dimensional models at different atmospheric O2 levels, _R. Soc. open sci_, 10230056230056, [https://doi.org/10.1098/rsos.230056](https://doi.org/10.1098/rsos.230056).
* G J Cooke, D R Marsh, C Walsh, S Rugheimer, G L Villanueva, Variability due to climate and chemistry in observations of oxygenated Earth-analogue exoplanets, Monthly Notices of the Royal Astronomical Society, Volume 518, Issue 1, January 2023, Pages 206–219, [https://doi.org/10.1093/mnras/stac2604](https://doi.org/10.1093/mnras/stac2604).
* Cooke G. J., Marsh D. R., Walsh C., Black B. and Lamarque J.-F. 2022, A revised lower estimate of ozone columns during Earth’s oxygenated history, _R. Soc. open sci_, 9211165211165, [https://doi.org/10.1098/rsos.211165](https://doi.org/10.1098/rsos.211165).

## Varied O<sub>2</sub> simulations (paleoclimates / Earth-analogue exoplanets)

These simulations were used to calculate new ozone estimates from 0.1% the present atmospheric level (PAL) of O<sub>2</sub> to 150% PAL of O<sub>2</sub> - see [Cooke et al. 2022](https://doi.org/10.1098/rsos.211165). They were then used to predict time-varying direct imaging of Earth-analogue exoplanets - see [Cooke et al. 2023](https://doi.org/10.1093/mnras/stac2604).

Instructions on how to set up these cases are found in the [_O2\_Earth\_analogues\_folder](/O2_Earth_analogues).

A list of cases (where PAL means present atmospheric level) are described in the table below:

<ins>**Table 1**</ins>

| Simulation name   | Description |
| ---------------   | ----------- |
| [Pre-industrial](/O2_Earth_analogues/cases/PI_baseline)    | Standard pre-industrial BWma1850 case |
| [150% PAL](/O2_Earth_analogues/cases/150pc_PAL_O2) | 150% PAL of O<sub>2</sub> |
| [50% PAL](/O2_Earth_analogues/cases/50pc_PAL_O2)           | 50% PAL of O<sub>2</sub>  |
| [10% PAL](/O2_Earth_analogues/cases/10pc_PAL_O2)           | 10% PAL of O<sub>2</sub>  |
| [10% PAL CH<sub>4</sub> em10](/O2_Earth_analogues/cases/10pc_PAL_O2_CH4_em10)    | 10% PAL of O<sub>2</sub> 5000 Tg/yr CH<sub>4</sub> flux |
| [10% PAL CH<sub>4</sub> em1](/O2_Earth_analogues/cases/10pc_PAL_O2_CH4_em1)      | 10% PAL of O<sub>2</sub> 500 Tg/yr CH<sub>4</sub> flux  |
| [10% PAL CH<sub>4</sub> em0.1](/O2_Earth_analogues/cases/10pc_PAL_O2_CH4_em0.1)  | 10% PAL of O<sub>2</sub> 50 Tg/yr CH<sub>4</sub> flux   |
| [5% PAL](/O2_Earth_analogues/cases/5pc_PAL_O2)            | 5% PAL O<sub>2</sub>      |
| [1% PAL](/O2_Earth_analogues/cases/1pc_PAL_O2)            | 1% PAL of O<sub>2</sub>   |
| [1% PAL CH<sub>4</sub> em10](/O2_Earth_analogues/cases/1pc_PAL_O2_CH4_em10)    | 1% PAL of O<sub>2</sub> 5000 Tg/yr CH<sub>4</sub> flux |
| [1% PAL CH<sub>4</sub> em1](/O2_Earth_analogues/cases/1pc_PAL_O2_CH4_em1)      | 1% PAL of O<sub>2</sub> 500 Tg/yr CH<sub>4</sub> flux  |
| [1% PAL CH<sub>4</sub> em0.1](/O2_Earth_analogues/cases/1pc_PAL_O2_CH4_em0.1)  | 1% PAL of O<sub>2</sub> 50 Tg/yr CH<sub>4</sub> flux   |
| [1% PAL 2 Ga YS](/O2_Earth_analogues/cases/)  | 1% PAL of O<sub>2</sub> 2 Gyr younger Sun |
| [1% PAL 2 Ga YS 4x CO<sub>2</sub>](/O2_Earth_analogues/cases/)  | 1% PAL of O<sub>2</sub> 2 Gyr younger Sun with 1120 ppmv CO<sub>2</sub> |
| [1% PAL 2 Ga YS 21 hr](/O2_Earth_analogues/cases/)  | 1% PAL of O<sub>2</sub> 2 Gyr younger Sun with 21 hour day |
| [0.5% PAL](/O2_Earth_analogues/cases/0.5pc_PAL_O2)            | 0.5% PAL of O<sub>2</sub>   |
| [0.1% PAL](/O2_Earth_analogues/cases/0.1pc_PAL_O2)            | 0.1% PAL of O<sub>2</sub>   |

Since these simulations were performed, we (Dan, Greg, and Doug Kinnison) have updated the source code to include absorption in the Schumann-Runge bands, and also added in previously missing cross sections. Whilst these updates do not make an appreciable difference for the pre-industrial atmosphere, it is now recommended to include these updates when investigating low-O<sub>2</sub> atmospheres (for Earth, this occurs at 1\% PAL and lower). For the understanding of why this is important, see [Ji et al., 2023](https://doi.org/10.1098/rsos.230056). The following simulations were performed for [Ji et al., 2023](https://doi.org/10.1098/rsos.230056).

<ins>**Table 2**</ins>

| Simulation name   | Description |
| ---------------   | ----------- |
| 10% PAL 10x lower halogens    | 10% PAL of O<sub>2</sub> with 10x lower amounts of halogen molecules, including CH<sub>3</sub>Cl and CH<sub>3</sub>Br |
| 1% PAL 10x lower halogens     | 1% PAL of O<sub>2</sub> with 10x lower amounts of halogen molecules, including CH<sub>3</sub>Cl and CH<sub>3</sub>Br |
| 1% PAL 1e9x lower halogens    | 1% PAL of O<sub>2</sub> with 1e9x lower amounts of halogen molecules, including CH<sub>3</sub>Cl and CH<sub>3</sub>Br |
| 10% PAL 1D LBC and SRB    | 10% PAL of O<sub>2</sub> with Kasting 1D model lower boundary conditions and absorption for HO<sub>2</sub>O and CO<sub>2</sub> in the Schumann-Runge Bands included |
| 1% PAL 1D LBC and SRB   | 1% PAL of O<sub>2</sub> with Kasting 1D model lower boundary conditions and absorption for HO<sub>2</sub>O and CO<sub>2</sub> in the Schumann-Runge Bands included |
| 0.1% PAL 1D LBC and SRB    | 0.1% PAL of O<sub>2</sub> with Kasting 1D model lower boundary conditions and absorption for HO<sub>2</sub>O and CO<sub>2</sub> in the Schumann-Runge Bands included |

## Tidally locked exoplanet scenarios

Instructions on how to set up these cases are found in the [_Tidally\_locked\_exoplanets\_folder](/Tidally_locked_exoplanets).

<ins>**Table 3**</ins>

| Simulation name       | Description |
| --------------------- | ----------- |
| [TP-1e P19 PI](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_P19_PI) | P19 spectrum, PI composition, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e P19 PI SPL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPL_P19_PI)  | P19 spectrum, PI composition, tidally locked TP-1e, substellar point over Africa |
| [TP-1e P19 PI noTL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_P19_PI_noTL)  | P19 spectrum, PI composition, tidally locked TP-1e, not tidally locked (1 day rotation) 
| [TP-1e P19 10% PAL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_P19_10pc) | P19 spectrum, 10% PAL of O<sub>2</sub>, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e P19 1% PAL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_P19_1pc) | P19 spectrum, 1% PAL of O<sub>2</sub>, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e P19 0.1% PAL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_P19_01pc) | P19 spectrum, 0.1% PAL of O<sub>2</sub>, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e W21 PI](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_W21_PI) | W21 spectrum, PI composition, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e W21 PI SPL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPL_W21_PI) | W21 spectrum, PI composition, tidally locked TP-1e, substellar point over Africa |
| [TP-1e W21 PI noTL](/Tidally_locked_exoplanets/cases/TP_1_e_W21_PI_noTL) | W21 spectrum, PI composition, tidally locked TP-1e, not tidally locked (1 day rotation) |
| [TP-1e W21 10% PAL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_W21_10pc) | W21 spectrum, 10% PAL of O<sub>2</sub>, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e W21 1% PA](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_W21_1pc) | W21 spectrum, 1% PAL of O<sub>2</sub>, tidally locked TP-1e, substellar point over Pacific Ocean |
| [TP-1e W21 0.1% PAL](/Tidally_locked_exoplanets/cases/TP_1_e_SSPO_W21_01pc) | W21 spectrum, 0.1% PAL of O<sub>2</sub>, tidally locked TP-1e, substellar point over Pacific Ocean |
| [PCb PI](/Tidally_locked_exoplanets/cases/PC_b_SSPO_PI) | GJ551 MUSCLES spectrum, PI composition, tidally locked PC b, substellar point over Pacific Ocean |
| [PCb PI SPL](/Tidally_locked_exoplanets/cases/PC_b_SSPL_PI) | GJ551 MUSCLES spectrum, PI composition, tidally locked PC b, substellar point  over Africa |
| [PCb 10% PAL](/Tidally_locked_exoplanets/cases/PC_b_SSPO_10pc) | GJ551 MUSCLES spectrum, 10% PAL of O<sub>2</sub>, tidally locked PC b, substellar point over Pacific Ocean |
| [PCb 1% PAL](/Tidally_locked_exoplanets/cases/PC_b_SSPO_1pc)  | GJ551 MUSCLES spectrum, 1% PAL of O<sub>2</sub>, tidally locked PC b, substellar point over Pacific Ocean |
| [PCb 0.1% PAL](/Tidally_locked_exoplanets/cases/PC_b_SSPO_01pc) | GJ551 MUSCLES spectrum, 0.1%PAL of O<sub>2</sub>, tidally locked PC b, substellar point over Pacific Ocean |

The tidally locked simulations cover several cases for known M dwarf terrestrial exoplanets. The basic modifications include changing the rotation rate of the model, changing the radius and gravitational acceleration, and fixing the solar zenith angle (in order to fix the substellar point). The substellar point is placed either in the middle of the Pacific ocean at 180&deg; longitude, or at 30&deg; longitude in Africa. Additionally, depending on the exoplanet in question, the solar file will need to be changed and scaled to the irradiance that the planet recieves. The Stellar Wind and Irradiance Model ([SWIM](https://github.com/jack-eddy-symposium/exoplanetary-impact/tree/main/SWIM)) has been developed for this purpose. The TRAPPIST-1e (TP-1e) simulations were originally started in 2020. In 2021, the Mega-MUSCLES spectrum (W21) became available, so this spectrum was introduced in order to compare the differences that arise in the simulations between the two spectra. 
