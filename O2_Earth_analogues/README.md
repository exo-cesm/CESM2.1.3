# 1. O<sub>2</sub> Earth-analogue and paleoclimate case descriptions

## Publications

The simulations in section 1 were used in the following publications:

Cooke GJ, Marsh DR, Walsh C, Black B, Lamarque J-F. 2022, A revised lower estimate of ozone columns during Earth’s oxygenated history. R. Soc. Open Sci. 9: 211165. https://doi.org/10.1098/rsos.211165.

G J Cooke, D R Marsh, C Walsh, S Rugheimer, G L Villanueva, Variability due to climate and chemistry in observations of oxygenated Earth-analogue exoplanets, Monthly Notices of the Royal Astronomical Society, 2022;, stac2604, https://doi.org/10.1093/mnras/stac2604.

## Simulations

We define O<sub>2</sub> analogues as Earth-like planets that have an oxygenated atmosphere. Simulations between 1000 times less than the present atmospheric level (PAL) of oxygen (the present atmospheric level of O<sub>2</sub> is 21% by volume) and up to 1.5 times PAL are included in this definition.

The simulations that have been performed are detailed in the table below:

<ins>**Table 1**</ins>

| Simulation name   | O<sub>2</sub> mixing ratio [PAL] | CH<sub>4</sub> lower boundary        | CO<sub>2</sub> lower boundary | Solar spectrum |
| ---------------   | ---------------------- | -------------------------- | ------------------- | -------------- |
| Pre-industrial    | 1.000                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 150% PAL          | 1.500                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 50% PAL           | 0.500                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 10% PAL           | 0.100                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 5% PAL            | 0.050                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL            | 0.010                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 0.5% PAL          | 0.005                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 0.1% PAL          | 0.001                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL CH4 em10   | 0.010                  | Flux of 10x PD emissions   | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL CH4 em1    | 0.010                  | Flux of 1x PD emissions    | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL CH4 em0.1  | 0.010                  | Flux of 0.1x PD emissions  | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL YS         | 0.010                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | 2 Ga Sun        |
| 1% PAL YS 4xCO2   | 0.010                  | Fixed at 0.8 ppmv          | Fixed at 1120 ppmv  | 2 Ga Sun        |

## Branching from run restarts

Currently, the restart files are stored on ARC4 at /nobackup/Alternative_Earths/restarts/. Permanent storage will be in place eventually.

To branch from one of these restart files, you can do it manually, or use the buildcase scripts in each case folder. 

## Fixed lower boundary modifications

The standard globally averaged lower boundary conditions are given in the following table:

<ins>**Table 2**</ins>

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

Boundary conditions can be changed. For example, the lower boundary for O<sub>2</sub>  and CO<sub>2</sub> were scaled using the following example operations on the \*cam.r.\* (restart) file and the  \*cam.i.\* (initial condition)file for the pre-industrial atmosphere:

ncap2 -O -s "CO2=CO2\*4.0" $infile $outfile #this multiplies the CO<sub>2</sub> field by 4
ncap2 -O -s "O2=O2\*0.1" $infile $outfile #this multiplies the O<sub>2</sub> field by 0.1

Adittionally, the flbc_file needed to be modified in user_nl_cam to ensure the correct fixed mixing ratio lower boundary condition was applied:
&chem_surfvals_nl
 flbc_file              = '/glade/scratch/gregcooke/LBC_files/LBC_4xCO2.nc'
/
For the simulations where O<sub>2</sub> was at 0.5% PAL or lower, O<sub>2</sub> needed to be inlcuded in the lower boundary file and fixed at the surface. The total pressure was not changed for any of the simulations. N<sub>2</sub> increases when O<sub>2</sub> decreases, for example, such that each simulation has the same atmospheric surface pressure.
For O<sub>2</sub> flbc_modifications in user_nl_cam:
&chem_surfvals_nl
 flbc_list       = 'CH3BR', 'CH3CL', 'CH4', 'O2',
         'CO2', 'H2', 'N2O', 'OCS', 'CFC11','CFC12','CFC11eq'
flbc_file = '/glade/scratch/gregcooke/LBC_files/LBC_0.001xO2.nc'
/

## Methane lower boundary

Changes to methane were made in user_nl_cam by removing it from flbc_list and specifing it in an emissions file:

&chem_inparm
srf_emis_specifier = 'CH4 -> /glade/scratch/gregcooke/LBC_files/CH4_1_PD_emissions.nc'
/
&chem_surfvals_nl
flbc_list              = 'CCL4', 'CF2CLBR', 'CF3BR', 'CFC11', 'CFC113', 'CFC12', 'CH3BR', 'CH3CCL3', 'CH3CL',
         'CO2', 'H2', 'HCFC22', 'N2O', 'CFC114', 'CFC115', 'HCFC141B', 'HCFC142B', 'CH2BR2', 'CHBR3',
         'H2402', 'OCS', 'CFC11eq'
/

## Solar spectrum changes

The standard solar spectrum in WACCM6 is SolarForcingCMIP6piControl_c160921.nc

We modified this solar file by using the Evolution Of Solar Flux IDL code from https://depts.washington.edu/naivpl/content/models.
This model is described in Claire et al. 2012: http://iopscience.iop.org/article/10.1088/0004-637X/757/1/95/meta

In the case directory, user_nl_cam was modified to include the new solar file for the Sun 2 Gyr ago:
&solar_data_opts
solar_irrad_data_file = "/glade/u/home/gregcooke/SolarForcingCMIP6piControl_c160921_2.0Ga.nc"
/

## Location of restart files

Restart files are currently located at /nobackup/Alternative\_Earths/restarts/ on the ARC4 supercomputer at the University of Leeds.

# 2. Since [Cooke et al. 2022a](https://doi.org/10.1098/rsos.211165), additional simulations have been perfomed

Perturbations have been tested on the 10% PAL and 1% PAL simulations:

<ins>**Table 3**</ins>

| Simulation name   | O<sub>2</sub> mixing ratio [PAL] | CH<sub>4</sub> lower boundary        | CO<sub>2</sub> lower boundary | Solar spectrum | Lower boundary conditions |
| ---------------   | ---------------------- | -------------------------- | ------------------- | -------------- | -----------------|
| 10% PAL CH4 em10   | 0.100                  | Flux of 10x PD emissions    | Fixed at 280 ppmv   | Modern day Sun  | N/A |
| 10% PAL CH4 em1 | 0.100                  | Flux of 1x PD emissions  | Fixed at 280 ppmv   | Modern day Sun  | N/A |
| 10% PAL CH4 em0.1 | 0.100                  | Flux of 0.1x PD emissions  | Fixed at 280 ppmv   | Modern day Sun  | N/A |
| 10% PAL 0.1x CH3CL | 0.100                  | Fixed at 0.8 ppmv  | Fixed at 280 ppmv   | Modern day Sun  | f<sub>CH<sub>3</sub>Cl</sub> = 4.57e-11 |
| 1% PAL 0.1x CH3CL   | 0.100                  | Fixed at 0.8 ppmv  | Fixed at 280 ppmv   | Modern day Sun  | f<sub>CH<sub>3</sub>Cl</sub> = 4.57e-11 |
| 1% PAL 1e-9 All BR & CL sources | 0.100                  | Fixed at 0.8 ppmv  | Fixed at 280 ppmv   | Modern day Sun  | All Cl and Br sources in [Table 2](#Fixed-lower-boundary-modifications) reduced by a factor of 1 billion |

where f<sub>X</sub> deotes the mixing ratio of species X.

### Simulations with code updates

CO<sub>2</sub> and H<sub>2</sub>O absorption has now been added to the Schumann-Runge bands. We performed perturbations on the 1% PAL simulation by adding in the absorption, changing the lower boundary conditions, and then both changes at the same time.

<ins>**Table 4**</ins>

| Simulation name   | O<sub>2</sub> mixing ratio [PAL] | CH<sub>4</sub> lower boundary | CO<sub>2</sub> lower boundary | Solar spectrum |
| ---------------   | ---------------------- | -------------------------- | ------------------- | --------------- |
| 1% PAL SRB        | 0.010                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL LBC        | 0.010                  | See table below            | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL SRB + LBC  | 0.010                  | See table below            | Fixed at 280 ppmv   | Modern day Sun  |

For the simulations labelled with LBC, the following lower boundary conditions were used:

<ins>**Table 5**</ins>

| Chemical | Mixing ratio | Flux [molecules cm<sup>-3</sup>]       |
| -------- | ------------ | -------------------------------------- |
| CCL<sub>4</sub>               | 2.50e-14 | N/A                   |
| CF<sub>2</sub>CLBR            | 4.45e-15 | N/A                   |
| CFC11eq                       | 3.21e-11 | N/A                   |
| CH<sub>2</sub>BR<sub>2</sub>  | 1.20e-12 | N/A                   |
| CH<sub>3</sub>BR              | 5.30e-12 | N/A                   |
| CH<sub>3</sub>CL              | N/A      | 2.92e8                |
| CH<sub>4</sub>                | 8.08e-07 | N/A                   |
| CHBR<sub>3</sub>              | 1.20e-12 | N/A                   |
| CO<sub>2</sub>                | 0.000284 | N/A                   |
| H<sub>2</sub>                 | N/A      | 3e-4                  |
| OCS                           | 3.37e-10 | N/A                   |
| N<sub>2</sub>O                | N/A      | 1.01e9                |
| CH<sub>4</sub>                | N/A      | 1.14e11               |
| CO                            | N/A      | 2.21e11               |
| Total BROY                    | 1.13e-11 | N/A                   |
| Total CLOY                    | 4.57e-10 | N/A                   |
| Total FOY                     | 8.89e-15 | N/A                   |

