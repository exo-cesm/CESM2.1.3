# O<sub>2</sub> analogues case descriptions
We define O<sub>2</sub> analogues as Earth-like planets that have an oxygenated atmosphere. Simulations between 1000 times less than the present atmospheric level (PAL) of oxygen (the present atmospheric level of O<sub>2</sub> is 21% by volume) and up to 1.5 times PAL are included in this definition.

The simulations that have been performed are detailed in the table below:

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
| 10% PAL CH4 em1   | 0.100                  | Flux of 1x PD emissions    | Fixed at 280 ppmv   | Modern day Sun  |
| 10% PAL CH4 em0.1 | 0.100                  | Flux of 0.1x PD emissions  | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL CH4 em10   | 0.010                  | Flux of 10x PD emissions   | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL CH4 em1    | 0.010                  | Flux of 1x PD emissions    | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL CH4 em0.1  | 0.010                  | Flux of 0.1x PD emissions  | Fixed at 280 ppmv   | Modern day Sun  |
| 1% PAL YS         | 0.010                  | Fixed at 0.8 ppmv          | Fixed at 280 ppmv   | 2 Ga Sun        |
| 1% PAL YS 4xCO2   | 0.010                  | Fixed at 0.8 ppmv          | Fixed at 1120 ppmv  | 2 Ga Sun        |


## Fixed lower boundary modifications

The standard globally averaged lower boundary conditions are given in the following table:

| Chemical | Mixing ratio |
| -------- | ------------ |
| CCL4 | 2.50e-14 |
| CF2CLBR | 4.45e-15 |
| CFC11eq | 3.21e-11 |
| CH2BR2 | 1.20e-12 |
| CH3BR | 5.30e-12 |
| CH3CL | 4.57e-10 |
| CH4 | 8.08e-07 |
| CHBR3 | 1.20e-12 |
| CO2 | 0.000284 |
| H2 | 5.00e-07 |
| OCS | 3.37e-10 |
| N2O | 2.73e-07 |
| Total BROY | 1.13e-11 |
| Total CLOY | 4.57e-10 |
| Total FOY | 8.89e-15 |


Boundary conditions can be cahnged. For example, the lower boundary for O<sub>2</sub>  and CO<sub>2</sub> were scaled using the following example operations on the \*cam.r.\* (restart) file and the  \*cam.i.\* (initial condition)file for the pre-industrial atmosphere:

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

## Publications

Cooke GJ, Marsh DR, Walsh C, Black B, Lamarque J-F. 2022, A revised lower estimate of ozone columns during Earth’s oxygenated history. R. Soc. Open Sci. 9: 211165. https://doi.org/10.1098/rsos.211165

## Location of restart files

Restart files are currently located at /glade/scratch/gregcooke/restarts/ on the NCAR Cheyenne supercomputer.
