# 1% PAL Young Sun (2Ga) with 21 hr rotation rate case setup instructions

## If setting the case up on ARC4:

run buildcase\_1pc\_O2\_2Ga\_YS\_21hr\_ARC4. A case will be created with the name . Change the name in buildcase\_1pc\_O2\_2Ga\_YS\_21hr\_ARC4 if you would like a different case name.

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, in the case directory, ./xmlchange the project to be blank (no input) and the job queue to be 40core-192G.q.

A user\_nl\_cam example will be copied in to the case directory.

The following modified files will be placed in SourceMods/src.cam/:
chemistry.F90; mo\_tgcm\_ubc.F90; upper\_bc.F90; mo\_jshort.F90; mo\_photo.F90

These source mods are important for modifying the upper boundary conditions, and for updating the absorption in the Schumann-Runge bands.

## If setting the case from scratch

The resolution used was 2.5 by 1.875 degrees (f19_g17)

Example case setup command:
./create_newcase --compset BWma1850 --res f19_g17 --case ~/b.e21.BWma1850.f19_g17.1pc_o2.2Ga_Young_Sun.001 --run-unsupported

The case should be run out for 5 days. Then, the surface mixing ratio of O<sub>2</sub> should be scaled to 0.0021. This can be done using the following command on the restart file and initial condition file (\*cam.r.\* and \*cam.i.\*).

ncap2 -O -s "O2=O2\*0.01" $infile $outfile

Ensure that the following modified files are placed in SourceMods/src.cam/:
chemistry.F90; mo_tgcm_ubc.F90; upper_bc.F90

The 2Ga Young Sun solar file (SolarForcingCMIP6piControl_c160921_2.0Ga.nc) also needs to be pointed to in user_nl_cam
