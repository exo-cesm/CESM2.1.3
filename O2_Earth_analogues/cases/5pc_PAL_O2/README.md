# 5% PAL case setup instructions

## If setting the case up on ARC4:

run buildcase_5pc_O2_ARC4. A case will be created with the name b.e21.BWma1850.f19_g17.5pc_o2.my_case.001. Change the name in buildcase_5pc_O2_ARC4 if you would like a different name.

The restart file for the 5% PAL case from the Cooke et al 2022 Royal Society Open Science paper (https://doi.org/10.1098/rsos.211165) is currently at /nobackup/Alternative_Earths/restarts/Five_pc_O2/Cooke_et_al_2022/

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, n the case directory, ./xmlchange the project to be blank (no input) and the job queue to be 40core-192G.q.

A user_nl_cam example will be copied in to the case directory.

The following modified files will be placed in SourceMods/src.cam/:
chemistry.F90; mo_tgcm_ubc.F90; upper_bc.F90; mo_jshort.F90; mo_photo.F90

These source mods are important for modifying the upper boundary conditions, and for updating the absorption in the Schumann-Runge bands. Note that Cooke et al. 2022 did not have the updated absorption in the Schumann-Runge bands.

## If starting from scratch:

The resolution used was 2.5 by 1.875 degrees (f19_g17)

Example case setup command:
./create_newcase --compset BWma1850 --res f19_g17 --case ~/b.e21.BWma1850.f19_g17.5pc_o2.001 --run-unsupported

The case should be run out for 5 days. Then, the surface mixing ratio of O<sub>2</sub> should be scaled to 0.0105. This can be done using the following command on the restart file and initial condition file (\*cam.r.\* and \*cam.i.\*).

ncap2 -O -s "O2=O2\*0.05" $infile $outfile

Ensure that the following modified files are placed in SourceMods/src.cam/:
chemistry.F90; mo_tgcm_ubc.F90; upper_bc.F90; mo_jshort.F90; mo_photo.F90


