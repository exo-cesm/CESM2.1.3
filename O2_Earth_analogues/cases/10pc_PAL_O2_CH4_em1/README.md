# 10% PAL CH4 em1 case setup instructions

The resolution used was 2.5 by 1.875 degrees (f19_g17)

Example case setup command:
./create_newcase --compset BWma1850 --res f19_g17 --case ~/b.e21.BWma1850.f19_g17.10pc_o2.CH4_PD_flux.001 --run-unsupported

The case should be run out for 5 days. Then, the surface mixing ratio of O<sub>2</sub> should be scaled to 0.021. This can be done using the following command on the restart file and initial condition file (\*cam.r.\* and \*cam.i.\*).

ncap2 -O -s "O2=O2\*0.1" $infile $outfile

Ensure that the following modified files are placed in SourceMods/src.cam/:
chemistry.F90; mo_tgcm_ubc.F90; upper_bc.F90

See the example user_nl_cam for how to include the present day methane emission file: CH4_1_PD_emissions.nc

