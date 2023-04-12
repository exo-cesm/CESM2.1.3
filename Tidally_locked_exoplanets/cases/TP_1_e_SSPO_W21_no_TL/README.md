# Proxima Centauri b (with the substellar point over ocean) case setup instructions

## If setting the case up on ARC4:

run buildcase_PC_b_SSPO. A case will be created with the name b.e21.BWma1850.f19_g17.PC_b_Ocn.my_case.001. Change the name in buildcase_PC_b_SSPO if you would like a different name.

The restart file for this case is currently at /nobackup/Alternative_Earths/restarts/PC_b_ocn/0027-06-04/

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, n the case directory, ./xmlchange the project to be blank (no input) and the job queue to be 40core-192G.q.

A user_nl_cam example will be copied in to the case directory. This contains a changed solar file, based on the [Mega-MUSCLES](https://archive.stsci.edu/prepds/muscles/) Proxima Centauri b spectrum (see GJ551). It also specifies a different rotation rate, gravity, and radius.

Included source mods fix the solar zenith angle so the planet is considered to be tidally locked. The src.share source code modification ensures the ocean rotates at the same rate as the atmosphere and the planetary radius and gravity are consistent with that in user_nl_cam.
