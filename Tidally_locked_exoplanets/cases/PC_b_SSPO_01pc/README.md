# Proxima Centauri b (0.1% PAL, with the substellar point over ocean) case setup instructions

## If setting the case up on ARC4:

run buildcase_PC_b_SSPO_01pc. A case will be created with the name b.e21.BWma1850.f19_g17.PC_b_Ocn.0.1pc_o2.my_case.001. Change the name in buildcase_PC_b_SSPO_01pc if you would like a different name.

The restart file for this case is currently at /nobackup/Alternative_Earths/restarts/PCb/01pc_PAL/0315-09-16/

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, in the case directory, ./xmlchange the project to be blank (no input) and the job queue to be 40core-192G.q.

A user_nl_cam example will be copied in to the case directory. This contains a changed solar file, based on the [Mega-MUSCLES](https://archive.stsci.edu/prepds/muscles/) Proxima Centauri b spectrum (see GJ551). It also specifies a different rotation rate, gravity, and radius.

Included source mods fix the solar zenith angle so the planet is considered to be tidally locked. New physics for absorption by carbon dioxide and water vapour in the Schumann-Runge bands are included. 
