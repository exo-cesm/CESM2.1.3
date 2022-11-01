# TRAPPIST-1 e case (with the substellar point over ocean) case setup instructions
## Uses the [Peacock et al. 2019](https://doi.org/10.3847/1538-4357/aaf891) model 1A version 1 spectrum 

## If setting the case up on ARC4:

run buildcase_TP_1_e_SSPO_P19. A case will be created with the name =b.e21.BWma1850.f19_g17.TP1_e_SSPO_P19.my_case.001. Change the name in buildcase_TP_1_e_SSPO_P19O if you would like a different name.

The restart file for this case is currently at /nobackup/Alternative_Earths/restarts/TP1_e_ocn/0173-02-25/

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, in the case directory, perform the ./xmlchange command and change the project to be blank (no input) and the job queue to be 40core-192G.q:

_./xmlchange PROJECT=_ 

_./xmlchange JOB\_QUEUE=40core-192G.q_

A user_nl_cam example will be copied in to the case directory. This contains a changed solar file, based on the [Peacock et al. 2019](https://doi.org/10.3847/1538-4357/aaf891) model 1A version 1 TRAPPIST-1 spectrum. It has then been scaled to TRAPPIST-1 e. It also specifies a different rotation rate, gravity, and radius.

Included source mods fix the solar zenith angle so the planet is considered to be tidally locked. The src.share source code modification ensures the ocean rotates at the same rate as the atmosphere and the planetary radius and gravity are consistent with that in user_nl_cam.
