# TRAPPIST-1 e (P19 spectrum with the substellar point over ocean) case setup instructions

## If setting the case up on ARC4:

run buildcase_TP_1_e_SSPO_P19. A case will be created with the name b.e21.BWma1850.f19_g17.TRAPPIST1_e_P19_PI.my_case.001. Change the name in buildcase_TP_1_e_SSPO_P19 if you would like a different name.

The restart file for this case is currently at /nobackup/Alternative_Earths/restarts/TP1_e/P19_SSPO_PI/0335-08-25/

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, n the case directory, ./xmlchange the project to be blank (no input) and the job queue to be 40core-192G.q.

A user_nl_cam example will be copied in to the case directory. This contains a changed solar file, based on the [Peacock et al. (2019)](https://archive.stsci.edu/hlsp/hazmat) TRAPPIST-1 spectrum. It also specifies a different rotation rate, gravity, and radius.

Included source mods fix the solar zenith angle so the planet is considered to be tidally locked. 

