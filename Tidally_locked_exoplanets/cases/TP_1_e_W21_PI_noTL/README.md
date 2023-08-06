# TRAPPIST-1 e (W21 spectrum, not tidally locked, 1 day rotation rate) case setup instructions

## If setting the case up on ARC4:

run buildcase_TP_1_e_SSPO_W21_no_TL. A case will be created with the name b.e21.BWma1850.f19_g17.TRAPPIST1_e.no_T_lock.W21.my_case.001. Change the name in buildcase_TP_1_e_SSPO_W21_no_TL if you would like a different name.

The restart file for this case is currently at /nobackup/Alternative_Earths/restarts/TP1_e/W21_SSPO_PI_noTL/0293-01-10/

A case will be created and the project will be planet with the job queue planet.q. If switching to the main queue, in the case directory, ./xmlchange the project to be blank (no input) and the job queue to be 40core-192G.q.

A user_nl_cam example will be copied in to the case directory. This contains a changed solar file, based on the [Wilson et al. (2021)](https://zenodo.org/record/4556130#.Y_82yezP39E) TRAPPIST-1 spectrum. It also specifies a different rotation rate, gravity, and radius.

Included source mods fix the solar zenith angle so the planet is considered to be tidally locked. 

