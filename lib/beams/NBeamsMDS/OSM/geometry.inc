c----------------------FILE geometry.inc------------------------------

c/    This file contains declarations for the variables that 
c/    characterize the straightline geometry in the calculation
c/    of the transmission probabilities.

      integer inu
      real L_i, L_j, L_perp, l_mfp, theta_ij, alpha_j, delta
      common /rcgeom/ inu, L_i, L_j, L_perp, l_mfp, theta_ij,
     .   alpha_j, delta
