	SUBROUTINE EXPINT(Z,E1,E2,E3,E4)
C		CALCULATES EXPONENTIAL INTEGRAL FUNCTIONS
	E0 = EXP(-1.*Z)/Z
      IF(Z.GT.1.)  GOTO 205
      E1 = -LOG(Z)-.57721566+.999999193*Z-.24991055*(Z**2)+.05519968*
     2     (Z**3)-.00976004*(Z**4)+.00107857*(Z**5)
      GOTO 206
205   E1 = (EXP(-Z)/Z)*(Z**2+2.334733*Z+.250621)/(Z**2+3.330657*Z+
     2     1.681534)
206   E2 = (EXP(-Z) - Z*E1)/1.
      E3 = (EXP(-Z) - Z*E2)/2.
      E4 = (EXP(-Z) - Z*E3)/3.
      E5 = (EXP(-Z) - Z*E4)/4.
      E6 = (EXP(-Z) - Z*E5)/5.
      E7 = (EXP(-Z) - Z*E6)/6.
	  RETURN
	  END
