	SUBROUTINE ATTEN
	INCLUDE 'SOLDIV.FI'
	DELRAD = AMINOR/55.
	Z3 = XNSOL*SVTOTSOL*(DELN/EPDIV) + XNBAR*SVTOTBAR*DELTB
 	Z4 = Z3 + XNC(56)*SVTOTPED*0.5*DELRAD
      CALL EXPINT(Z3,E13,E23,E33,E43)
	CALL EXPINT(Z4,E14,E24,E34,E44)   
      FNCOLD(56)=(GAMOUTSOL/XNC(56))*(E23-E24)/(XNC(56)*SVTOTPED*DELRAD)
	IZAP = 0
      DO 100 I = 1,55
	IF(IZAP.EQ.1) GOTO 101
	J = 56-I
	E23 = E24
	Z4 = Z4 + XNC(J)*SVTOTPED*DELRAD
	CALL EXPINT(Z4,E14,E24,E34,E44) 
	FNCOLD(J)=(GAMOUTSOL/XNC(J))*(E23-E24)/(XNC(J)*SVTOTPED*DELRAD)
	IF(FNCOLD(J).LT.1.E-9) IZAP = 1
100	CONTINUE
101	CONTINUE 
	RETURN
	END

		
