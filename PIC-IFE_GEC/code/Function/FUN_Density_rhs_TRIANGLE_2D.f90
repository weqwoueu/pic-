SUBROUTINE FUN_Density_rhs_TRIANGLE_2D(rho_, phi_, region, f_, x, y, n)

USE IFE_MAIN_PARAM
USE IFE_RHS_PARAM
USE Object_Data_2D
USE IFE_Data

IMPLICIT NONE

INTEGER									n
REAL(8), DIMENSION(n), INTENT(IN)	::	rho_, phi_
INTEGER, INTENT(IN)					::	region
REAL(8), DIMENSION(n), INTENT(OUT)	::	f_
REAL(8), DIMENSION(n), INTENT(IN)	::	x, y

INTEGER	i
REAL(8)                                 Xi, Yi

REAL(8)									beta_minus, beta_plus, y0
!$ ab.ZWZ for analysing error convergence ======== \\
REAL(8) :: a11, a12, a21, a22, z, r, r_1
!$ ab.ZWZ for analysing error convergence ======== //

IF (n /=3) THEN
PRINT*, ' n/=3, Check FUN_Density_rhs_2D.f90, STOP'
STOP
ENDIF

beta_plus = Global_Beta(1);
IF (SIZE(Global_Beta)>1) beta_minus  = Global_Beta(2);

y0 = 1.0

!IF (region==vacuum) THEN
!	DO i=1,SIZE(phi_)
!		IF ( phi_(i)<=phi_bkgd(-vacuum) ) THEN 
!			f_(i) = rho_(i) - Rho_bkgd(-vacuum) *	&
!					DEXP( (phi_(i)-phi_bkgd(-vacuum))/Te_bkgd(-vacuum) )
!		ELSE
!			f_(i) = rho_(i) - Rho_bkgd(-vacuum) *	&
!					( One+(phi_(i)-phi_bkgd(-vacuum))/Te_bkgd(-vacuum) )
!		END IF
!	END DO
!ELSE
!	DO i=1,SIZE(phi_)
!		IF ( phi_(i)<=phi_bkgd(-vacuum) ) THEN 
!			f_(i) = Zero
!		ELSE
!			f_(i) = Zero
!		END IF
!	END DO
!ENDIF

!DO i=1,SIZE(phi_)
!	f_(i) = -4.0
!!	f_(i) = -2.0
!!	f_(i) = -10.0
!!	f_(i) =  0.0
!ENDDO


!IF (region==vacuum) THEN
!	DO i=1,SIZE(phi_)
!!		f_(i) = -4.0/beta_minus
!!		f_(i) = -4.0 + (x(i)*y(i) - x(i)*y(i))
!!	    IF (x(i)==0) THEN
!!		   Xi = SmallValue*Smallvalue
!!		ELSE
!!		   Xi = x(i)
!!		ENDIF
!	    IF (y(i)==0) THEN
!		   Yi = SmallValue
!		ELSE
!		   Yi = y(i)
!		ENDIF
!!		f_(i) = -6.0 + 2*y0/Yi
!
!		f_(i) = exp(phi_(i)) - exp(x(i) + y(i))
!
!	END DO
!ELSE
!	DO i=1,SIZE(phi_)
!!		f_(i) = -4.0/beta_plus
!!		f_(i) = -4.0 + (x(i)*y(i) - x(i)*y(i))
!	    IF (y(i)==0) THEN
!		   Yi = SmallValue
!		ELSE
!		   Yi = y(i)
!		ENDIF
!!		f_(i) = -6.0 + 2*y0/Yi
!
!		f_(i) = exp(phi_(i)) - exp(x(i) + y(i))
!	END DO
!
!ENDIF
!
!DO i=1,SIZE(phi_)
!	f_(i) =  0.0
!ENDDO

!$ ab.ZWZ for analysing error convergence ======== \\
!a11 = 2.
!a12 = 2.
!a21 = 2.
!a22 = 5.
!DO i=1,SIZE(phi_)
!    !f_(i) = -4.0*(1+x(i)*x(i)+y(i)*y(i))*EXP(x(i)*x(i)+y(i)*y(i))                                  !$ Cartesian  isotropy
!    f_(i) = -4.0*(1+(y(i)-1)/2/y(i)+x(i)*x(i)+(y(i)-1)*(y(i)-1))*EXP(x(i)*x(i)+(y(i)-1)*(y(i)-1))  !$ Cylindrical isotropy
!    !f_(i) = -(14+8*x(i)**2+16.0*x(i)*y(i)+20*y(i)**2)*EXP(x(i)**2+y(i)**2)                           !$ Cartesian  anisotropy
!    !z = x(i)
!    !r = y(i)
!    !r_1 = r - 1 
!    !f_(i) = -( a11*(2+4.*z**2) + a12*(4.*z*r_1) + a21*(2.*z/r+4.*z*r_1) + a22*(2.*r_1/r+2+4.*r_1**2) )*EXP(z**2+r_1**2) !$ Cylindrical  anisotropy
!END DO
!$ ab.ZWZ for analysing error convergence ======== //

IF (region == -1) THEN
	DO i=1,SIZE(phi_)
!        f_(i) = rho_(i) +(x(i)**2)+(y(i)**2 +(0.5-1)*(0.61**2))-phi_(i)
		f_(i) = rho_(i)
        !f_(i) = -4
!        f_(i) = rho_(i)+(x(i)**2)+(y(i)**2)-phi_(i)
!		f_(i) = 0
!        f_(i) = rho_(i)+(x(i)**2)+(y(i)**2 +(0.1-1)*(0.61**2))-phi_(i)
!		f_(i) = Zero
	END DO
ELSE
	DO i=1,SIZE(phi_)
!		f_(i) = Zero
!		f_(i) = rho_(i)
!        f_(i) = rho_(i)+(x(i)**2)+(y(i)**2)/2-phi_(i)
!        f_(i) = rho_(i)+(x(i)**2)+(y(i)**2)-phi_(i)
		f_(i) = rho_(i)
        !f_(i) = -4
!        f_(i) = rho_(i)+(x(i)**2)+(y(i)**2)/10-phi_(i)
	END DO
ENDIF



END