SUBROUTINE FUN_sphere_u_ErrAna(	x, y, region,	&
							dind_x, dind_y, f)

!USE Objects
USE IFE_Data
USE IFE_RHS_PARAM

IMPLICIT NONE

REAL(8), DIMENSION(:), INTENT(IN)	::	x, y !, z
INTEGER, INTENT(IN)					::	region, dind_x, dind_y !, dind_z
REAL(8), DIMENSION(:), INTENT(OUT)	::	f


REAL(8), DIMENSION(:), ALLOCATABLE	::	r
REAL(8)									beta_minus, beta_plus



beta_minus = Global_Beta(1);
beta_plus  = Global_Beta(2);

ALLOCATE(r(SIZE(x)))

IF (region/=-1) THEN
	IF		(dind_x==1) THEN

		f = 5 * x * (x**2 + y**2)**1.5
	ELSEIF	(dind_y==1) THEN

		f = 5 * y * (x**2 + y**2)**1.5

	ELSE

		f = (x**2 + y**2)**2.5		
	END IF
ELSEIF (region==-1) THEN
	IF		(dind_x==1) THEN

		f = 5 * x * (x**2 + y**2)**1.5

	ELSEIF	(dind_y==1) THEN

		f = 5 * y * (x**2 + y**2)**1.5	

	ELSE

		f = (x**2 + y**2)**2.5
        
	END IF
END IF

DEALLOCATE(r)

END