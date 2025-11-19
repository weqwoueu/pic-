MODULE 	DSMC_Data_2D
	
IMPLICIT NONE

REAL(8), DIMENSION(:,:), ALLOCATABLE	::	CCG         !,CS

REAL(8), DIMENSION(:), ALLOCATABLE	::	CC          !CC is the real volume of the cell
INTEGER(4), DIMENSION(:), ALLOCATABLE	::	IR             

REAL(8)                                ::	SP(2)
    
REAL(8)                                ::	SPM(5)

REAL, PARAMETER                       ::	neutraltime=2.e1
INTEGER, PARAMETER                    ::	atomt=100
 

!SPM(1) Åö×²½ØÃæ
!SPM(2) the reference temperature 
!SPM(3) the viscosity temperature power law 
!SPM(4) the reciprocal of the VSS scattering parameter
!SPM(5) the gamma function of (5/2-viscosity-temperature power law)

END MODULE 
