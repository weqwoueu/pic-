MODULE MCC_Main_Param
    
IMPLICIT NONE

INTEGER(4), PARAMETER                 	:: energy_max = 500
REAL(8), PARAMETER	                    :: co_bohm=1./64.,co_bohm1=1./64.

REAL(8)                                 :: SV_MAX_E
REAL(8)                                 :: SELS(500),SEXC(500),SION(500)
REAL(8)                                 :: iatom

INTEGER(4), PARAMETER                 	:: Mergenum=50

REAL(8), PARAMETER                 	    :: atom_flux=3e-6       !原子质量流量3mg

REAL(8), PARAMETER                      :: vaz=200     !!原子的定向速度

!用于原子当流体时或原子当背景时给定原子分布
!REAL(8), DIMENSION(:,:), ALLOCATABLE	::	gden

!REAL(8), DIMENSION(:,:), ALLOCATABLE	::	bohm,elastic,excite,ionize



	!REAL, PARAMETER       :: mfactor=1.
!	REAL, PARAMETER       :: gfactor=100.
!	REAL, PARAMETER       :: pfactor=1600
    REAL(8), PARAMETER       :: mul=240. !2015/10/22 LC
	REAL(8), PARAMETER       :: pfactor=1 !2014/12/7 LC
	REAL(8), PARAMETER       :: sfactor=1                !taccogna方法：减小推力器尺寸
	!REAL, PARAMETER       :: kb=1.38E-23              !Boltz constant
	REAL(8), PARAMETER       :: EPSILON0=8.85E-12        !Vaccum permitivity
!	REAL, PARAMETER       :: EPSILON0=1        !Vaccum permitivity
!	REAL, PARAMETER       :: E=1.6022E-19             !Electron charge
!	REAL, PARAMETER       :: E=1             !Electron charge
!	REAL, PARAMETER       :: Me=9.1094E-31 !*gfactor            !Electron mass
!!	REAL, PARAMETER       :: Me=1            !Electron mass
!	REAL, PARAMETER       :: Mi=2.19E-25*mfactor      !Xe ion mass
!	REAL, PARAMETER       :: NERO=4.4169657E17/sfactor                !number density reference1E18
	REAL(8), PARAMETER       :: NERO=1.0E17   !!5.12E14  !!!4.0E17/sfactor !2014/12/7 LC
!	REAL, PARAMETER       :: NERO=1                !number density reference
!	REAL, PARAMETER       :: Fnum=1E9                !number density reference

!!! bjw
	!REAL, PARAMETER       :: T_Ref=5                  !Temperature reference eV
	!REAL, PARAMETER       :: Phi_Ref=5                !Voltage reference V

	REAL(8), PARAMETER       :: t_wall=0.026             !eV   300K
	
	REAL(8), PARAMETER       :: phiwall=0

    REAL(8)                    omigp       !=1E-11
    REAL(8)                    omigc 
!	REAL(8), PARAMETER	  :: B_Ref=1.776663259E-3  !2014/12/7 LC    
    !REAL                    B_Ref
    REAL(8)                   LMDD 
    REAL(8)                    Vel_Ref 
    REAL(8)                    EPSILON_Ref
!    REAL                    Fnum 2014/12/7 LC




END MODULE