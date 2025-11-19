MODULE Constant_Variable_2D

IMPLICIT NONE

!!! ************************ bjw add for impic 2019-6-3 **********************************************
	REAL(8), PARAMETER       :: mfactor=1
	REAL(8), PARAMETER       :: pfactor=1600
	REAL(8), PARAMETER       :: kb=1.3807d-23             !Boltz constant
	REAL(8), PARAMETER       :: EPSILON0=8.8542d-12         !Vaccum permitivity 单位：F/m
	REAL(8), PARAMETER       :: E=1.6022d-19             !Electron charge
	REAL(8), PARAMETER       :: Me=9.1095d-31            !Electron mass
    !REAL(8), PARAMETER       :: MI=2.19E-25    !Xe ion mass
    REAL(8), PARAMETER       :: MI=6.67E-27     !He ion mass
    !REAL(8), PARAMETER       :: MI=6.64E-26     !Ar ion mass(kg)
	!REAL(8), PARAMETER       :: NERO=1                !number density reference
    
    !!! 无量纲参数
    REAL(8)  :: m_ref
    REAL(8)  :: q_ref
    REAL(8)  :: epsilon_ref
    REAL(8)  :: n_ref
    REAL(8)  :: T_ref
    REAL(8)  :: Phi_ref
    REAL(8)  :: Efield_ref
    REAL(8)  :: Bfield_ref
    REAL(8)  :: L_ref
    REAL(8)  :: time_ref
    REAL(8)  :: v_ref

    REAL(8)  :: Density_ref          !!! 密度参考值(m-3)
    REAL(8)  :: Temperature_ref      !!! 温度参考值(eV)
    REAL(8), PARAMETER :: eVtoK=11605.d0            !!!! 1eV = 11605 K  (1eV=kT,k为玻尔兹曼常数)
    INTEGER(4)                              ::  ParticlePerGrid   !!! 每个网格内的粒子数（电子）
    REAL(8)                                 ::  affp_bjw(3)           !!! affp(ipf(jj)):每种粒子的模拟粒子数代表的真实粒子数
!!! ************************ bjw add for impic 2019-6-3 **********************************************    
    
	!REAL, PARAMETER       :: mfactor=1
	!REAL, PARAMETER       :: pfactor=1600
	!REAL, PARAMETER       :: kb=1.38E-23              !Boltz constant
	!REAL, PARAMETER       :: EPSILON0=1        !Vaccum permitivity
	!REAL, PARAMETER       :: E=1             !Electron charge
	!REAL, PARAMETER       :: Me=1            !Electron mass
	!REAL, PARAMETER       :: MI=2.19E-25*mfactor      !Xe ion mass
	!REAL, PARAMETER       :: NERO=1                !number density reference
	!REAL, PARAMETER       :: T_Ref=1                  !Temperature reference eV
	!REAL, PARAMETER       :: Phi_Ref=1                !Voltage reference V
    
    
END MODULE