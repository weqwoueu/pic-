MODULE IMPIC_Data_2D

!!! ************************ bjw add for impic 2019-6-3 **********************************************
    
IMPLICIT NONE

!!!! BJW ADD for implicit PIC
!REAL(8), DIMENSION(:,:), ALLOCATABLE       :: A_bar, A_bar_n_minus_1, A_bar_n, A_n_plus_1
REAL(8), DIMENSION(3)       :: A_bar, A_bar_n_minus_1, A_bar_n, A_n_plus_1
REAL(8), DIMENSION(:,:), ALLOCATABLE       :: Chi, OneAndChi
!REAL(8), DIMENSION(:,:), ALLOCATABLE       :: A_n_plus_1

LOGICAL      ::  IMPIC_index

!!!! BJW add for Magnetic(B)
LOGICAL      ::  Bfiled_index
REAL(8), DIMENSION(:,:,:), ALLOCATABLE     ::  Bfield   !! Bfield(1:3,i,j)
REAL(8)      :: Omega(3), OmegaT
REAL(8)      :: TransB(3,3)
REAL(8), DIMENSION(:,:,:), ALLOCATABLE     ::  TransChi   !! Bfield(1:9,i,j)


END MODULE