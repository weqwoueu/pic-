MODULE Field_2D

IMPLICIT NONE

! Field Vectors
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	efx, efy, efz, bfx, bfy, bfz
REAL(8), DIMENSION(:,:), ALLOCATABLE	::	bt !$ ab.ZWZ

! Normalization
REAL(8)										c	! Simulation speed of light
REAL(8), DIMENSION(:), ALLOCATABLE		::	affp, affp_cell

! Static External Field (Uniform)
REAL(8)										efx_ext, efy_ext, efz_ext
REAL(8)										bfx_ext, bfy_ext, bfz_ext

! Potential and Charge for ES Code
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	phi, phibkgd, rho, rhoback, rhobeam
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	rho_s
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	rho_n !$ ab.ZWZ
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	rho_p !$ ab.ZWZ

! Kinectic Energy 
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	Ek_s !$ ab.ZWZ
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	Ek_tot !$ ab.ZWZ


! Current Density for ES Code (Diagnostic Purpose)
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	vj_sx, vj_sy, vj_sz    !!!!!电流密度
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	vj_totx, vj_toty, vj_totz
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	vx_s, vy_s, vz_s          !!!!!节点速度

REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE	::	vj_s1, vj_s2
REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE	::	t_s

! Reference 
! Density Reference for Particle Species
REAL(8), DIMENSION(:), ALLOCATABLE		::	dens0
!!! ************************ bjw add for impic 2019-6-3 **********************************************
REAL(8), DIMENSION(:), ALLOCATABLE		::	theta0    !!! 入射角度

! Density Reference for Background and Analytical Beam Ions  
REAL(8)										ainfty, a_beam0, a_beam0_o

! Density Potential Reference for Boltzman Electron Hybrid
REAL(8)										den0_ref, Te_ref, phi0_ref
REAL(8)										den0_left, den0_right, den0_right0, Te_left, Te_right, phi0_left, phi0_right

! Quasi-neutral Zone
REAL(8)										f_upqn, f_downqn, Ie_x,Ie_y,Ie_z 
INTEGER										k_upqn, k_downqn

REAL(8), DIMENSION(:,:), ALLOCATABLE		::	Aver_phi, Aver_rho
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	Aver_rho_s
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	Aver_Ek_s       !$ ab.ZWZ
REAL(8), DIMENSION(:,:,:), ALLOCATABLE	::	Aver_Ek_tot     !$ ab.ZWZ

END MODULE

