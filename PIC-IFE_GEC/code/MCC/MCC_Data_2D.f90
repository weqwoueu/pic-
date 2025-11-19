MODULE MCC_Data_2D
   
!! Jinwei Bai 
!! Purpose:		define the parameter structrue for MCC. 
!! Last Update:	2019-3-10     
    
IMPLICIT NONE    

!IONIZATION RATE
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	ionize

!bohm collision rate
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	bohm

!elastic collision rate
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	elastic

!excite collision rate
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	excite

!ion recombinate to atom
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	recom

!wall conduct rate
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	cond
!
!!Capacitance
!REAL(8), DIMENSION(:,:), ALLOCATABLE		::	capa

!atom density for fluid model!用于原子当流体时或原子当背景时给定原子分布

REAL(8), DIMENSION(:,:), ALLOCATABLE		::	gden

REAL(8), DIMENSION(:), ALLOCATABLE		::	ENER  !粒子的真实能量


!!!! bjw add ,查表法得到截面
REAL(8), DIMENSION(:,:), ALLOCATABLE		::	section_ela, section_exc, section_ion 

END MODULE