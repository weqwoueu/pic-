!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: Input_2D.f90                                       C
!
!  Purpose: Input data and normalization
!                                                                      C
!  Author: Yuchuan Chu                                Date: 03-May-12  C
!  Comments:                                                           C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
SUBROUTINE Input_2D(xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects, bc_point_1, bc_point_2, bc_index, bc_value)

USE Object_Data_2D
USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
USE TimeControl
USE Object_2D
USE Constant_Variable_2D
USE IMPIC_Data_2D
IMPLICIT NONE

REAL(8)		                                    ::  xmin, xmax, ymin, ymax, zmin, zmax
INTEGER		                                    ::  nnx, nny, nnz, N_Objects, N_Boundary
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
REAL(8), DIMENSION(:,:), POINTER	                ::	bc_point_1, bc_point_2
INTEGER, DIMENSION(:), POINTER	                ::	bc_index
REAL(8), DIMENSION(:), POINTER		            ::	bc_value

INTEGER                                            ::  i, jj, iface, iobject, Idummy

INTEGER                                            ::  N_PIC_Boundary
!=====================================LY modification, 2022-6-13, Checking Objects====================================
Real(8) :: xDist, yDist
Real(8) :: Focusx1, Focusy1, Focusx2, Focusy2
Real(8) :: FocusDist
!=====================================LY modification, 2022-6-13, Checking Objects====================================

!==================Add for periodic boundary conditions===================
REAL(8),DIMENSION(:,:),	POINTER		  ::	ALQ
INTEGER num
!=========================================================================

!!!!!!!!!!!!!!!!!!!!Normalization!!!!!!!!!!!!!!!!!!!!!!!!!!!
!REAL(8)        ::  LMDD
!REAL(8)        ::  Vel_Ref
!REAL(8)        ::  Time_Ref                                
!
!!!!!!!!!!!!!!!!!!!!Initialize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!LMDD = SQRT(EPSILON0*kb*11600*T_Ref/(E*E*NERO))
!Vel_Ref = SQRT(kb*T_Ref*11600/Me)
!Time_Ref = LMDD / Vel_Ref
!
!!!!temp
!LMDD = 1
!Vel_Ref = 1
!Time_Ref = 1



!!!!!!!!!!!!!!!!Setup_IFE_Mesh!!!!!!!!!!!!!!!!!!!!!!
WRITE(6,*)
WRITE(6,*)'READ "mesh.inp"'
OPEN(1, ACTION='READ', FILE='./INPUT/mesh.inp')
READ(1,*) xmin, ymin, zmin
READ(1,*) xmax, ymax, zmax
READ(1,*) nnx, nny, nnz
READ(1,*) (hx(i), i = 1, 2)
READ(1,*) hz
!IF(hx(1) /= hx(2)) THEN
!    stop
!ENDIF
CLOSE(1)

!-------initialize the Vert_o-----
Vert_o(1) = xmin - hx(1)
Vert_o(2) = ymin - hx(2)
Vert_o(3) = zmin - hz
!-------------

!!! ************************ bjw add for impic 2019-6-3 **********************************************
nx = nnx
ny = nny
nz = nnz
!!! ************************ bjw add for impic 2019-6-3 **********************************************

WRITE(6,*)
WRITE(6,*)'READ "object.inp"'

OPEN(1, ACTION='READ', FILE='./INPUT/object.inp')
READ(1,*) N_Objects, vacuum
IF( N_Objects /= 0) THEN
  WRITE(6,*)
  WRITE(6,*)'Number of Object=',N_Objects,' vacuum=',vacuum
  ALLOCATE(objects(N_Objects))
  DO i = 1, N_Objects
    objects(i)%Shape		=0
    objects(i)%Axis			=0
    objects(i)%Dimensions	=0
    objects(i)%Locations	=0
    objects(i)%Regions		=0
    objects(i)%Phi			=0
    objects(i)%Eps			=0
    READ(1,*) objects(i)%Shape
    READ(1,*) objects(i)%Axis
    READ(1,*) objects(i)%Dimensions(:)
    READ(1,*) objects(i)%Locations(1,:)
    READ(1,*) objects(i)%Locations(2,:)
    !=========LY modification, 2022-6-13========
    READ(1,*) objects(i)%Locations(3,:)
    READ(1,*) objects(i)%Locations(4,:)
    !=========LY modification, 2022-6-13========
    READ(1,*) objects(i)%Regions
    READ(1,*) objects(i)%Direction
    READ(1,*) objects(i)%Phi
    READ(1,*) objects(i)%Eps
    READ(1,*) objects(i)%Erosion
    IF (objects(i)%Erosion > 0) THEN
      READ(1,*) objects(i)%Wall(1)%Shape
      READ(1,*) objects(i)%Wall(1)%Channelwall
      READ(1,*) objects(i)%Wall(1)%Limits(1,1:2)
      READ(1,*) objects(i)%Wall(1)%Limits(2,1:2)
      READ(1,*) objects(i)%Wall(1)%stepx
    ELSE
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
      READ(1,*)
    ENDIF
    WRITE(6,*)
    WRITE(6,*) 'objects(i)%Shape	        =',	objects(i)%Shape
    WRITE(6,*) 'objects(i)%Axis				=', objects(i)%Axis
    WRITE(6,*) 'objects(i)%Dimensions(:)	=', objects(i)%Dimensions(:)
    WRITE(6,*) 'objects(i)%Locations(1,:)	=', objects(i)%Locations(1,:)
    WRITE(6,*) 'objects(i)%Locations(2,:)   =',	objects(i)%Locations(2,:)
    !=========LY modification, 2022-6-13=========
    WRITE(6,*) 'objects(i)%Locations(3,:)	=', objects(i)%Locations(3,:)
    WRITE(6,*) 'objects(i)%Locations(4,:)   =',	objects(i)%Locations(4,:)
    !=========LY modification, 2022-6-13=========
    WRITE(6,*) 'objects(i)%Regions	        =',	objects(i)%Regions
    WRITE(6,*) 'objects(i)%Direction	        =',	objects(i)%Direction
    WRITE(6,*) 'objects(i)%Phi				=',	objects(i)%Phi
    WRITE(6,*) 'objects(i)%Eps				=',	objects(i)%Eps
    WRITE(6,*) 'objects(i)%Erosion			=',	objects(i)%Erosion
    WRITE(6,*)
  END DO
ELSE
  WRITE(6,*)
  WRITE(6,*)'No Object'
  WRITE(6,*)
  !DO i = 1, 30
  ! READ(1,*)
  !ENDDO
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
  !READ(1,*)
ENDIF



READ(1,*) N_Boundary
WRITE(6,*) 'N_Boundary	        =',	N_Boundary
num = 0
IF( N_Boundary /= 0) THEN
   ALLOCATE(bc_index(N_Boundary), bc_value(N_Boundary), bc_point_1(2,N_Boundary), bc_point_2(2,N_Boundary))
   bc_index	  = -10
   bc_value	  = 0.0
   bc_point_1(1,:) = xmin - 10. 
   bc_point_1(2,:) = ymin - 10.
   bc_point_2(1,:) = xmax + 10. 
   bc_point_2(2,:) = ymax + 10.

   DO i = 1, N_Boundary
	 READ(1, *)bc_index(i), bc_value(i), bc_point_1(1,i), bc_point_1(2,i),bc_point_2(1,i), bc_point_2(2,i)
!WRITE(6,*) bc_index(i), bc_value(i), bc_point_1(1,i), bc_point_1(2,i),bc_point_2(1,i), bc_point_2(2,i)

!================================================================================================================
     IF(bc_index(i)==-1)	THEN
	    num=num+1
	 ENDIF
!================================================================================================================

   
   ENDDO
ENDIF
WRITE(*,*) num
!===================================Add for periodic boundary conditions============================================ 
IF(num>0)THEN
    ALLOCATE(ALQ(2,N_Boundary))
	DO i=1,N_Boundary
    	READ(1,*) ALQ(1,i),ALQ(2,i)
	ENDDO
ENDIF
!===================================================================================================================

!=================================periodic boundary conditions===================
IF(num>0)THEN
    OPEN(2, ACTION='WRITE', FILE='ALQ.msh')
	DO i=1,N_Boundary
	    WRITE(2,*) ALQ(1,i),ALQ(2,i)	
	ENDDO
	CLOSE(2)
    DEALLOCATE(ALQ)
ENDIF
!=================================================================================


CLOSE(1)        !$ close object.inp


!=====================================LY modification, 2022-6-13, Checking Objects====================================
Do i = 1, N_Objects
  If (objects(i)%Shape == 5) Then   !There denote that the object is ellipse.
    xDist = objects(i)%Dimensions(1)    !radius of x axis of the ellipse
    yDist = objects(i)%Dimensions(2)    !radius of y axis of the ellipse
    Focusx1 = objects(i)%Locations(1,1) !Fx1
    Focusy1 = objects(i)%Locations(1,2) !Fy1
    Focusx2 = objects(i)%Locations(2,1) !Fx2
    Focusy2 = objects(i)%Locations(2,2) !Fy2
    FocusDist = DSQRT((Focusx2-Focusx1)**2 + (Focusy2-Focusy1)**2)  !2C: F1F2 of Distance
    
    If ((xDist**2)-(yDist**2)-(FocusDist/2)**2 < 1.0D-8) Then !There we need set a reasonable limited number.
      Write(*,*) 'Ellipse parameter is True!'
    Else
      Write(*,*) 'Ellipse parameter is not meet a**2-b**2-c**2=0'
      Write(*,*) 'a**2-b**2-c**2 = ', (xDist**2)-(yDist**2)-(FocusDist/2)**2
      Write(*,*) 'Please check Ellipse input parameter!'
      Stop
    End If
  End If
End Do
!=====================================LY modification, 2022-6-13, Checking Objects====================================


!!!!!!!!!!!!!!!!!!!!!!Start_2D!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(6,*)
WRITE(6,*) 'Start: Read in PIC input file [pic.inp]'
OPEN(5, ACTION = 'READ', FILE = './INPUT/pic.inp')

READ(5,*) 
READ(5,*) irestart, ilap
!!! ************************ bjw add for impic 2019-6-3 **********************************************  
READ(5,*) IMPIC_index
READ(5,*) Bfiled_index
READ(5,*) Density_ref
READ(5,*) Temperature_ref
READ(5,*) ParticlePerGrid
READ(5,*) affp_bjw(1:3)           !!! 每颗模拟粒子代表的真实粒子数 
!!! ************************ bjw add for impic 2019-6-3 **********************************************  
READ(5,*) N_part_tot
READ(5,*) OutputRemovedParticles, OutputGlobalMoment
READ(5,*) ispe_tot
READ(5,*) nt, dt
READ(5,*) Erostart
READ(5,*) n_updatefld
READ(5,*) n_pdiag
ALLOCATE(n_stride(ispe_tot))
READ(5,*) (n_stride(jj), jj=1, ispe_tot)
READ(5,*) n_engydiag
READ(5,*) n_gdiag
READ(5,*) n_dump
READ(5,*) index_x, index_y, index_z
READ(5,*) xcenter, ycenter, zcenter
READ(5,*) c

ALLOCATE(qs(ispe_tot), xm(ispe_tot), qm(ispe_tot))
DO jj=1, ispe_tot
	READ(5,*) qs(jj), xm(jj)
END DO


ALLOCATE(affp(ispe_tot), affp_cell(ispe_tot), dens0(ispe_tot))

DO jj=1, ispe_tot
	affp(jj) = One
	dens0(jj)= Zero
END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO iface = 1,4
	f_periodic(iface) = .TRUE.
	f_zeroe(iface)    = .FALSE.
END DO

! particle
DO iface = 1,6
	periodic(iface) = .TRUE.
	pabsorb(iface)  = .FALSE.
	preflect(iface) = .FALSE.
	pemit(iface)    = .FALSE.
END DO

! 1 field
READ(5,*) (f_periodic(iface),  iface = 1,4)
READ(5,*) (   f_zeroe(iface),  iface = 1,4)
! 2 particle
READ(5,*) (periodic(iface),  iface = 1,6)
READ(5,*) ( pabsorb(iface),  iface = 1,6)
READ(5,*) (preflect(iface),  iface = 1,6)
READ(5,*) (   pemit(iface),  iface = 1,6)

READ(5,*) (phiouter(iface),iface = 1,4)

READ(5,*) den0_ref, Te_ref, phi0_ref


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SetupPartInject_2D!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(.NOT. ALLOCATED(ipf))		ALLOCATE(ipf(ispe_tot))
IF(.NOT. ALLOCATED(sdis))		ALLOCATE(sdis(ispe_tot))
IF(.NOT. ALLOCATED(N_inject))	ALLOCATE(N_inject(ispe_tot))
IF(.NOT. ALLOCATED(N_inject_x))	ALLOCATE(N_inject_x(ispe_tot))
IF(.NOT. ALLOCATED(N_inject_y))	ALLOCATE(N_inject_y(ispe_tot))
IF(.NOT. ALLOCATED(Length))		ALLOCATE(Length(ispe_tot))
IF(.NOT. ALLOCATED(R_start))	ALLOCATE(R_start(ispe_tot))
IF(.NOT. ALLOCATED(Z_start))	ALLOCATE(Z_start(ispe_tot))
IF(.NOT. ALLOCATED(vd))			ALLOCATE(vd(ispe_tot,3))
IF(.NOT. ALLOCATED(vt))			ALLOCATE(vt(ispe_tot,3))
IF(.NOT. ALLOCATED(tmpj))		ALLOCATE(tmpj(ispe_tot))
IF(.NOT. ALLOCATED(theta0))     ALLOCATE(theta0(ispe_tot))

READ(5,*) ispe_inject
DO jj = 1, ispe_inject
	READ(5,*) ipf(jj)
    READ(5,*) dens0(ipf(jj))  !!! 真实密度
    READ(5,*) theta0(ipf(jj))  !!! 入射角度
    READ(5,*) tmpj(ipf(jj))   !!! 真实温度(eV)
    !READ(5,*) N_inject_x(ipf(jj)),N_inject_y(ipf(jj))
    READ(5,*) vd(ipf(jj),1), vd(ipf(jj),2), vd(ipf(jj),3)	
    READ(5,*) vt(ipf(jj),1), vt(ipf(jj),2), vt(ipf(jj),3)
	READ(5,*) Length(ipf(jj))		
	READ(5,*) Z_start(ipf(jj))		
	READ(5,*) R_start(ipf(jj))		
END DO

READ(5,*)  Te_Sec

CLOSE(5) !$ ab.ZWZ 2021/9/11



!!! ************************ bjw add for impic 2019-6-3 **********************************************  
!!!!!!!!!!!!!!!!!!!!!!!!!normalize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 质量参考
m_ref = Me
 
!!! 电荷量参考
q_ref = E

!!! 介电常数参考
epsilon_ref = EPSILON0

!!! 粒子数密度参考
n_ref = Density_ref


!!! 温度参考
T_ref = Temperature_ref * eVtoK   !(1eV=11605K)

!!! 长度参考(德拜长度)
L_ref = SQRT(EPSILON0 * kB * T_ref / n_ref / E / E)


!$ ====== For mcc benchmark ===== \\
!L_ref = 0.01/32             !$ a cell's length equals to Debye length (implicit)
!L_ref = 0.02/128             !$ (explicit)
!L_ref = 2*0.01/32             !$ 2 cell's length equals to Debye length
!L_ref = 0.067/128             !> turner2013 1st
!n_ref = EPSILON0 * kB * T_ref / L_ref / L_ref / E / E
!$ ====== For mcc benchmark ===== //

!!! 时间参考(电子回旋频率的倒数) --- zwz: 电子振荡频率吧
time_ref = SQRT(EPSILON0 * Me / n_ref / E / E)
WRITE(6,*) ' ***** time_ref=',time_ref

!!! 速度参考
!v_ref = L_ref / time_ref
v_ref = SQRT(kB * T_ref / Me)

!!! 电势参考
Phi_ref = kB * T_ref / q_ref

!!! 电势参考
Efield_ref = Phi_ref / L_ref

!!! 磁感应强度参考
Bfield_ref = m_ref / q_ref / time_ref


!! ====== For mcc benchmark ===== \\
!L_ref = 0.01/32             !$ a cell's length equals to Debye length
!L_ref = 2*0.01/32             !$ 2 cell's length equals to Debye length
!n_ref = EPSILON0 * kB * T_ref / L_ref / L_ref / E / E
print*, 'n_ref=',n_ref
!dt = 1/(1600*13.56D6)/time_ref
!dt = 5D-11/time_ref
!dt = 1D-10/time_ref
!dt = 2D-10/time_ref             !> (implicit)
!dt = 1.25D-11/time_ref         !> (explicit)
!Frequency = 15.0D6
Frequency = 13.56D6
!dt = (1/(400*Frequency))/time_ref

!Frequency = 27.12D6
!Frequency = 40.68D6
!Frequency = 60.0D6
Period = (1/Frequency) / (dt*time_ref)
PStep = INT(Period)
Period_res = Period - Pstep

!nt = INT(5120/13.56D6/time_ref/dt)
!nt = INT(1500/Frequency/time_ref/dt)
!nt = INT(4000/Frequency/time_ref/dt)
!nt = INT(1000/Frequency/time_ref/dt)
!nt = INT(1280/Frequency/time_ref/dt)
!! ====== For mcc benchmark ===== //


OPEN(100, ACTION = 'WRITE', FILE = './OUTPUT/physics_parameter.inp')
WRITE(100,*) '********** Real physics values **********'
WRITE(100,*) '###Domain:'
WRITE(100,*) '[Xmin,Xmax]=[',xmin*L_ref,',',xmax*L_ref,']'
WRITE(100,*) '[Ymin,Ymax]=[',ymin*L_ref,',',ymax*L_ref,']'
WRITE(100,*) '[Zmin,Zmax]=[',zmin*L_ref,',',zmax*L_ref,']'
WRITE(100,*) 

DO jj=1, ispe_tot
	WRITE(100,*) '###Partical index:', jj
    WRITE(100,*) 'mass:', xm(jj)
    WRITE(100,*) 'charge quantity:', qs(jj)
    WRITE(100,*) 'density:', dens0(ipf(jj))
    WRITE(100,*) 'temperature(eV):', tmpj(ipf(jj))   
    WRITE(100,*) 'drift velocity:', vd(ipf(jj),1), vd(ipf(jj),2), vd(ipf(jj),3)	
    WRITE(100,*) 'thermal velocity', vt(ipf(jj),1), vt(ipf(jj),2), vt(ipf(jj),3)
    WRITE(100,*) 
END DO


WRITE(100,*) '********** Dimensionless values **********'
WRITE(100,*) '###Reference values:'
WRITE(100,*) 'mass ref', m_ref
WRITE(100,*) 'charge quantity ref', q_ref 
WRITE(100,*) 'epsilon ref', epsilon_ref
WRITE(100,*) 'density ref', n_ref
WRITE(100,*) 'temperature ref', T_ref
WRITE(100,*) 'length ref', L_ref
WRITE(100,*) 'time ref', time_ref
WRITE(100,*) 'velocity ref', v_ref
WRITE(100,*) 'potential ref', Phi_ref
WRITE(100,*) 'Electric field ref', Efield_ref
WRITE(100,*) 'magnetic field ref', Bfield_ref
WRITE(100,*) 

WRITE(100,*) '###Domain:'
WRITE(100,*) '[Xmin,Xmax]=[',xmin,',',xmax,']'
WRITE(100,*) '[Ymin,Ymax]=[',ymin,',',ymax,']'
WRITE(100,*) '[Zmin,Zmax]=[',zmin,',',zmax,']'
WRITE(100,*) 





ALLOCATE(t_parameter(ispe_tot))
t_parameter = 0.5 * xm * v_ref**2 / E


DO jj=1, ispe_tot
	qs(jj)  =   qs(jj) / q_ref
	xm(jj)  =   xm(jj) / m_ref
	qm(jj)  =   qs(jj) / xm(jj)
END DO

DO i=1,4
    phiouter(i) = phiouter(i) / Phi_Ref
ENDDO

DO jj = 1, ispe_inject
    dens0(ipf(jj)) = dens0(ipf(jj)) / n_ref
    vd(ipf(jj),:)	= vd(ipf(jj),:) / v_ref			
	tmpj(ipf(jj)) = tmpj(ipf(jj)) * eVtoK / T_Ref
    vt(ipf(jj),:) = SQRT(tmpj(ipf(jj))/xm(jj))
END DO

WRITE(*,*) 'T_Ref', T_Ref
WRITE(*,*) 'L_ref', L_ref
WRITE(*,*) 'time_ref', time_ref
WRITE(*,*) 'v_ref', v_ref		
WRITE(*,*) 'Phi_Ref', Phi_Ref
WRITE(*,*) 'Efield_ref', Efield_ref
WRITE(*,*) 'Bfield_ref', Bfield_ref

WRITE(*,*) 'dens0', dens0
WRITE(*,*) 'qs', qs
WRITE(*,*) 'xm', xm
WRITE(*,*) 'qm', qm
WRITE(*,*) 'tmpj', tmpj
WRITE(*,*) 'vd', vd(ipf(1),:), vd(ipf(2),:)
WRITE(*,*) 'vt', vt(ipf(1),:), vt(ipf(2),:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! ************************ bjw add for impic 2019-6-3 **********************************************  

DO jj=1, ispe_tot
	WRITE(100,*) '###Partical index:', jj
    WRITE(100,*) 'mass:', xm(jj)
    WRITE(100,*) 'charge quantity:', qs(jj)
    WRITE(100,*) 'charge mass ratio:', qm(jj)
    WRITE(100,*) 'density:', dens0(ipf(jj))
    WRITE(100,*) 'temperature(K):', tmpj(ipf(jj))   
    WRITE(100,*) 'drift velocity:', vd(ipf(jj),1), vd(ipf(jj),2), vd(ipf(jj),3)	
    WRITE(100,*) 'thermal velocity', vt(ipf(jj),1), vt(ipf(jj),2), vt(ipf(jj),3)
    WRITE(100,*) 
END DO
CLOSE(100)


!WRITE(6,*)
!WRITE(6,*)'READ "boundary.inp"'
!
!OPEN(1, ACTION='READ', FILE='boundary.inp')
!READ(1,*) N_PIC_Boundary
!IF( N_PIC_Boundary /= 0) THEN
!    WRITE(6,*)
!    WRITE(6,*)'Number of Boundary=',N_PIC_Boundary
!    ALLOCATE(boundaries(N_PIC_Boundary))
!    DO i = 1, N_PIC_Boundary
!	    boundaries(i)%Shape		=0
!	    boundaries(i)%Axis		=0
!	    boundaries(i)%Length   	=0
!	    boundaries(i)%Locations	=0
!	    boundaries(i)%BoundType =0
!	    READ(1,*)
!	    READ(1,*) boundaries(i)%Shape
!        READ(1,*) boundaries(i)%Axis
!!        READ(1,*) boundaries(i)%Length
!        READ(1,*) boundaries(i)%Locations(1,:)
!        READ(1,*) boundaries(i)%Locations(2,:)
!        READ(1,*) boundaries(i)%BoundType
!        READ(1,*) boundaries(i)%Obound
!        READ(1,*) boundaries(i)%Erosion
!!		IF (boundaries(i)%Erosion > 0) THEN
!!			READ(1,*) boundaries(i)%Wall(1)%Shape
!!			READ(1,*) boundaries(i)%Wall(1)%Channelwall
!!		ENDIF 
!   	    WRITE(6,*) boundaries(i)%Shape
!        WRITE(6,*) boundaries(i)%Axis
!        WRITE(6,*) boundaries(i)%Length
!        WRITE(6,*) boundaries(i)%Locations(1,:)
!        WRITE(6,*) boundaries(i)%Locations(2,:)
!        WRITE(6,*) boundaries(i)%BoundType  
!        WRITE(6,*) boundaries(i)%Erosion
!    END DO
!ENDIF
!
!CLOSE(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





MaxObjects = N_Objects
vacuumRegion = vacuum
IF(.NOT.ALLOCATED(phi_object))			ALLOCATE(phi_object(MaxObjects),eps_object(MaxObjects))
IF(.NOT.ALLOCATED(shape_object))		ALLOCATE(shape_object(MaxObjects),axis_object(MaxObjects))
IF(.NOT.ALLOCATED(radius2_object))		ALLOCATE(radius2_object(MaxObjects))
IF(.NOT.ALLOCATED(f_object_l))			ALLOCATE(f_object_l(MaxObjects,2),f_object_r(MaxObjects,2))
IF(.NOT.ALLOCATED(i_object_l))			ALLOCATE(i_object_l(MaxObjects,2),i_object_r(MaxObjects,2))
IF(.NOT.ALLOCATED(region_object))		ALLOCATE(region_object(MaxObjects))
!!!! **************** bjw add 2019.9.30 *****************
IF(.NOT.ALLOCATED(direction_object))    ALLOCATE(direction_object(MaxObjects))
!!!! **************** bjw add 2019.9.30 *****************
IF(.NOT.ALLOCATED(dimension_object))	ALLOCATE(dimension_object(MaxObjects,2))

DO iobject = 1, MaxObjects
	shape_object(iobject) = objects(iobject)%shape		! Object shape
	axis_object(iobject)  = objects(iobject)%axis
	
	dimension_object(iobject,:) = objects(iobject)%Dimensions(:)
	IF(shape_object(iobject)>1) THEN
		radius2_object(iobject) = dimension_object(iobject,2)**2
	ELSE
		radius2_object(iobject) = dimension_object(iobject,1)**2
	END IF		

	f_object_l(iobject,:) = objects(iobject)%Locations(1,:)   ! Object left location

	f_object_r(iobject,:) = objects(iobject)%Locations(2,:)   ! Object left location

    !!!! **************** bjw add 2019.9.30 *****************
    !region_object(iobject) = objects(iobject)%Regions(1)      ! Object regions
    IF (objects(iobject)%Direction == -1) THEN    !!! bjw add 2019.9.30  inner
        direction_object(iobject) = objects(iobject)%Direction
	    region_object(iobject) = objects(iobject)%Regions(1)      ! Object regions
    ELSEIF (objects(iobject)%Direction == -2) THEN  !!! bjw add 2019.9.30  outer
        direction_object(iobject) = objects(iobject)%Direction
        region_object(iobject) = objects(iobject)%Regions(2)      ! Object regions
    ENDIF
    !!!! **************** bjw add 2019.9.30 *****************
    
	Idummy = objects(iobject)%Regions(2)

	phi_object(iobject) = objects(iobject)%Phi                ! Grid potential

	eps_object(iobject) = objects(iobject)%Eps                ! Grid electric permittivity
END DO



OPEN(10, ACTION = 'WRITE', FILE = './OUTPUT/normalize.inp')
WRITE(10,*) xmin, ymin
WRITE(10,*) nnx, nny
CLOSE(10)


END
