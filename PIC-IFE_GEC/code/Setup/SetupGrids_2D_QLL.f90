SUBROUTINE SetupGrids_2D_QLL(delta, xmin, xmax, ymin, ymax, zmin, zmax, nnx, nny, nnz)

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: SetupGrids_2D.f90                                  C
!
!  Purpose: Setup the PIC mesh
!                                                                      C
!  Reviewer: Yuchuan Chu                              Date: 03-May-12  C
!  Comments: modified for normalization                                C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

! 2000-0308:
! Remember that the field solver iterates from j=1,nr, and k=1,nz
! so the grid lay out is different for fixed Phi and open boundaries
! 2 types of grid lay out for the 2 boundary conditions
! I: if fix potential 
! then, phi(0) and phi(nz+1) are fixed 
! II: if open boundary
! then, phi(1) and phi(nr) are determined by the E=0 condition
! here
! for symmetric surface: i need to have the particle boundary surface
! and the vertex concide
! so here everything starts at (r,z)=(0,0).

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
USE Cell_Volume_Data !$ ab.ZWZ 2021/7/11
USE IMPIC_Data_2D, ONLY: Bfiled_index
USE IFE_Data
!=========LY modification, 2022-7-25=========
Use IFE_MAIN_PARAM, Only: SmallValue
!=========LY modification, 2022-7-25=========

IMPLICIT NONE


INTEGER		::  delta 
REAL(8)   :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER   ::  nnx, nny, nnz

INTEGER		i, j
REAL(8)		R, delta_R
REAL(8)		dummy(2)

!=========LY modification, 2022-7-25=========
Integer :: MeshNodeNum, num
Real(8) :: x, y
Real(8) :: hx_partition, hy_partition
Integer :: k
Real(8), Parameter :: Pi_LY = 3.14159265358979D0
!=========LY modification, 2022-7-25=========

WRITE(6,*)
WRITE(6,*) 'SetupGrids: Read in PIC mesh [mesh.inp]'

OPEN(1, ACTION = 'READ', FILE = './INPUT/mesh.inp')
READ(1,*) (Vert_o(i), i = 1,3)
READ(1,*)  xmax, ymax, zmax
! nx, ny, nz							: Number of mesh points in x, y, z
READ(1,*) nnx, nny, nnz

! hx(1:3)								: Grid resolution in each direction
READ(1,*) (hx(i), i = 1, 2)
READ(1,*) hz
CLOSE(1)
! xmin, ymin, xmax, ymax				: Domain limits
Vert_o(1) = xmin
Vert_o(2) = ymin
Vert_o(3) = zmin

nx = nnx
ny = nny
nz = nnz


hxi = (1._8)/hx
hzi = (1._8)/hz
! Origin for the PIC domain
Vert_o(1:2) = Vert_o(1:2) - hx
Vert_o(3) = Vert_o(3) - hz

WRITE(6,*) 'Origin, Vert_o = ', (Vert_o(i), i = 1,3)

WRITE(6,*) 'Grid spacing hx(123) = ', (hx(i),i=1,2),hz
WRITE(6,*) 'Inverse spacing = ', (hxi(i),i=1,2), hzi

! Calculate  domain length
an_gridt(1) = (nx-1)*hx(1)
an_gridt(2) = (ny-1)*hx(2)
an_gridt(3) = (nz-1)*hz
WRITE(6,*) 'nx,ny,nz: ', nx,ny,nz
WRITE(6,*) 'hx: ', hx, hz
WRITE(6,*) 'an_gridt: ', an_gridt 

! Domain boundary
f_left_wall(1:2) = Vert_o(1:2) + hx
f_left_wall(3) = Vert_o(3) + hz
! Set right side domain boundary
f_right_wall = f_left_wall + an_gridt

!Note:
!f_left_wall(1)--xmin--left boundary
!f_left_wall(2)--ymin--bottom boundary
!f_right_wall(1)--xmax--right boundary
!f_right_wall(2)--ymax--top boundary

WRITE(6,*) 'Left Walls  : ', (f_left_wall(i), i = 1,3) 
WRITE(6,*) 'Right Walls : ', (f_right_wall(i), i = 1,3) 

! Initialize the grids 
IF (.NOT.ALLOCATED(VertX)) ALLOCATE(VertX(2,0:nx+1,0:ny+1))
	DO j = 0, ny+1
		DO i = 0, nx+1
			VertX(1,i,j) = Vert_o(1) + i*hx(1)
			VertX(2,i,j) = Vert_o(2) + j*hx(2)
		END DO
	END DO

IF (.NOT.ALLOCATED(i_grid_flag))	ALLOCATE(i_grid_flag(0:nx+1,0:ny+1))	
! Initialize the grid types 
	DO j = 0, ny+1
		DO i = 0, nx+1
			i_grid_flag(i,j) = 0
		END DO
  END DO

! Initialize the field arrays
!IF (.NOT.ALLOCATED(efx))	ALLOCATE(efx(0:nx+1,0:ny+1),efy(0:nx+1,0:ny+1),efz(0:nx+1,0:ny+1))
!IF (Bfiled_index) Then
!    IF (.NOT.ALLOCATED(bfx))	ALLOCATE(bfx(0:nx+1,0:ny+1),bfy(0:nx+1,0:ny+1),bfz(0:nx+1,0:ny+1))
!    IF (.NOT.ALLOCATED(bt))     ALLOCATE(bt(0:nx+1,0:ny+1)) !$ ab.ZWZ
!End If
!IF (.NOT.ALLOCATED(phi))	ALLOCATE(phi(0:nx+1,0:ny+1),rho(0:nx+1,0:ny+1))
!IF (.NOT.ALLOCATED(rho_s))	ALLOCATE(rho_s(0:nx+1,0:ny+1,ispe_tot))
!IF (.NOT.ALLOCATED(rho_n))	ALLOCATE(rho_n(0:nx+1,0:ny+1,ispe_tot))                 !$ ab.ZWZ
!IF (.NOT.ALLOCATED(rho_p))	ALLOCATE(rho_p(0:nx+1,0:ny+1,ispe_tot))                 !$ ab.ZWZ
!IF (.NOT.ALLOCATED(Ek_s))	ALLOCATE(Ek_s(0:nx+1,0:ny+1,ispe_tot))                  !$ ab.ZWZ
!Ek_s = 0.
!IF (.NOT.ALLOCATED(Ek_tot))	ALLOCATE(Ek_tot(0:nx+1,0:ny+1,ispe_tot))                !$ ab.ZWZ
!Ek_tot = 0.
!IF (.NOT.ALLOCATED(Aver_phi))	ALLOCATE(Aver_phi(0:nx+1,0:ny+1),Aver_rho(0:nx+1,0:ny+1))
!IF (.NOT.ALLOCATED(Aver_rho_s))	ALLOCATE(Aver_rho_s(0:nx+1,0:ny+1,ispe_tot))

!IF (.NOT.ALLOCATED(Aver_Ek_s))	ALLOCATE(Aver_Ek_s(0:nx+1,0:ny+1,ispe_tot))         !$ ab.ZWZ
!Aver_Ek_s = 0.
!IF (.NOT.ALLOCATED(Aver_Ek_tot))	ALLOCATE(Aver_Ek_tot(0:nx+1,0:ny+1,ispe_tot))   !$ ab.ZWZ
!Aver_Ek_tot = 0.

!=========LY modification for Multi-Layer-Grid, 2022-7-25=========
! wsy revise rho 2022 4 18
IF (.NOT.ALLOCATED(phi))	        ALLOCATE(phi(Field_Size, 1))
IF (.NOT.ALLOCATED(rho))	        ALLOCATE(rho(Field_Size, 1))
IF (.NOT.ALLOCATED(rho_s))	      ALLOCATE(rho_s(Field_Size, 1, ispe_tot))
IF (.NOT.ALLOCATED(rho_n))	      ALLOCATE(rho_n(Field_Size, 1, ispe_tot))
IF (.NOT.ALLOCATED(rho_p))	      ALLOCATE(rho_p(Field_Size, 1, ispe_tot))
IF (.NOT.ALLOCATED(Ek_s))	        ALLOCATE(Ek_s(Field_Size, 1, ispe_tot))
IF (.NOT.ALLOCATED(Ek_tot))       ALLOCATE(Ek_tot(Field_Size, 1, ispe_tot))

IF (.NOT.ALLOCATED(Aver_phi))	    ALLOCATE(Aver_phi(Field_Size, 1))
IF (.NOT.ALLOCATED(Aver_rho))     ALLOCATE(Aver_rho(Field_Size, 1))
IF (.NOT.ALLOCATED(Aver_rho_s))	  ALLOCATE(Aver_rho_s(Field_Size, 1, ispe_tot))
IF (.NOT.ALLOCATED(Aver_Ek_s))	  ALLOCATE(Aver_Ek_s(Field_Size, 1, ispe_tot))
IF (.NOT.ALLOCATED(Aver_Ek_tot))	ALLOCATE(Aver_Ek_tot(Field_Size, 1, ispe_tot))

IF (.NOT.ALLOCATED(efx))	        ALLOCATE(efx(Field_Size,1),efy(Field_Size,1),efz(Field_Size,1))
IF (.NOT.ALLOCATED(bfx))          ALLOCATE(bfx(Field_Size,1),bfy(Field_Size,1),bfz(Field_Size,1))

Do i = 1, Field_Size
  phi(i, 1) = 0._8
  rho(i, 1) = 0._8
  rho_s(i, 1, :) = 0._8
  rho_n(i, 1, :) = 0._8
  rho_p(i, 1, :) = 0._8
  Ek_s(i, 1, :) = 0._8
  Ek_tot(i, 1, :) = 0._8
  
  Aver_phi(i, 1) = 0._8
  Aver_rho(i, 1) = 0._8
  Aver_rho_s(i, 1, :) = 0._8
  Aver_Ek_s(i, 1, :) = 0._8
  Aver_Ek_tot(i, 1, :) = 0._8
  
  efx(i, 1) = 0._8
  efy(i, 1) = 0._8
  efz(i, 1) = 0._8
  bfx(i, 1) = 0._8
  bfy(i, 1) = 0._8
  bfz(i, 1) = 0._8
End Do

!=========LY modification for Multi-Layer Mesh, 2022-7-25=========
!Here, we need to calculate the actual operating area of every mesh node for particle charge interpolation.
MeshNodeNum = Size(P_average,2)
If (.NOT.Allocated(Cell_volume_zwz))  Allocate(Cell_volume_zwz(1:2, 1:MeshNodeNum))
Cell_volume_zwz = 0.0

Do i = 1, MeshNodeNum
  
  num = 0
  Cell_volume_zwz(1,i) = i
  
  Do j = 3, 6
    If (INT(P_average(j,i)) /= 0) Then
      num = num + 1
    End If
  End Do
  
  x = P_flag(1,INT(P_average(3,i)))
  y = P_flag(2,INT(P_average(3,i)))
  
  If (num == 1) Then
  !The mesh node only correspond 1 basis node:
  !CG--domain corner node, domain edge node, interior node.
  !DG--domain corner node.
    
    If (INT(P_flag(3, INT(P_average(3,i)))) == 0) Then  !CG basis node.
      
      If ((Abs(x-xmin)<SmallValue.AND.Abs(y-ymin)<SmallValue) .OR. &
          (Abs(x-xmin)<SmallValue.AND.Abs(y-ymax)<SmallValue) .OR. &
          (Abs(x-xmax)<SmallValue.AND.Abs(y-ymin)<SmallValue) .OR. &
          (Abs(x-xmax)<SmallValue.AND.Abs(y-ymax)<SmallValue)) Then
      !CG--domain corner node.
        hx_partition = HP(1,HT(2,INT(P_flag(5,INT(P_average(3,i)))))) - HP(1,HT(1,INT(P_flag(5,INT(P_average(3,i))))))
        hy_partition = HP(2,HT(4,INT(P_flag(5,INT(P_average(3,i)))))) - HP(2,HT(1,INT(P_flag(5,INT(P_average(3,i))))))
        
        If (delta == 0) Then  !2D Cartesian coordinates.
          
          Cell_volume_zwz(2, i) = hx_partition * hy_partition / 4.0
          
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
          
          If ((Abs(x-xmin)<SmallValue.AND.Abs(y-ymin)<SmallValue) .OR. &
              (Abs(x-xmax)<SmallValue.AND.Abs(y-ymin)<SmallValue)) Then
            !bottom-left corner and bottom-right corner.
            R = ymin
            Cell_volume_zwz(2,i) = 2.0 * Pi_LY * (R+hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
          Elseif ((Abs(x-xmin)<SmallValue.AND.Abs(y-ymax)<SmallValue) .OR. &
                  (Abs(x-xmax)<SmallValue.AND.Abs(y-ymax)<SmallValue)) Then
            !top-left corner and top-right corner.
            R = ymax
            Cell_volume_zwz(2,i) = 2.0 * Pi_LY * (R-hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
          End If

        End If
        
      Elseif (Abs(x-xmin)<SmallValue .OR. Abs(x-xmax)<SmallValue .OR. &
              Abs(y-ymin)<SmallValue .OR. Abs(y-ymax)<SmallValue) Then
      !CG--domain edge node.
        hx_partition = HP(1,HT(2,Int(P_flag(5,Int(P_average(3,i)))))) - HP(1,HT(1,Int(P_flag(5,Int(P_average(3,i))))))
        hy_partition = HP(2,HT(4,Int(P_flag(5,Int(P_average(3,i)))))) - HP(2,HT(1,Int(P_flag(5,Int(P_average(3,i))))))
        
        If (delta == 0) Then  !2D Cartesian coordinates.
          
          Cell_volume_zwz(2, i) = hx_partition * hy_partition / 2.0
          
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
          
          If (Abs(y-ymin)<SmallValue) Then
            !bottom edge.
            R = ymin
            Cell_volume_zwz(2, i) = 2.0 * Pi_LY * (R+hy_partition/4.0) * (hy_partition/2.0) * hx_partition
          Elseif (Abs(y-ymax)<SmallValue) Then
            !top edge.
            R = ymax
            Cell_volume_zwz(2, i) = 2.0 * Pi_LY * (R-hy_partition/4.0) * (hy_partition/2.0) * hx_partition
          Elseif (Abs(x-xmin)<SmallValue .OR. Abs(x-xmax)<SmallValue) Then
            !left edge and right edge.
            R= y
            Cell_volume_zwz(2, i) = 2.0 * Pi_LY * R * hy_partition * hx_partition / 2.0
          End If
          
        End If
        
      Else
      !CG--interior basis node.
        hx_partition = HP(1,HT(2,Int(P_flag(5,Int(P_average(3,i)))))) - HP(1,HT(1,Int(P_flag(5,Int(P_average(3,i))))))
        hy_partition = HP(2,HT(4,Int(P_flag(5,Int(P_average(3,i)))))) - HP(2,HT(1,Int(P_flag(5,Int(P_average(3,i))))))
        
        If (delta == 0) Then  !2D Cartesian coordinates.
          
          Cell_volume_zwz(2, i) = hx_partition * hy_partition / 1.0
          
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
          
          R = y
          Cell_volume_zwz(2, i) = 2.0 * Pi_LY * R * hy_partition * hx_partition / 1.0
          
        End If
        
      End IF
      
    Elseif (Int(P_flag(3, Int(P_average(3,i)))) == 1) Then  !DG basis node.
      
      If ((Abs(x-xmin)<SmallValue.AND.Abs(y-ymin)<SmallValue) .OR. &
          (Abs(x-xmin)<SmallValue.AND.Abs(y-ymax)<SmallValue) .OR. &
          (Abs(x-xmax)<SmallValue.AND.Abs(y-ymin)<SmallValue) .OR. &
          (Abs(x-xmax)<SmallValue.AND.Abs(y-ymax)<SmallValue)) Then
      !DG--domain corner basis node.
        hx_partition = HP(1,HT(2,Int(P_flag(5,Int(P_average(3,i)))))) - HP(1,HT(1,Int(P_flag(5,Int(P_average(3,i))))))
        hy_partition = HP(2,HT(4,Int(P_flag(5,Int(P_average(3,i)))))) - HP(2,HT(1,Int(P_flag(5,Int(P_average(3,i))))))
        
        If (delta == 0) Then  !2D Cartesian coordinates.
          
          Cell_volume_zwz(2, i) = hx_partition * hy_partition / 4.0
          
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
          
          If ((Abs(x-xmin)<SmallValue.AND.Abs(y-ymin)<SmallValue) .OR. &
              (Abs(x-xmax)<SmallValue.AND.Abs(y-ymin)<SmallValue)) Then
            !bottom-left corner and bottom-right corner.
            R = ymin
            Cell_volume_zwz(2, i) = 2.0 * Pi_LY * (R+hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
          Elseif ((Abs(x-xmin)<SmallValue.AND.Abs(y-ymax)<SmallValue) .OR. &
                  (Abs(x-xmax)<SmallValue.AND.Abs(y-ymax)<SmallValue)) Then
            !top-left corner and top-right corner.
            R = ymax
            Cell_volume_zwz(2, i) = 2.0 * Pi_LY * (R-hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
          End If

        End If
        
      Else
        Write(6,*) 'error condition: 1 DG, please revise code.'
        Write(6,*) 'error locate on SetupGirds_2D_QLL.f90'
        Stop
      End If
      
    End If
    
  Elseif (num == 2) Then
  !The mesh node correspond 2 basis node:
  !1 CG and 1 DG--domain edge node, interior node.
  !2 DG--domain edge node, interior node.
    
    If ((Int(P_flag(3, Int(P_average(3,i))))==0.AND.Int(P_flag(3, Int(P_average(4,i))))==1) .OR. &
        (Int(P_flag(3, Int(P_average(3,i))))==1.AND.Int(P_flag(3, Int(P_average(4,i))))==0)) Then  !1CG, 1DG basis node.
      
      If (Abs(x-xmin)<SmallValue.OR.Abs(x-xmax)<SmallValue.OR.Abs(y-ymin)<SmallValue.OR.Abs(y-ymax)<SmallValue) Then
      !1 CG and 1 DG--domain edge node.
        
        Do k = 3, 4
          
          hx_partition = HP(1,HT(2,Int(P_flag(5,Int(P_average(k,i)))))) - HP(1,HT(1,Int(P_flag(5,Int(P_average(k,i))))))
          hy_partition = HP(2,HT(4,Int(P_flag(5,Int(P_average(k,i)))))) - HP(2,HT(1,Int(P_flag(5,Int(P_average(k,i))))))
          
          If (delta == 0) Then  !2D Cartesian coordinates.
            
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
            
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            
            If (Abs(y-ymin)<SmallValue) Then
              !bottom edge.
              R = ymin
              Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * &
                                      (R+hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
            Elseif (Abs(y-ymax)<SmallValue) Then
              !top edge.
              R = ymax
              Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * &
                                      (R-hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
            Elseif (Abs(x-xmin)<SmallValue .OR. Abs(x-xmax)<SmallValue) Then
              !left edge and right edge.
              R = y
              Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 2.0 / 2.0
            End If
            
          End If
        End Do
        
      Else 
      !1 CG and 1 DG--interior node.
        If (Int(P_flag(3, Int(P_average(3,i))))==0) Then  !P_average(3,i)--CG, P_average(4,i)--DG
          
          hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(3,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(3,i))))))
          hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(3,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(3,i))))))
          
          If (delta == 0) Then  !2D Cartesian coordinates.
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition * 3.0 / 4.0
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            R = y
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition * 3.0 / 4.0
          End If
          
          hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(4,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(4,i))))))
          hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(4,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(4,i))))))
          
          If (delta == 0) Then  !2D Cartesian coordinates.
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            R = y
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
          End If
          
        Else  !P_average(3,i)--DG, P_average(4,i)--CG
          
          hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(3,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(3,i))))))
          hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(3,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(3,i))))))
          
          If (delta == 0) Then  !2D Cartesian coordinates.
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            R = y
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
          End If
          
          hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(4,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(4,i))))))
          hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(4,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(4,i))))))
          
          If (delta == 0) Then  !2D Cartesian coordinates.
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition * 3.0 / 4.0
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            R = y
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition * 3.0 / 4.0
          End If
        End If
      End If
      
    Elseif (Int(P_flag(3, Int(P_average(3,i))))==1.AND.Int(P_flag(3, Int(P_average(4,i))))==1) Then !2DG basis node.

      If (Abs(x-xmin)<SmallValue.OR.Abs(x-xmax)<SmallValue.OR.Abs(y-ymin)<SmallValue.OR.Abs(y-ymax)<SmallValue) Then
      !2DG--domain edge node.
        Do k = 3, 4
          
          hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
          hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
          
          If (delta == 0) Then  !2D Cartesian coordinates.
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            If (Abs(y-ymin)<SmallValue) Then
              !bottom edge.
              R = ymin
              Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * &
                                      (R+hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
            Elseif (Abs(y-ymax)<SmallValue) Then
              !top edge.
              R = ymax
              Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * &
                                      (R-hy_partition/4.0) * (hy_partition/2.0) * hx_partition / 2.0
            Elseif (Abs(x-xmin)<SmallValue .OR. Abs(x-xmax)<SmallValue) Then
              !left edge and right edge.
              R = y
              Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 2.0 / 2.0
            End If
          End If
        End Do
      Else
      !2DG--interior node.
        Do k = 3, 4
          hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
          hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
          If (delta == 0) Then  !2D Cartesian coordinates.
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
          Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            R = y
            Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
          End If
        End Do
      End If
      
    Else
      Write(6,*) 'error condition: num=2, please revise code.'
      Write(6,*) 'error locate on SetupGirds_2D_QLL.f90'
      Stop
    End If
    
  Elseif (num == 3) Then
  !The mesh node correspond 3 basis node:
  !1CG and 2 DG--interior node.
    
    If (Int(P_flag(3, Int(P_average(3,i))))==0) Then  !P_average(3,i)--CG, P_average(4,i) and P_average(5,i)--DG
      hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(3,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(3,i))))))
      hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(3,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(3,i))))))
      
      If (delta == 0) Then  !2D Cartesian coordinates.
        Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 2.0
      Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
        R = y
        Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 2.0
      End If
      
      Do k = 4, 5
        hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        If (delta == 0) Then  !2D Cartesian coordinates.
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
          R = y
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
        End If
      End Do
      
    Elseif (Int(P_flag(3, Int(P_average(4,i))))==0) Then  !P_average(4,i)--CG, P_average(3,i) and P_average(5,i)--DG
      hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(4,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(4,i))))))
      hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(4,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(4,i))))))
      
      If (delta == 0) Then  !2D Cartesian coordinates.        
        Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 2.0
      Elseif (delta == 1) Then  !2D axis-symmetric coordinates.        
        R = y
        Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 2.0
      End If
      
      Do k = 3, 5, 2
        hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        If (delta == 0) Then  !2D Cartesian coordinates.          
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.          
          R = y
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
        End If
      End Do
      
    Elseif (Int(P_flag(3, Int(P_average(5,i))))==0) Then  !P_average(5,i)--CG, P_average(3,i) and P_average(4,i)--DG
      hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(5,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(5,i))))))
      hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(5,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(5,i))))))
      
      If (delta == 0) Then  !2D Cartesian coordinates.        
        Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 2.0
      Elseif (delta == 1) Then  !2D axis-symmetric coordinates.        
        R = y
        Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 2.0
      End If
      
      Do k = 3, 4
        hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        If (delta == 0) Then  !2D Cartesian coordinates.          
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.          
          R = y
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
        End If
      End Do
      
    Else
      Write(6, *) 'error condition:1 CG and 2 DG, please revise code' 
      Write(6, *) 'error locate on SetupGirds_2D_QLL.f90'
      Stop
    End If
    
  Elseif (num == 4) Then
  !The mesh node correspond 4 basis node:
  !1CG and 3 DG--interior node.
  !4DG--interior node.
    
    If (Int(P_flag(3, Int(P_average(3,i))))==1 .AND. Int(P_flag(3, Int(P_average(4,i))))==1 .AND. &
        Int(P_flag(3, Int(P_average(5,i))))==1 .AND. Int(P_flag(3, Int(P_average(6,i))))==1) Then
    !4DG--interior node.
      Do k = 3, 6
        hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        If (delta == 0) Then  !2D Cartesian coordinates.          
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.          
          R = y
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
        End If
      End Do
      
    Else
    !1 CG and 3 DG--interior node.
      Do k = 3, 6
        hx_partition = HP(1, HT(2, Int(P_flag(5, Int(P_average(k,i)))))) - HP(1, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        hy_partition = HP(2, HT(4, Int(P_flag(5, Int(P_average(k,i)))))) - HP(2, HT(1, Int(P_flag(5, Int(P_average(k,i))))))
        If (delta == 0) Then  !2D Cartesian coordinates.          
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + hx_partition * hy_partition / 4.0
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.          
          R = y
          Cell_volume_zwz(2, i) = Cell_volume_zwz(2, i) + 2.0 * Pi_LY * R * hy_partition * hx_partition / 4.0
        End If
      End Do
    End If
    
  Else
    Write(6,*) 'error num condition, please revise code.', num
    Write(6,*) 'error locate on "SetupGrids_2D_QLL.f90"'
    Stop
  End If
  
End Do

!=========Checking Actual action area=========
Open(1, ACTION = 'WRITE', FILE = './OUTPUT/CellVolume.dat')
  Write(1,*) 'TITLE = "CellVolume Plot"'
  Write(1,"(A150)") 'VARIABLES = "x" "y" "Cell_volume"'
  Do i = 1, MeshNodeNum
    Write(1,"(E15.6,' ',E15.6,' ',E15.6)") P_average(1,i), P_average(2,i), Cell_volume_zwz(2,i)
  End Do
Close(1)
!=========LY modification for Multi-Layer Mesh, 2022-7-25=========

!=========Old Code=========
!$ =========== ab.ZWZ to setup Cell_Volume =========== \\ 
!-------------  unuse-------------
!IF (.NOT.ALLOCATED(Cell_volume_zwz))	ALLOCATE(Cell_volume_zwz(1:nx,1:ny))
!
!IF(delta == 0)THEN
!	Cell_volume_zwz(:,:) = hx(2)*hx(1) 
!    
!    DO i = 1, nx
!        Cell_volume_zwz(i,1) = Cell_volume_zwz(i,1)/2   ! bottom
!        Cell_volume_zwz(i,ny) = Cell_volume_zwz(i,ny)/2 ! top
!    ENDDO
!    DO j = 1, ny  
!        Cell_volume_zwz(1,j) = Cell_volume_zwz(1,j)/2   !left
!        Cell_volume_zwz(nx,j) = Cell_volume_zwz(nx,j)/2 !right
!    ENDDO
!    
!ELSEIF(delta == 1)THEN
!!----- except bottom and top
!	DO j = 2, ny-1
!		R = VertX(2,1,j)
!		DO i = 1, nx
!			Cell_volume_zwz(i,j) = (2.*PI*R*hx(2))*hx(1)
!		ENDDO
!    ENDDO
!    delta_R=hx(2)
!    DO i = 1, nx    !$ bottom
!        R = VertX(2,i,1)
!        Cell_volume_zwz(i,1) = (2.*PI*(R+hx(2)/4))*hx(2)/2.*hx(1)
!    ENDDO
!    DO i = 1, nx    !$ top
!        R = VertX(2,i,ny)
!        Cell_volume_zwz(i,ny) = (2.*PI*(R-hx(2)/4))*hx(2)/2.*hx(1)      
!    ENDDO
!    
!    !--- left and right
!    DO j = 1, ny 
!        Cell_volume_zwz(1,j) = Cell_volume_zwz(1,j)/2
!        Cell_volume_zwz(nx,j) = Cell_volume_zwz(nx,j)/2
!    ENDDO  
!!$ =========== ab.ZWZ to setup Cell_volume_zwz =========== //
!ENDIF
!
!OPEN(1, ACTION = 'WRITE', FILE = './OUTPUT/CellVolume.dat')
!    WRITE(1,*) 'TITLE = "Field Plot"'
!    WRITE(1,"(A150)") 'VARIABLES = "x" "y" "Cell_volume"'!$ modified by ZWZ  2020/10/28
!    WRITE(1,"(' ZONE I = ',I6,', J= ',I6)") nx, ny
!	DO j=1,ny
!		DO i=1,nx
!            WRITE(1,"(E15.6,' ',E15.6,' ',E15.6)") VertX(1:2,i,j),Cell_volume_zwz(i,j)      
!		END DO
!	END DO
!CLOSE(1)
!=========Old Code=========

END