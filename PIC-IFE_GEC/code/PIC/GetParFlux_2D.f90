SUBROUTINE GetParFlux_2D

! Updated:	11/9/2004 01:30 AM
! Purpose:	Get particle flux density for diagnostics.

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Particle_2D

IMPLICIT NONE

INTEGER		i, j, nden_var, i_part, isp, ii, jj, kk     !!!!!!!!!!!Ie_x是x方向的电子电流
REAL(8)		rxp, ryp, rzp, dx, dy, dz, xcellmdx, ycellmdy, zcellmdz

WRITE(6,*) 'GetParFlux'

IF (.NOT.ALLOCATED(vj_sx)) THEN
	ALLOCATE(vj_sx(0:nx+1,0:ny+1,ispe_tot), vj_sy(0:nx+1,0:ny+1,ispe_tot), vj_sz(0:nx+1,0:ny+1,ispe_tot))
	ALLOCATE(vj_totx(0:nx+1,0:ny+1), vj_toty(0:nx+1,0:ny+1), vj_totz(0:nx+1,0:ny+1))
END IF

! zerc zero out the flux density density (including the guard cells)
vj_sx = Zero
vj_sy = Zero
vj_sz = Zero
vj_totx = Zero
vj_toty = Zero
vj_totz = Zero
Ie_x = Zero
Ie_y = Zero
Ie_z = Zero
            
DO i_part = 1, ntot

! 1st find out the species
	isp = INT(part(i_part,7))

! information on particle position on the grid
	rxp=(part(i_part,1) - Vert_o(1))*hxi(1)
	i=rxp
	dx=rxp-i

	ryp=(part(i_part,2) - Vert_o(2))*hxi(2)
	j=ryp
	dy=ryp-j

!	rzp=(part(i_part,3) - Vert_o(3))*hxi(3)
!	k=rzp
!	dz=rzp-k

! cellx-dx, celly-dy, cellz-dz
	xcellmdx = 1. -dx
	ycellmdy = 1. -dy
!	zcellmdz = 1. -dz

! deposit velocity

 	vj_sx(i,j,isp)  = vj_sx(i,j,isp) + xcellmdx*ycellmdy*part(i_part,4)
 	vj_sy(i,j,isp)  = vj_sy(i,j,isp) + xcellmdx*ycellmdy*part(i_part,5)
 	vj_sz(i,j,isp)  = vj_sz(i,j,isp) + xcellmdx*ycellmdy*part(i_part,6)

	vj_sx(i+1,j,isp) = vj_sx(i+1,j,isp) + dx*ycellmdy*part(i_part,4)
 	vj_sy(i+1,j,isp) = vj_sy(i+1,j,isp) + dx*ycellmdy*part(i_part,5)
 	vj_sz(i+1,j,isp) = vj_sz(i+1,j,isp) + dx*ycellmdy*part(i_part,6)

 	vj_sx(i,j+1,isp) = vj_sx(i,j+1,isp) + xcellmdx*dy*part(i_part,4)
 	vj_sy(i,j+1,isp) = vj_sy(i,j+1,isp) + xcellmdx*dy*part(i_part,5)
 	vj_sz(i,j+1,isp) = vj_sz(i,j+1,isp) + xcellmdx*dy*part(i_part,6)

 	vj_sx(i+1,j+1,isp) = vj_sx(i+1,j+1,isp) + dx*dy*part(i_part,4)
 	vj_sy(i+1,j+1,isp) = vj_sy(i+1,j+1,isp) + dx*dy*part(i_part,5)
 	vj_sz(i+1,j+1,isp) = vj_sz(i+1,j+1,isp) + dx*dy*part(i_part,6)

 	vj_sx(i,j,isp) = vj_sx(i,j,isp) + xcellmdx*ycellmdy*part(i_part,4)
 	vj_sy(i,j,isp) = vj_sy(i,j,isp) + xcellmdx*ycellmdy*part(i_part,5)
 	vj_sz(i,j,isp) = vj_sz(i,j,isp) + xcellmdx*ycellmdy*part(i_part,6)

 	vj_sx(i,j+1,isp) = vj_sx(i,j+1,isp) + xcellmdx*dy*part(i_part,4)
 	vj_sy(i,j+1,isp) = vj_sy(i,j+1,isp) + xcellmdx*dy*part(i_part,5)
 	vj_sz(i,j+1,isp) = vj_sz(i,j+1,isp) + xcellmdx*dy*part(i_part,6)

 	vj_sx(i+1,j,isp) = vj_sx(i+1,j,isp) + dx*ycellmdy*part(i_part,4) 
 	vj_sy(i+1,j,isp) = vj_sy(i+1,j,isp) + dx*ycellmdy*part(i_part,5) 
 	vj_sz(i+1,j,isp) = vj_sz(i+1,j,isp) + dx*ycellmdy*part(i_part,6) 

 	vj_sx(i+1,j+1,isp) = vj_sx(i+1,j+1,isp) + dx*dy*part(i_part,4)
 	vj_sy(i+1,j+1,isp) = vj_sy(i+1,j+1,isp) + dx*dy*part(i_part,5)
 	vj_sz(i+1,j+1,isp) = vj_sz(i+1,j+1,isp) + dx*dy*part(i_part,6)

! end particle loop
END DO

! normalize the current density and find the total current density
DO isp=1,ispe_tot

		DO j=0, ny+1
			DO i=0, nx+1

				vj_sx(i,j,isp)=vj_sx(i,j,isp)*affp_cell(isp)*qs(isp)
				vj_sy(i,j,isp)=vj_sy(i,j,isp)*affp_cell(isp)*qs(isp)
				vj_sz(i,j,isp)=vj_sz(i,j,isp)*affp_cell(isp)*qs(isp)
!!!!!!!!!!---------------------------QLL 节点速度
!
!                IF(rho_s(i,j,k,isp)/=0)THEN
!				vx_s(i,j,k,isp)=vj_sx(i,j,k,isp)/rho_s(i,j,k,isp)
!				vy_s(i,j,k,isp)=vj_sy(i,j,k,isp)/rho_s(i,j,k,isp)
!				vz_s(i,j,k,isp)=vj_sz(i,j,k,isp)/rho_s(i,j,k,isp)
!                ENDIF
!!!!!!!!!!---------------------------QLL

!  summer over all species
!  total  density
				vj_totx(i,j) = vj_totx(i,j) + vj_sx(i,j,isp)
				vj_toty(i,j) = vj_toty(i,j) + vj_sy(i,j,isp)
				vj_totz(i,j) = vj_totz(i,j) + vj_sz(i,j,isp)
			END DO          
		END DO

END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11-------------------------------QLL

	DO j = 1, ny
		DO i = 1, nx

        Ie_x = Ie_x+vj_sx(i,j,1)
        Ie_y = Ie_y+vj_sy(i,j,1)
        Ie_z = Ie_z+vj_sz(i,j,1)
        
        ENDDO
    ENDDO



OPEN(1234,ACTION='WRITE',FILE='Vxyz.dat')
WRITE(1234,*) 'TITLE = "Field Plot"'
WRITE(1234,*) 'VARIABLES = "x" "y" "Jx" "Jy" "Jz" '
WRITE(1234,456) nx, ny


	DO jj = 1, ny
		DO ii = 1, nx

        WRITE(1234,567)VertX(1:2,ii,jj),vj_sx(ii,jj,1),vj_sy(ii,jj,1),vj_sz(ii,jj,1)
        
        ENDDO
    ENDDO

!OPEN(1234,ACTION='WRITE',FILE='Vxyz_ion.dat')
!WRITE(1234,*) 'TITLE = "Field Plot"'
!WRITE(1234,*) 'VARIABLES = "x" "y" "Jx" "Jy" "Jz" '
!WRITE(1234,456) nx, ny
!
!
!	DO jj = 1, ny
!		DO ii = 1, nx
!
!        WRITE(1234,567)VertX(1:2,ii,jj),vj_sx(ii,jj,1),vj_sy(ii,jj,1),vj_sz(ii,jj,1)
!        
!        ENDDO
!    ENDDO
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11-------------------------------QLL


456 FORMAT (' ZONE I = ',I6,', J= ',I6)
567 FORMAT (F8.4,' ',F8.4,' ',F8.4,' ',F14.7,' ',F14.7,' ',F14.7)
CLOSE(1234)
      
END