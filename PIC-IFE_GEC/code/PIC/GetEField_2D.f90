SUBROUTINE GetEfield_2D

! Updated:	11/28/2004 11:30 PM
! Purpose:	Get electric field
! this one is from 3D-Plume
! difference between this one and the one in optics:
! this one calls EfBC.f for boundary , to handel fixed Phi cases
! when domain is defined at i=0, and nx+1
! the one in optics is hard wired for E=0 at symmetric surface in x and y
! and fixed phi at z at k=0, nz+1
! get_efield.f 
! get efield then call get EfBC
! obtain the vector Efield from potential Phi
! calculate the E-field from phi E=-gradient Phi
! for inner grid points 
! use center difference scheme to get 2nd order accuracy
! eq(3-28) on p.43 of Anderson, Tannehill, Pletcher CFD text

USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D

IMPLICIT NONE

INTEGER		i, j
REAL(8)        :: pp1,pp2,pp3,pp4,pp5,pp6,qq1,qq2,qq3,qq4,qq5,qq6,x,y


WRITE(6,*) 
WRITE(6,*) 'GetEfield'

!pp1 = -0.08685
!pp2 = 0.9716
!pp3 = -4.186
!pp4 = 8.78
!pp5 = -9.331
!pp6 = 4.456
!qq1 = -10.02
!qq2 = 37.85
!qq3 = -63.84
!qq4 = 40.61

!efx(i,:) = (pp1*x**5 + pp2*x**4 + pp3*x**3 + pp4*x**2 + pp5*x + pp6) / &
!            (x**4 + qq1*x**3 + qq2*x**2 + qq3*x+ qq4)&
!            *100/Efield_ref 
        
    
	DO j=1, ny
		DO i=1, nx
				efx(i,j)= (phi(i-1,j) - phi(i+1,j))/Two/hx(1)
				efy(i,j)= (phi(i,j-1) - phi(i,j+1))/Two/hx(2)
            
                
                !y = (j-1)*hx(2) * 7.4341E-5 *100   !! cm L_ref
                !efy(i,j)=  efy(i,j) + (pp1*y**5 + pp2*y**4 + pp3*y**3 + pp4*y**2 + pp5*y + pp6) / &
                !                      (y**4 + qq1*y**3 + qq2*y**2 + qq3*y+ qq4)&
                !                       *100/134523.318       !Efield_ref
                
		END DO
	END DO

! note 2003-0325: the above are sufficient if  
! the domain is such that
!  X--X--X--.....-X--X
!     |           |
!  0  1          nx  nx+1
! in this case, i only need Efield at i=1 to nx
! 
! call EfBC is needed only when i need to use Efield at i=0, and nx+1  
! if the domain is set such that
!  X--X--X--.....-X--X
!    |              |
!  0  1          nx  nx+1
! or
!  X--X--X--.....-X--X
!  |                 |
!  0  1          nx  nx+1

! for boundary grids at i=0 and nx+1
!  use eq(3-29) and (3-30) on p.44 of Anderson, Tannehill, Pletcher CFD 
! to get 2nd order accuracy

CALL EfBC_2D

!efx = efx + 0.1487

END





