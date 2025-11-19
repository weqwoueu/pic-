SUBROUTINE SetupBfield_circle(xmin, xmax, ymin, ymax, zmin, zmax, nnx, nny, nnz, VertX)
  
!!!! bjw add 2019.10.23

!USE PIC_MAIN_PARAM_2D
!USE Domain_2D
USE Field_2D
USE Constant_Variable_2D
!USE Particle_2D

IMPLICIT NONE

REAL(8)        :: xmin, xmax, ymin, ymax, zmin, zmax
INTEGER        :: nnx, nny, nnz
REAL(8)		   :: VertX(2,0:nnx+1,0:nny+1)


INTEGER		   :: i, j

REAL(8)        :: x0,y0,x,y,r,cosTheta,sinTheta


REAL(8)        :: p1,p2,p3,p4,p5,p6,q1,q2,q3,q4,q5,q6
REAL(8)        :: pp1,pp2,pp3,pp4,pp5,pp6,qq1,qq2,qq3,qq4,qq5,qq6


!p1 = -3.589
!p2 = 23.34
!p3 = -41.34
!p4 = 45.23
!p5 = -17.13
!p6 = 1.462
!q1 = -7.573
!q2 = 23.23
!q3 = -35.17
!q4 = 28.05
!q5 = -7.276
!!!!注：该方程为拟合曲线，x<0-3.5cm>,y<0-20.5mT>
!!fB(x) = (p1*x**5 + p2*x**4 + p3*x**3 + p4*x**2 + p5*x + p6) / &
!!        (x**5 + q1*x**4 + q2*x**3 + q3*x**2 + q4*x + q5)  
!
!!pp1 = 0.01136
!!pp2 = -0.02071
!!pp3 = -0.07823
!!pp4 = 0.003637
!!pp5 = 0.6333
!!qq1 = -9.831
!!qq2 = 36.51
!!qq3 = -60.84
!!qq4 = 38.81
!!qq5 = -0.7761
!!!!注：该方程为拟合曲线，x<0-3.5cm>,y<0-6 100V/m>
!!fE(x) = (pp1*x**4 + pp2*x**3 + pp3*x**2 + pp4*x + pp5) /
!!    (x**5 + qq1*x**4 + qq2*x**3 + qq3*x**2 + qq4*x + qq5)
!
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
!!fE(x) = (pp1*x**5 + pp2*x**4 + pp3*x**3 + pp4*x**2 + pp5*x + pp6) /
!!    (x**4 + qq1*x**3 + qq2*x**2 + qq3*x+ qq4)
!
!
!DO i = 1,nnx
!    
!    DO j = 1,nny
!        
!        x = VertX(1,i,j) * L_ref *100   !! cm
!        y = VertX(2,i,j) * L_ref *100
!        
!        bfx(i,j) = 0.0
!        bfy(i,j) = 0.0
!        !bfz(i,:) = (p1*x**5 + p2*x**4 + p3*x**3 + p4*x**2 + p5*x + p6) / &
!        !            (x**5 + q1*x**4 + q2*x**3 + q3*x**2 + q4*x + q5) &
!        !            /1000/Bfield_ref 
!        bfz(i,j) = (p1*y**5 + p2*y**4 + p3*y**3 + p4*y**2 + p5*y + p6) / &
!                    (y**5 + q1*y**4 + q2*y**3 + q3*y**2 + q4*y + q5) &
!                    /1000/Bfield_ref 
!    
!    
!        
!        !efx(i,:) = (pp1*x**4 + pp2*x**3 + pp3*x**2 + pp4*x + pp5) / &
!        !           (x**5 + qq1*x**4 + qq2*x**3 + qq3*x**2 + qq4*x + qq5) &
!        !           *100/Efield_ref 
!          
!        !efx(i,:) = (pp1*x**5 + pp2*x**4 + pp3*x**3 + pp4*x**2 + pp5*x + pp6) / &
!        !            (x**4 + qq1*x**3 + qq2*x**2 + qq3*x+ qq4)&
!        !            *100/Efield_ref
!        
!        !efy(i,:) = 0.0
!        
!        !efx(i,j) = 0.0  
!        !efy(i,j) = (pp1*y**5 + pp2*y**4 + pp3*y**3 + pp4*y**2 + pp5*y + pp6) / &
!        !            (y**4 + qq1*y**3 + qq2*y**2 + qq3*y+ qq4)&
!        !            *100/Efield_ref
!        !
!        !efz(i,j) = 0.0    
!    
!    END DO
!    
!END DO


!x0 = (xmin + xmax) / 2.0 
!y0 = (ymin + ymax) / 2.0 
!
!!!! cos Theta = (x-x0)/r , sin Theta = (y-y0)/r
!!!! if cos Theta > 0 and sin Theta > 0 , 第一象限
!!!! if cos Theta < 0 and sin Theta > 0 , 第二象限
!!!! if cos Theta < 0 and sin Theta < 0 , 第三象限
!!!! if cos Theta > 0 and sin Theta < 0 , 第四象限
!
DO i = 1,nnx
    DO j = 1, nny
        x = VertX(1,i,j)
        y = VertX(2,i,j)

        IF (y<=340.0) THEN
            bfx(i,j) = 0.0
            bfy(i,j) = 0.0
            bfz(i,j) = 170.0*EXP(-5.0*(y/340.0-1.0)*(y/340.0-1.0))/10000.0/Bfield_ref
        ELSE
            bfx(i,j) = 0.0
            bfy(i,j) = 0.0
            bfz(i,j) = 170.0*EXP(-1.5*(y/340.0-1.0)*(y/340.0-1.0))/10000.0/Bfield_ref
        END IF
   	END DO
END DO     
   


OPEN(1, ACTION = 'WRITE', FILE = 'Bfield.dat')
WRITE(1,*) 'TITLE = "Field Plot"'
WRITE(1,*) 'VARIABLES = "x" "y" "bf" "bfx" "bfy" "bfz"'
WRITE(1,40) nnx, nny
	DO j=1,nny
		DO i=1,nnx
			WRITE(1,50) VertX(1:2,i,j),SQRT(bfx(i,j)*bfx(i,j)+bfy(i,j)*bfy(i,j)+bfz(i,j)*bfz(i,j)), bfx(i,j), bfy(i,j), bfz(i,j)
		END DO
	END DO
40 FORMAT (' ZONE I = ',I6,', J= ',I6)
50 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
CLOSE(1)


OPEN(1, ACTION = 'WRITE', FILE = 'Efield.dat')
WRITE(1,*) 'TITLE = "Field Plot"'
WRITE(1,*) 'VARIABLES = "x" "y" "ef" "efx" "efy" "efz"'
WRITE(1,40) nnx, nny
	DO j=1,nny
		DO i=1,nnx
			WRITE(1,50) VertX(1:2,i,j),SQRT(efx(i,j)*efx(i,j)+efy(i,j)*efy(i,j)+efz(i,j)*efz(i,j)), efx(i,j), efy(i,j), efz(i,j)
		END DO
	END DO
!40 FORMAT (' ZONE I = ',I6,', J= ',I6)
!50 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
CLOSE(1)




END        
