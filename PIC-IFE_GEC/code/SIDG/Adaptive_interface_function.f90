SUBROUTINE Adaptive_interface_function(x, y, objects, temp)

USE Object_Data_2D

IMPLICIT NONE

!TYPE(ObjectType), DIMENSION(:), INTENT(IN) :: objects
TYPE(ObjectType), INTENT(IN) :: objects
INTEGER             :: i, N_Objects
REAL(8), INTENT(IN) :: x, y
REAL(8)             :: temp, pi, r, x_c, y_c, x_r, y_r
PARAMETER (pi = 3.141592653589793)
!=========LY modification, 2022-7-25=========
Real(8) :: Fx1, Fy1, Fx2, Fy2
Real(8) :: slope, theta
Real(8) :: thetaA, thetaB, thetaC
!=========LY modification, 2022-7-25=========

i = 0
r = 0.0
x_c = 0.0
y_c = 0.0
x_r = 0.0
y_r = 0.0
temp = 0.0 


!N_Objects = SIZE(objects, 1)
!DO i = 1, N_Objects
  IF (objects%Shape==1 .AND. objects%Axis==0) THEN
    
    r = objects%Dimensions(1)
    x_c = objects%Locations(1,1)
    y_c = objects%Locations(1,2)
    temp = (x-x_c)**2 + (y-y_c)**2 -r**2
    
  ELSEIF (objects%Shape==5) THEN    !ellipse, LY modification, 2022-5-16
    
    !=========LY modification, 2022-5-30=========
    x_r = objects%Dimensions(1)   !radius of x axis of the ellipse
    y_r = objects%Dimensions(2)   !radius of y axis of the ellipse
    Fx1 = objects%Locations(1,1)  !focus1 x of the ellipse
    Fy1 = objects%Locations(1,2)  !focus1 y of the ellipse
    Fx2 = objects%Locations(2,1)  !focus2 x of the ellipse
    Fy2 = objects%Locations(2,2)  !focus2 y of the ellipse
    
    x_c = (Fx1+Fx2)/2.0
    y_c = (Fy1+Fy2)/2.0
    !=========LY modification, 2022-5-30=========
    
    !=========LY modification, 2022-6-1=========
    If (ABS(Fy1-Fy2) < 1.0D-8) Then    !horizontal ellipse, focus on x axis.
      temp = ((x-x_c)**2)/(x_r**2) + ((y-y_c)**2)/(y_r**2)
    Elseif (ABS(Fx1-Fx2) < 1.0D-8) Then    !vertical ellipse, focus on y axis.
      temp = ((x-x_c)**2)/(x_r**2) + ((y-y_c)**2)/(y_r**2)
    Else    !inclined ellipse
      slope = (Fy2-Fy1) / (Fx2-Fx1)
      If (ABS(slope) < 1.0D-8) Then   !There denote slope=0, but we already deal with the condition in above.
        Write(*,*) 'The slope = ', slope
        Stop
      End If
      
      If (slope > 1.0D-8) Then
        theta = DATAN(slope)
      Elseif (slope < -1.0D-8) Then
        theta = DATAN(slope) + pi
      End If
      thetaA = (((DCOS(theta))**2)/(x_r**2) + ((DSIN(theta))**2)/(y_r**2))
      thetaB = (1/(x_r**2) - 1/(y_r**2)) * DSIN(2*theta)
      thetaC = (((DSIN(theta))**2)/(x_r**2) + ((DCOS(theta))**2)/(y_r**2))
      
      temp = thetaA*x*x + thetaB*x*y + thetaC*y*y - (2*thetaA*x_c+thetaB*y_c)*x - (thetaB*x_c+2*thetaC*y_c)*y + &
              (thetaA*x_c*x_c+thetaB*x_c*y_c+thetaC*y_c*y_c)
    End If
    !=========LY modification, 2022-6-1=========

  END IF
!END DO
!temp = x**2+y**2-(pi/6.28)**2

END SUBROUTINE