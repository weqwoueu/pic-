SUBROUTINE Function_True(delta, index, x, y, derivative_degree_x, derivative_degree_y, r)
!!! bjw add
USE IFE_Data   
IMPLICIT NONE
INTEGER, INTENT(IN) :: delta !$ ab.ZWZ 2021/7/7
REAL(8)     x, y, r, index, r0
INTEGER     derivative_degree_x, derivative_degree_y


!r0 = 0.31415926
r0 = 0.500253607
!r0 = 0.4

IF (delta == 0) THEN    
    !IF (index == Global_Beta(1)) THEN   !!! outside
    IF (x*x+y*y - r0*r0 > 0.0) THEN   !!! outside  
    
        IF (derivative_degree_x == 0 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(1) * (x**2 + y**2)
            !r = EXP(x*x+y*y) / 10.0 !Global_Beta(1)
            !r = ((x*x+y*y)**(1.5))/10.0 + (0.9)*(r0**3)
            r = ((x*x+y*y)**(2.0))/10.0
    
        ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(1) * (2*x)
            !r = 2.0*x*EXP(x*x+y*y) / 10.0 !Global_Beta(1)
            !r = (0.3)*x*((x*x+y*y)**(0.5))
            r = 4.0*x*(x*x+y*y)/10.0
    
        ELSEIF (derivative_degree_x == 0 .AND. derivative_degree_y ==1 ) THEN
    
            !r = 1.0/Global_Beta(1) * (2*y)
            !r = 2.0*y*EXP(x*x+y*y) / 10.0 !Global_Beta(1)
            !r = (0.3)*y*((x*x+y*y)**(0.5))
            r = 4.0*y*(x*x+y*y)/10.0
    
        END IF
    
    !ELSEIF (index == Global_Beta(2)) THEN  !!!inside
    ELSEIF (x*x+y*y - r0*r0 <= 0.0) THEN   !!! inside    
    
        IF (derivative_degree_x == 0 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(2) * (x**2 + y**2) + &
            !   (1.0/Global_Beta(1) - 1.0/Global_Beta(2))* r0 * r0
            !r = EXP(x*x+y*y) / Global_Beta(2) + (1.0/Global_Beta(1) - 1.0/Global_Beta(2)) * EXP(r0 * r0)
            !r = EXP(x*x+y*y) / 1.0 + (1.0/10.0 - 1.0/1.0) * EXP(r0 * r0)
            !r = ((x*x+y*y)**(1.5))
            r = ((x*x+y*y)**(2.0)) - (0.9)*(r0**4)
    
        ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(2) * (2*x)
            !r = 2.0*x*EXP(x*x+y*y) / 1.0 !Global_Beta(2)
            !r = 3*x*((x*x+y*y)**(0.5))
            r = 4.0*x*(x*x+y*y)
    
        ELSEIF (derivative_degree_x == 0 .AND. derivative_degree_y ==1 ) THEN
    
            !r = 1.0/Global_Beta(2) * (2*y)
            !r = 2.0*y*EXP(x*x+y*y) / 1.0 !Global_Beta(2)
            !r = 3*y*((x*x+y*y)**(0.5))
            r = 4.0*y*(x*x+y*y)
    
        END IF
    
    END IF
ELSEIF(delta == 1) THEN
    IF (x*x+(y-1)*(y-1) - r0*r0 > 0.0) THEN   !!! outside  
    
        IF (derivative_degree_x == 0 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(1) * (x**2 + (y-1)**2)
            !r = EXP(x*x+(y-1)*(y-1)) / 10.0 !Global_Beta(1)
            r = EXP(x*x+(y-1)*(y-1)) / Global_Beta(1)
    
        ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(1) * (2*x)
            !r = 2.0*x*EXP(x*x+(y-1)*(y-1)) / 10.0 !Global_Beta(1)
            r = 2.0*x*EXP(x*x+(y-1)*(y-1)) / Global_Beta(1)
    
        ELSEIF (derivative_degree_x == 0 .AND. derivative_degree_y ==1 ) THEN
    
            !r = 1.0/Global_Beta(1) * (2*(y-1))
            !r = 2.0*(y-1)*EXP(x*x+(y-1)*(y-1)) / 10.0 !Global_Beta(1)
            r = 2.0*(y-1)*EXP(x*x+(y-1)*(y-1)) / Global_Beta(1)
    
        END IF
    
    !ELSEIF (index == Global_Beta(2)) THEN  !!!inside
    ELSEIF (x*x+(y-1)*(y-1) - r0*r0 <= 0.0) THEN   !!! inside    
    
        IF (derivative_degree_x == 0 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(2) * (x**2 + (y-1)**2) + &
            !   (1.0/Global_Beta(1) - 1.0/Global_Beta(2))* r0 * r0
            !r = EXP(x*x+(y-1)*(y-1)) / 1.0 + (1.0/10.0 - 1.0/1.0) * EXP(r0 * r0)
            r = EXP(x*x+(y-1)*(y-1)) / Global_Beta(2) + (1.0/Global_Beta(1) - 1.0/Global_Beta(2)) * EXP(r0 * r0)
    
        ELSEIF (derivative_degree_x == 1 .AND. derivative_degree_y ==0 ) THEN
    
            !r = 1.0/Global_Beta(2) * (2*x)
            !r = 2.0*x*EXP(x*x+(y-1)*(y-1)) / 10.0 !Global_Beta(2)
            r = 2.0*x*EXP(x*x+(y-1)*(y-1)) / Global_Beta(2)
    
        ELSEIF (derivative_degree_x == 0 .AND. derivative_degree_y ==1 ) THEN
    
            !r = 1.0/Global_Beta(2) * (2*(y-1))
            !r = 2.0*(y-1)*EXP(x*x+(y-1)*(y-1)) / 10.0 !Global_Beta(2)
            r = 2.0*(y-1)*EXP(x*x+(y-1)*(y-1)) / Global_Beta(2)
    
        END IF
    
    END IF
ENDIF

!WRITE(*,*) 'Function_True', index, x, y, derivative_degree_x, derivative_degree_y, r

END