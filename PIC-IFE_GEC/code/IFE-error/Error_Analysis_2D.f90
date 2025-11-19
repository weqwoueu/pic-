SUBROUTINE Error_Analysis_2D(delta, nnx)
!! bjw add    
!! Phi----silulate results    

!USE IFE_Param
USE IFE_INTERFACE, ONLY: Generate_Gauss_reference_triangle, Generate_Gauss_reference, &
                         Compute_IFE_Norm_error
IMPLICIT NONE
INTEGER, INTENT(IN) :: delta, nnx !$ mb.ZWZ 2021/7/7
INTEGER     ::	it
CHARACTER*20	fname, sname 
CHARACTER*30    filename

INTEGER                                     Gauss_point_number
REAL(8),DIMENSION(:),POINTER            ::  Gauss_coefficient_reference
REAL(8),DIMENSION(:,:),POINTER          ::  Gauss_point_reference
INTEGER                                     Gauss_point_number_triangle
REAL(8),DIMENSION(:),POINTER            ::  Gauss_coefficient_reference_triangle
REAL(8),DIMENSION(:,:),POINTER          ::  Gauss_point_reference_triangle

REAL(8)      :: error_L2, error_H1x, error_H1y, error_H1

Gauss_point_number = 9
Gauss_point_number_triangle = 9

ALLOCATE(Gauss_coefficient_reference(Gauss_point_number), Gauss_point_reference(Gauss_point_number,2))
ALLOCATE(Gauss_coefficient_reference_triangle(Gauss_point_number_triangle), &
         Gauss_point_reference_triangle(Gauss_point_number_triangle,2))

CALL Generate_Gauss_reference( Gauss_point_number, Gauss_coefficient_reference, Gauss_point_reference)
CALL Generate_Gauss_reference_triangle( Gauss_point_number_triangle, Gauss_coefficient_reference_triangle, &
                                        Gauss_point_reference_triangle)


CALL Compute_IFE_Norm_error(delta, Gauss_coefficient_reference, Gauss_point_reference, Gauss_point_reference_triangle, &
                                       Gauss_coefficient_reference_triangle, 0, 0, error_L2)   !!! L2

CALL Compute_IFE_Norm_error(delta, Gauss_coefficient_reference, Gauss_point_reference, Gauss_point_reference_triangle, &
                                       Gauss_coefficient_reference_triangle, 1, 0, error_H1x)   !!! H1x

CALL Compute_IFE_Norm_error(delta, Gauss_coefficient_reference, Gauss_point_reference, Gauss_point_reference_triangle, &
                                       Gauss_coefficient_reference_triangle, 0, 1, error_H1y)   !!! H1y

error_H1 = DSQRT(error_H1x**2 + error_H1y**2)

!WRITE(fname,999) it
WRITE(sname,999) nnx-1
999 FORMAT(I4.4)
!filename = 'Norm_Error_'//TRIM(sname)//'_'//TRIM(fname)//'.txt'
filename = 'Norm_Error_'//TRIM(sname)//'.txt'
OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))

!OPEN(1, ACTION = 'WRITE', FILE = 'Norm_Error.txt')
    WRITE(1,*) 'Number of Element =', nnx-1
    !WRITE(1,*) 'Number of Time =', it
    WRITE(1,*) 'L2 Error', 'H1 Error'
    WRITE(1,*) error_L2, error_H1
CLOSE(1)

WRITE(*,*) error_L2, error_H1

END