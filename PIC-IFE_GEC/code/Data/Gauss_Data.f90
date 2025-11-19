MODULE	Gauss_Data

IMPLICIT NONE

!=======
REAL(8), DIMENSION(9)                   ::  Gauss_Coefficient_Reference_Nine
REAL(8), DIMENSION(9,2)                 ::  Gauss_Point_Reference_Nine

REAL(8), DIMENSION(9)                   ::  Gauss_Coefficient_Local_Nine
REAL(8), DIMENSION(9,2)                 ::  Gauss_Point_Local_Nine

REAL(8), DIMENSION(4)                   ::  Gauss_Coefficient_Reference_Four
REAL(8), DIMENSION(4,2)                 ::  Gauss_Point_Reference_Four

REAL(8), DIMENSION(4)                   ::  Gauss_Coefficient_Local_Four
REAL(8), DIMENSION(4,2)                 ::  Gauss_Point_Local_Four

!=======
REAL(8), DIMENSION(9)                   ::  Gauss_Coefficient_Reference_Triangle_Nine
REAL(8), DIMENSION(9,2)                 ::  Gauss_Point_Reference_Triangle_Nine

REAL(8), DIMENSION(9)                   ::  Gauss_Coefficient_Local_Triangle_Nine
REAL(8), DIMENSION(9,2)                 ::  Gauss_Point_Local_Triangle_Nine

REAL(8), DIMENSION(4)                   ::  Gauss_Coefficient_Reference_Triangle_Four
REAL(8), DIMENSION(4,2)                 ::  Gauss_Point_Reference_Triangle_Four

REAL(8), DIMENSION(4)                   ::  Gauss_Coefficient_Local_Triangle_Four
REAL(8), DIMENSION(4,2)                 ::  Gauss_Point_Local_Triangle_Four

REAL(8), DIMENSION(3)                   ::  Gauss_Coefficient_Reference_Triangle_Three
REAL(8), DIMENSION(3,2)                 ::  Gauss_Point_Reference_Triangle_Three

REAL(8), DIMENSION(3)                   ::  Gauss_Coefficient_Local_Triangle_Three
REAL(8), DIMENSION(3,2)                 ::  Gauss_Point_Local_Triangle_Three

!----------------------- add for DG---------------
REAL(8), DIMENSION(8)                   ::  Gauss_Coefficient_Local_1D_Linear_Eight
REAL(8), DIMENSION(2,8)                 ::  Gauss_Point_Local_1D_Linear_Eight

REAL(8), DIMENSION(4)                   ::  Gauss_Coefficient_Local_1D_Linear_Four
REAL(8), DIMENSION(2,4)                 ::  Gauss_Point_Local_1D_Linear_Four

REAL(8), DIMENSION(2)                   ::  Gauss_Coefficient_Local_1D_Linear_Two
REAL(8), DIMENSION(2,2)                 ::  Gauss_Point_Local_1D_Linear_Two
!----------------------- add for DG---------------

!=======
REAL(8), DIMENSION(4)					::	Gauss_coefficient_reference_1D_Four
REAL(8), DIMENSION(4)					::	Gauss_point_reference_1D_Four	

REAL(8), DIMENSION(4)					::	Gauss_coefficient_Local_1D_Four
REAL(8), DIMENSION(4)					::	Gauss_point_Local_1D_Four					

REAL(8), DIMENSION(8)					::	Gauss_coefficient_reference_1D_Eight
REAL(8), DIMENSION(8)					::	Gauss_point_reference_1D_Eight	

REAL(8), DIMENSION(8)					::	Gauss_coefficient_Local_1D_Eight
REAL(8), DIMENSION(8)					::	Gauss_point_Local_1D_Eight

REAL(8), DIMENSION(2)					::	Gauss_coefficient_reference_1D_Two
REAL(8), DIMENSION(2)					::	Gauss_point_reference_1D_Two	

REAL(8), DIMENSION(2)					::	Gauss_coefficient_Local_1D_Two
REAL(8), DIMENSION(2)					::	Gauss_point_Local_1D_Two


END MODULE