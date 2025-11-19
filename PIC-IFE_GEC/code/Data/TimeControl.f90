MODULE TimeControl

IMPLICIT NONE

REAL(8)										dt
INTEGER(4)									nt, n_updatefld,Erostart
INTEGER(4)									n_pdiag, n_engydiag, n_gdiag, n_dump, ilap, irestart
INTEGER(4), DIMENSION(:), ALLOCATABLE	::	n_stride
INTEGER(4)									index_x, index_y, index_z 
REAL(8)										xcenter, ycenter, zcenter 
!> ab.ZWZ for RF power source
Real(8)                                 ::  Frequency
Real(8)                                 ::  Period, Period_res
Integer(4)                              ::  Pstep, AverStep
Integer(4)                              ::  nPeriod = 0, iStep_previous = 0
END MODULE