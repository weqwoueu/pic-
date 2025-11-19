SUBROUTINE SetEleRef_2D

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: SetEleRef_2D.f90                                   C
!
!  Purpose: Set electron density reference
!                                                                      C
!  Reviewer: Yuchuan Chu                              Date: 06-May-12  C
!  Comments: modified for normalization data input in input_2D.f90     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

USE Field_2D

IMPLICIT NONE

WRITE(6,*)
WRITE(6,*) 'SetEleRef'

! READ electron  reference
WRITE(6,*) 'plasma reference'
WRITE(6,*) 'den0_ref=',den0_ref,' Te_ref=',Te_ref,' phi0_ref=',phi0_ref


END










