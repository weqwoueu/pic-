SUBROUTINE generate_DG_edge_flag(DGE, dimensions, DGP, DGT, DG_flag)
!----------------------- wsy add DG_flag 2021/7/30 -------------------------------------
!第一行为边的边界条件 0为内部边 -1为Dirichlet -2为Neuman -3为Robin
!第二行为边的外法向量的x坐标
!第三行为边的外法向量的y坐标
!第四行为边的跃度(Dirichlet边界的跃度为0)

IMPLICIT NONE

REAL :: top, bottom, left, right
Integer :: number_of_edges, n
INTEGER, DIMENSION(:,:), POINTER	::  DG_flag, DGE, DGT
REAL(8), DIMENSION(:,:), POINTER	::  DGP
REAL(8), DIMENSION(2,2), INTENT(IN)		::	dimensions

top = dimensions(2,2)
bottom = dimensions(2,1)
left = dimensions(1,1)
right = dimensions(1,2)
number_of_edges = SIZE(DGE, 2)
ALLOCATE(DG_flag(4, number_of_edges))

DO n = 1, number_of_edges
    IF (DGP(1,DGE(1,n)) == left .AND. DGP(1,DGE(2,n)) == left) THEN
        DG_flag(1,n)=1
        DG_flag(2,n)=-1
        DG_flag(3,n)=0
        DG_flag(4,n)=0
    ELSE IF (DGP(1,DGE(1,n)) == right .AND. DGP(1,DGE(2,n)) == right) THEN
        DG_flag(1,n)=1
        DG_flag(2,n)=1
        DG_flag(3,n)=0
        DG_flag(4,n)=0
    ELSE IF (DGP(2,DGE(1,n)) == bottom .AND. DGP(2,DGE(2,n)) == bottom) THEN
        DG_flag(1,n)=1
        DG_flag(2,n)=0
        DG_flag(3,n)=-1
        DG_flag(4,n)=0
    ELSE IF (DGP(2,DGE(1,n)) == top .AND. DGP(2,DGE(2,n)) == top) THEN
        DG_flag(1,n)=1
        DG_flag(2,n)=0
        DG_flag(3,n)=1
        DG_flag(4,n)=0
    ELSE IF (DGP(1,DGE(1,n))==DGP(1,DGE(2,n))) THEN ! vertical inner edge
        DG_flag(1,n)=0
        DG_flag(2,n)=-1
        DG_flag(3,n)=0
        IF (DGP(1,DGE(1,n))==DGP(1,DGT(2,DGE(5,n)))) THEN
            DG_flag(4,n)=-1
        ELSE IF (DGP(1,DGE(1,n))==DGP(1,DGT(1,DGE(5,n)))) THEN
            DG_flag(4,n)=1
        END IF
    ELSE IF (DGP(2,DGE(1,n))==DGP(2,DGE(2,n))) THEN ! horizontal inner edge
        DG_flag(1,n)=0
        DG_flag(2,n)=0
        DG_flag(3,n)=-1
        IF (DGP(2,DGE(1,n))==DGP(2,DGT(1,DGE(5,n)))) THEN
            DG_flag(4,n)=1
        ELSE IF (DGP(2,DGE(1,n))==DGP(2,DGT(4,DGE(5,n)))) THEN
            DG_flag(4,n)=-1
        END IF
    END IF

END DO

END SUBROUTINE