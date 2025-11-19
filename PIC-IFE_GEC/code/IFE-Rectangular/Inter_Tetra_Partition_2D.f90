SUBROUTINE Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local, el_flag, vert, Intrs_pts, t_int_pt, p_int_pt)

IMPLICIT NONE

REAL(8), DIMENSION(:,:), INTENT(IN)	::	vert, Intrs_pts
INTEGER, DIMENSION(:), INTENT(IN)	::	tri_sign(2)
REAL(8), DIMENSION(:,:), POINTER 	::	p_int_pt, vertices
INTEGER, DIMENSION(:,:), POINTER	::	t_int_pt

INTEGER									n_Intrs_pts, n_tetra_pt, el_flag
INTEGER                                 pointer_reference_to_local(4)



n_Intrs_pts	= SIZE(Intrs_pts,2)

ALLOCATE(p_int_pt(2,4+n_Intrs_pts))
ALLOCATE(vertices(2,4))

vertices(:,1) = vert(:,pointer_reference_to_local(1))
vertices(:,2) = vert(:,pointer_reference_to_local(2))
vertices(:,3) = vert(:,pointer_reference_to_local(3))
vertices(:,4) = vert(:,pointer_reference_to_local(4))

p_int_pt(:,1:4)	= vertices(:,1:4)
p_int_pt(:,5)	= Intrs_pts(1,:)
p_int_pt(:,6)	= Intrs_pts(2,:)

IF (n_Intrs_pts==2) THEN
   
   n_tetra_pt	= 4
   ALLOCATE(t_int_pt(4,n_tetra_pt))
   
   IF (el_flag == 1) THEN
   ! T1 region   
   t_int_pt(1:3,1)	= (/1,	5,	6/)
   t_int_pt(4,1)	= 1

   ! T2 region
   t_int_pt(1:3,2)	= (/2,	6,	3/)
   t_int_pt(1:3,3)	= (/6,	3,	5/)
   t_int_pt(1:3,4)	= (/3,	5,	4/)

   t_int_pt(4,2:4)	= 2
   ELSEIF (el_flag == 2) THEN

   t_int_pt(1:3,1)	= (/1,	6,	4/)
   t_int_pt(1:3,2)	= (/4,	5,	6/)

   t_int_pt(4,1:2)	= 1

   ! T2 region

   t_int_pt(1:3,3)	= (/6,	3,	5/)
   t_int_pt(1:3,4)	= (/3,	6,	2/)

   t_int_pt(4,3:4)	= 2

   ENDIF
   
END IF

DEALLOCATE(vertices)    !$ ab.ZWZ 2021/9/1

END
