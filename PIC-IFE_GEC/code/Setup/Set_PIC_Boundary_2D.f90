SUBROUTINE Set_PIC_Boundary_2D(nnx,nny)

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: Set_PIC_Boundary_2D                                C
!
!  Purpose: Set PIC boundary and boundary type
!                                                                      C
!  Reviewer: Yuchuan Chu                            Date: 17-May-2012  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

USE Domain_2D
USE Particle_2D
USE Boundary_Data_2D
USE IFE_Data

IMPLICIT NONE

INTEGER                                                        ::  N_Boundary
INTEGER                                                        ::  i, j
INTEGER                                                        ::  nnx, nny
INTEGER                                                        ::  n_x, n_y, n_tot, n_tot_old
REAL(8)                                                        ::  x_start, y_start, x_end, y_end
INTEGER                                                        ::  num_of_interface_element

N_Boundary = SIZE(boundaries,1)
num_of_interface_element = SIZE(information_2,2)

IF( N_Boundary /= 0) THEN
    WRITE(6,*)
    WRITE(6,*)'Number of Boundary=',N_Boundary

    DO i = 1, N_Boundary
        WRITE(6,*)
        WRITE(6,*) 'boundaries(i)%Shape	            =',	boundaries(i)%Shape
        WRITE(6,*) 'boundaries(i)%Axis			    =', boundaries(i)%Axis
        WRITE(6,*) 'boundaries(i)%Length	        =', boundaries(i)%Length  
        WRITE(6,*) 'boundaries(i)%Locations(1,:)	=', boundaries(i)%Locations(1,:)
        WRITE(6,*) 'boundaries(i)%Locations(2,:)   =',	boundaries(i)%Locations(2,:) 
        WRITE(6,*) 'boundaries(i)%BoundType			=',	boundaries(i)%BoundType	 
        WRITE(6,*)
    END DO
ENDIF

n_x = nnx - 1
n_y = nny - 1
                     
!ALLOCATE(collect(ispe_tot, 0:n_x, 0:n_y))
!ALLOCATE(collectq(0:n_x,0:n_y))
!collectq = 0

!ALLOCATE(collect(ispe_tot, 6, num_of_interface_element))
!ALLOCATE(collectq(6,num_of_interface_element))
!collectq = 0
!
!DO j=1, ispe_tot
!    DO i=1, num_of_interface_element
!        collect(j,1:4,i) = information_2(3:6,i)
!        collectq(1:4,i) = information_2(3:6,i)
!    ENDDO
!ENDDO

END
