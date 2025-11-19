SUBROUTINE generate_Tier_initial_2D(P, T, CellMesh)
!LY Add for Multigrid Store, 2022-1-17
  
USE Cell_Data_2D
IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER :: P
INTEGER, DIMENSION(:,:), POINTER :: T
TYPE(CellDataType), DIMENSION(:), POINTER :: CellMesh

INTEGER :: n_element, i, j

n_element = SIZE(T,2)
!The variable "n_element" denote number of initial element without refinement.

ALLOCATE(CellMesh(n_element))
DO i = 1, n_element
  CellMesh(i)%Globalindex = 0
  CellMesh(i)%Localindex = 0
  CellMesh(i)%isSplitted = 0
  CellMesh(i)%Boundary(1:4) = 0.0
  CellMesh(i)%Tier = 0
  CellMesh(i)%Parent = 0
  CellMesh(i)%Child(1:4) = 0
  CellMesh(i)%Finalindex = 0
  CellMesh(i)%Nodeindex(1:4) = 0
  CellMesh(i)%FinalParent = 0 !LY modification for Speeding Particle Positioning, 2022-7-25
END DO

DO i = 1, n_element
  
  CellMesh(i)%Localindex = i
  
  IF (T(5,i) > 0) THEN
    CellMesh(i)%isSplitted = 1    !DG element and it need refine.
  END IF
  
  CellMesh(i)%Boundary(1) = P(1, T(1,i))    !Cell left boundary.
  CellMesh(i)%Boundary(2) = P(1, T(2,i))    !Cell right boundary.
  CellMesh(i)%Boundary(3) = P(2, T(1,i))    !Cell bottom boundary.
  CellMesh(i)%Boundary(4) = P(2, T(4,i))    !Cell top boundary.
  
  CellMesh(i)%Tier = 1    !Every element is locate on 1st Tier on initial mesh.
  
  CellMesh(i)%Finalindex = i   !If no mesh refine, the element final index is themselves.
  
  DO j = 1, 4
    CellMesh(i)%Nodeindex(j) = T(j, i)   !We store four node global index of element.
  END DO
  
  !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
  CellMesh(i)%FinalParent = i
  !=========LY modification for Speeding Particle Positioning, 2022-7-25=========
END DO

END SUBROUTINE generate_Tier_initial_2D