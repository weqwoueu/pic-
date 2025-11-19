MODULE Cell_Data_2D
  
!Purpose : Cell Data Module for Multigrid Store on Local Mesh Refine.
!Author : Li Yang, 2022-1-17
  
IMPLICIT NONE

!Cell Data Structure
TYPE CellDataType
  
  INTEGER :: Globalindex
  !Global index of element on parallel computation.

  INTEGER :: Localindex
  !Index of element on Cell Structure, no parallel.

  INTEGER :: isSplitted
  !Refine flag of element.
  ! = 1 : Yes
  ! = 0 : No

  REAL(8) :: Boundary(4)
  !Boundary information of element.
  !Boundary(1) = Left
  !Boundary(2) = Right
  !Boundary(3) = Bottom
  !Boundary(4) = Top

  INTEGER :: Tier
  !Refine layer flag of element.
  ! = 1 : Initial mesh, no refine.
  ! = 2 : 1st refine mesh.
  ! = 3 : 2st refine mesh.
  !etc.

  INTEGER :: Parent
  !The element's parent element index on above mesh structure.

  INTEGER :: Child(4)
  !If the element is refined or splitted, then it will generate four child element.
  !Child(1) = left_bottom :  left and bottom chile element.
  !Child(2) = left_top :     left and top chile element.
  !Child(3) = right_bottom : right and bottom chile element.
  !Child(4) = right_top :    right and top chile element.

  INTEGER :: Finalindex
  !When mesh refinement finish, the final index of element on end mesh.(HT index)
  
  INTEGER :: Nodeindex(4)
  
  !=========LY modification for Speeding of Particle positioning, 2022-7-25========
  INTEGER :: FinalParent
  !The 'FinalParent' denote the most original Parent element index in initial mesh.
  !=========LY modification for Speeding of Particle positioning, 2022-7-25========
  
END TYPE CellDataType

END MODULE Cell_Data_2D