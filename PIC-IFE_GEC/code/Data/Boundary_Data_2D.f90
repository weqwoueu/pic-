MODULE Boundary_Data_2D

INTEGER, PARAMETER	::	MaxNumOfWallsPerObject = 1
! Wall Data Structure
TYPE WallType
	INTEGER Shape
	! = 1 Polygonal line
	! = 2 Closed curve
	INTEGER Channelwall
	! = 0 // (for a Closed curve)
	! = 1 // dostream innerwall of the Hall channel
	! = 2 // dostream outerwall of the Hall channel
	! = 3 // innerwall of the Hall channel
	! = 4 // outerwall of the Hall channel
	INTEGER	nxp
	REAL(8), DIMENSION(:,:), POINTER	::	node
END TYPE


! Boundary Data Structure
TYPE BoundaryType

	INTEGER	Shape
	! = 1 : vertical line or horizontal line
	! = 2 : line
	! = 3 : circle
	INTEGER	Axis
	! = 0	No axis
	REAL(8)	Length
	! straight lines means length, circle means radius
	REAL(8)	Locations(2,2)
	! Locations(1,:) = (xc, yc),   for Circle center
	! Locations(1,:) = : start point
	! Locations(2,:) = : end point
    INTEGER    BoundType
    ! Boundary Type
    
    INTEGER    Obound
    ! = 0 : not object boundary
	! > 0 :  object index with this boundary
    
	INTEGER	Erosion
	! = 0 : stabilization
	! = 1 : evolution
	
	TYPE(WallType) Wall(MaxNumOfWallsPerObject)              

END TYPE

END MODULE