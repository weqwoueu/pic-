MODULE Object_Data_2D

! Purpose:		Object Data Module
! Last Update:	4/19/2004 05:09 PM

!!!!!!!!!!!!!!!!!!!!!!!chj add for erosion!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, PARAMETER	::	MaxNumOfWallsPerObject = 1
! Wall Data Structure
TYPE WallType
	INTEGER Shape
	! = 1 Polygonal line
	! = 2 Closed curve
	INTEGER Channelwall
	! = 0 // (for a Closed curve)
	! = 3 // upwall of the Hall channel
	! = 4 // downwall of the Hall channel
	REAL(8) Limits(2,2)
	! Locations(1,:) = (Radius, yc),   for Circle start point
	! Locations(1,:) = : Understructure(xmin, ymin)		for Rectangular Parallelopiped
	! Locations(2,:) = : Understructure(xmax, ymax)		for Rectangular Parallelopiped
	REAL(8)	stepx
	! stepx must be an aliquot part of either (xmax-xmin) or  (ymax-ymin)
	INTEGER	nxp
	REAL(8), DIMENSION(:,:), POINTER	::	node
END TYPE
!!!!!!!!!!!!!!!!!!!!!!!chj add for erosion!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Object Data Structure
TYPE ObjectType

	INTEGER	Shape
	! = 1 : Circle
	! = 2 : Understructure
	! = 3 : Box
	INTEGER	Axis
	! = 0	No axis (for a Circle this means complete Circle)
	! = 1 // x-axis	 : Understructure
	! = 2 // y-axis	 : Understructure
	! = 3 // z-axis	 : Understructure
	! = 4 _|_ x-axis : Understructure
	! = 5 _|_ y-axis : Understructure
	! = 6 _|_ z-axis : Understructure
	REAL(8)	Dimensions(2)
	! Dimensions(1)		= (r)			for Circle radius
	! Dimensions(1:2)	= : Understructure(r1, r2)		for Circular Right Cylinder
	! Dimensions(:)		= : Understructure(Lx, Ly)	for Rectangular Parallelopiped
	!REAL(8)	Locations(2,2)
  !=========LY modification, 2022-6-13=========
  Real(8) Locations(4,2)
  !=========LY modification, 2022-6-13=========
  
	! Locations(1,:) = (xc, yc),   for Circle center
	! Locations(1,:) = : Understructure(xc1, yc1, zc1), 
	! Locations(2,:) = : Understructure(xc2, yc2, zc2)	for Circular Right Cylinder
	! Locations(2,1) = : Understructurebase loaction	for Sphere (cut by a plane parallel to one of the axes planes and given by Axis)
	! Locations(1,:) = : Understructure(xmin, ymin, zmin)		for Rectangular Parallelopiped
	! Locations(2,:) = : Understructure(xmax, ymax, zmax)		for Rectangular Parallelopiped

	INTEGER	Regions(2)
	! = (inner region index, outer region index)
	! -1			vacuum
	! -2, -3, ...	physical objects
    
    !!!! **************** bjw add 2019.9.30 *****************
    INTEGER	Direction
	! -1			inner of Box or Cicle
	! -2, -3, ...	outer of Box or Cicle
    !!!! **************** bjw add 2019.9.30 *****************
    
	REAL(8) Phi
	! OBJECT POTENTIAL

	REAL(8) Eps
	! OBJECT EPSILON
	
	INTEGER	Erosion
	! = 0 : stabilization
	! = 1 : evolution

	TYPE(WallType) Wall(MaxNumOfWallsPerObject)

END TYPE

!INTEGER, PARAMETER	::	vacuum = -1
INTEGER		::	vacuum 

!TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
!TYPE(ObjectType), DIMENSION(:), ALLOCATABLE	:: objects

END MODULE