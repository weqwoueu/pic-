SUBROUTINE SetObjects_2D

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: SetObjects_2D.f90                                  C
!
!  Purpose: Setup the thruster object geometry and boundary conditions
!                                                                      C
!  Reviewer: Yuchuan Chu                              Date: 06-May-12  C
!  Comments: modified for normalization data input in input_2D.f90     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
USE PIC_MAIN_PARAM_2D
USE Domain_2D
USE Field_2D
USE Object_2D
!=========LY modification, 2022-5-19=========
Use IFE_Data, Only: HP, node_index, PhiPeriodicalNode
!=========LY modification, 2022-5-19=========

IMPLICIT NONE

INTEGER		iobject, i, j, k, idim
REAL(8)		x, r_temp

REAL(8), DIMENSION(2)	::	Xp, Xpo, Xt, Xb, Xc, X1, X2, Vn, Xtc
REAL(8), DIMENSION(2)	::	Xpc, Xp1, Xtb, Xpb, Xpt, XpbXtb
REAL(8)		Xtb2, XpbXtb2, XpbDotXtb, XptDotXtb, Xp1DotVn, Xpc2, XpcDotXtc, XpcDotXpc
INTEGER		coord1, coord2, coord3, Idummy
!=========LY modification, 2022-6-13=========
Real(8)	:: xx, yy
Real(8) :: r, x_c, y_c, x_r, y_r
Real(8) :: Fx1, Fy1, Fx2, Fy2
Real(8) :: temp, slope, theta, thetaA, thetaB, thetaC
Integer :: ni, nj, nindex
!=========LY modification, 2022-6-13=========
!=========LY modification for Phi-Periodical case, 2022-8-11=========
Integer :: count
Integer, Dimension(:), Allocatable :: PhiPeriodicalNodeTemp
Allocate(PhiPeriodicalNodeTemp(Size(HP,2)))
!=========LY modification for Phi-Periodical case, 2022-8-11=========

WRITE(6,*)
WRITE(6,*) 'SetObjects'

! There are MaxObjects

! 1) READ in the geometry information and phi info

! Total number of grids
WRITE(6,*) 'Total # of objects', MaxObjects

IF (MaxObjects > MaxTotObjects) THEN
	WRITE(6,*) 'Wrong Input. MaxTotObjects = ', MaxTotObjects
	STOP
END IF

DO iobject = 1, MaxObjects
	WRITE(6,*) 'Object =', iobject
	
	WRITE(6,*) '      Shape = ', shape_object(iobject)
	WRITE(6,*) '      Axis  = ', axis_object(iobject)

	WRITE(6,*) '      Dimension    = ', dimension_object(iobject,:)
  
	IF(shape_object(iobject)>1) THEN
		radius2_object(iobject) = dimension_object(iobject,2)**2
	ELSE
		radius2_object(iobject) = dimension_object(iobject,1)**2
	END IF		

	WRITE(6,*) '      Left Loc  = ', f_object_l(iobject,:)

	WRITE(6,*) '      Right Loc = ', f_object_r(iobject,:)

	WRITE(6,*) '      Region	= ', region_object(iobject)

	WRITE(6,*) '      Potential = ', phi_object(iobject)

	WRITE(6,*) '   Permittivity = ', eps_object(iobject)


END DO


! 2) Setup the geometry of the object and bc

! Object boundary and i-j-k location 
DO iobject=1,MaxObjects
	WRITE(6,*) 'object #',iobject
	DO idim=1,2
		WRITE(6,*) 'idim=',idim
		IF(f_object_l(iobject,idim) <= f_left_wall(idim)) THEN
			WRITE(6,*) 'f_l lt left wall, reset '
			f_object_l(iobject,idim) = f_left_wall(idim)
		END IF
		IF(f_object_r(iobject,idim) >= f_right_wall(idim)) THEN
			WRITE(6,*) 'f_r gt right wall, reset '
			f_object_r(iobject,idim) = f_right_wall(idim)
		END IF

		i_object_l(iobject,idim) = INT((f_object_l(iobject,idim)-Vert_o(idim))*hxi(idim))
		i_object_r(iobject,idim) = INT((f_object_r(iobject,idim)-Vert_o(idim))*hxi(idim))

		WRITE(6,*) 'i-j',i_object_l(iobject,idim), i_object_r(iobject,idim)
	END DO
END DO

! Set solid inner grid
DO iobject = 1,MaxObjects

	IF	(axis_object(iobject)==1) THEN
		coord1 = 1
		coord2 = 2
	ELSEIF	(axis_object(iobject)==2) THEN
		coord1 = 2
		coord2 = 1
	ENDIF


	IF (region_object(iobject) < vacuumRegion) THEN
	! Electrically conducting object
    WRITE(6,*) 'There is Fix potential of objects!'

		IF (shape_object(iobject)==3) THEN
		! Box
      !=========LY modification, 2022-7-25=========
      Do i = 1, Size(HP,2)
				xx = HP(1, i)
				yy = HP(2, i)
				
        If (xx >= f_object_l(iobject,1) .AND. xx <= f_object_r(iobject,1) .AND. &
            yy >= f_object_l(iobject,2) .AND. yy <= f_object_r(iobject,2)) then
          Phi(i, 1) =  phi_object(iobject)
          
          !!=========LY modification for Phi-Periodical case, 2022-8-11=========
          !!We set the first object is the Phi-Periodical object, so we can use it's information via objects(1).
          !If (iobject == 1) Then
          !  PhiPeriodicalNodeTemp(count) = i
          !  count = count + 1
          !End If
          !!=========LY modification for Phi-Periodical case, 2022-8-11=========
          
        End If
      End Do
      !=========LY modification, 2022-7-25=========
      
      !=========Old Code=========
				!DO j=i_object_l(iobject,2), i_object_r(iobject,2)
					!DO i=i_object_l(iobject,1), i_object_r(iobject,1)
					!	i_grid_flag(i,j) = 1
					!	phi(i,j) = phi_object(iobject)
					!END DO
     !   END DO
      !=========Old Code=========
        
    ELSEIF (shape_object(iobject)==1) THEN	
		! Circle.
      !=========LY modification, 2022-7-25=========
      r = dimension_object(iobject, 1)
      x_c = f_object_l(iobject, 1)
      y_c = f_object_l(iobject, 2)
      
      Do i = 1, Size(HP,2)
        xx = HP(1, i)
        yy = HP(2, i)
        
        
        !=========Original Fix potential code=========
        If ((xx-x_c)**2+(yy-y_c)**2 - r**2 < (1.0D-12)) Then
          Phi(i, 1) = phi_object(iobject)
        End If
        !=========Original Fix potential code=========
      End Do
      !=========LY modification, 2022-7-25=========
    
      !=========Old Code=========
				!DO j = 0, ny+1
				!	DO i = 0, nx+1
				!		Xp = VertX(:,i,j)
				!		Xc = f_object_l(iobject,:)
				!		Xpc = Xp - Xc
    !
				!		CALL Vec_Dot_3_2D(Xpc, Xpc, XpcDotXpc)
				!		!IF (XpcDotXpc<=radius2_object(iobject)+GTol) THEN
				!		!	i_grid_flag(i,j) = 1
				!		!	phi(i,j) = phi_object(iobject)
				!		!ENDIF
    !                    
    !                    !!!! **************** bjw add 2019.9.30 *****************
    !                    IF ((XpcDotXpc<=radius2_object(iobject)+GTol) .AND.(direction_object(iobject) == -1)) THEN
				!			i_grid_flag(i,j) = 1
				!			phi(i,j) = phi_object(iobject)
    !                    ELSEIF ((XpcDotXpc>=radius2_object(iobject)+GTol) .AND.(direction_object(iobject) == -2)) THEN
    !                        i_grid_flag(i,j) = 1
				!			phi(i,j) = phi_object(iobject)
				!		ENDIF
    !                    !!!! **************** bjw add 2019.9.30 *****************
				!	END DO
    !    END DO
      !=========Old Code=========
        
    Elseif (shape_object(iobject)==4 .OR. shape_object(iobject)==6) Then  !LY modification for quadrangle or triangle, arch, 2022-6-13
      !=========LY modification, 2022-7-25=========
      Do i = 1, Size(HP,2)
        If (node_index(i) < vacuumRegion) Then  !There denote the node inside objects.
          Phi(i,1) = phi_object(iobject)
        End If
      End Do
      !=========LY modification, 2022-7-25=========
      
    Elseif (shape_object(iobject)==5) Then  !LY modification for ellipse, 2022-6-13
      x_r = dimension_object(iobject, 1)    !radius of x axis of the ellipse
      y_r = dimension_object(iobject, 2)    !radius of y axis of the ellipse
      
      Fx1 = f_object_l(iobject, 1)          !Focus1 x of the ellipse
      Fy1 = f_object_l(iobject, 2)          !Focus1 y of the ellipse
      Fx2 = f_object_r(iobject, 1)          !Focus2 x of the ellipse
      Fy2 = f_object_r(iobject, 2)          !Focus2 y of the ellipse
      x_c = (Fx1+Fx2) / 2.0                 !Centre x of the ellipse
      y_c = (Fy1+Fy2) / 2.0                 !Centre y of the ellipse
        
      !=========LY modification, 2022-7-25=========
      Do i = 1, Size(HP,2)
        xx = HP(1, i)
        yy = HP(2, i)
        
        !=========LY modification, 2022-6-1=========
        If (ABS(Fx1-Fx2)<1.0D-8 .OR. ABS(Fy1-Fy2)<1.0D-8) Then    !vertical or horizontal ellipse.
          temp = ((xx-x_c)/x_r)**2+((yy-y_c)/y_r)**2
        Else    !inclined ellipse.
          slope = (Fy2-Fy1) / (Fx2-Fx1)
          If (ABS(slope) < 1.0D-8) Then   !There denote slope=0, but we already deal with the condition in above.
            Write(*,*) 'The slope = ', slope
            Stop
          End If
          
          If (slope > 1.0D-8) Then
            theta = DATAN(slope)
          Elseif (slope < -1.0D-8) Then
            theta = DATAN(slope) + pi
          End If
          thetaA = (((DCOS(theta))**2)/(x_r**2) + ((DSIN(theta))**2)/(y_r**2))
          thetaB = (1/(x_r**2) - 1/(y_r**2)) * DSIN(2*theta)
          thetaC = (((DSIN(theta))**2)/(x_r**2) + ((DCOS(theta))**2)/(y_r**2))

          temp = thetaA*xx*xx + thetaB*xx*yy + thetaC*yy*yy - (2*thetaA*x_c+thetaB*y_c)*xx - (thetaB*x_c+2*thetaC*y_c)*yy + &
                (thetaA*x_c*x_c+thetaB*x_c*y_c+thetaC*y_c*y_c)
        End If
        
        If (temp-1.0 < -(1.0D-8)) Then
          Phi(i, 1) = phi_object(iobject)
        End If
        !=========LY modification, 2022-6-1=========
      End Do
      !=========LY modification, 2022-7-25=========
      
		ELSEIF (shape_object(iobject)==11.AND.axis_object(iobject)<=2.AND.axis_object(iobject)>=1) THEN		
!		! Rectangular Plate parallel to an axis
				DO j = 0, ny+1
					DO i = 0, nx+1
						Xp(1:2) = VertX(:,i,j)
						X1(1:2) = f_object_l(iobject,:)
						X2(1:2) = f_object_r(iobject,:)
						IF (	(X2(1)-Xp(1))*(Xp(1)-X1(1)) >= -GTol						.AND.	&
								(X2(2)-Xp(2))*(Xp(2)-X1(2)) >= -GTol						.AND.	&
								DABS((Xp(coord2)-X1(coord2))*(X2(coord1)-X1(coord1)))<= GTol	) THEN
							i_grid_flag(i,j) = 1
							phi(i,j) = phi_object(iobject)
						ENDIF
					END DO
				END DO
		END IF
	END IF

END DO

!=========LY modification for Phi-Periodical case, 2022-8-11=========
!We use the 'PhiPeriodicalNode' to store the node index in Phi-Periodical objects.
Allocate(PhiPeriodicalNode(count-1))
Do i = 1, count-1
  PhiPeriodicalNode(i) = PhiPeriodicalNodeTemp(i)
End Do
!=========LY modification for Phi-Periodical case, 2022-8-11=========
END


