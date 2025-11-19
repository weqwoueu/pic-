Subroutine Line_Ellipse_Intersection_BL(line_ends, Ellipse_radius, Ellipse_Focus, inters_point, inters_flag)
  
  Use IFE_MAIN_PARAM
  Implicit None
  
  Real(8), Intent(In)		::	line_ends(2,2), Ellipse_radius(2), Ellipse_Focus(2,2)
  Real(8), Intent(Out)	::	inters_point(2)
  Integer, Intent(Out)	::	inters_flag
  
  Real(8), Dimension(2) :: root1, root2
  Real(8) :: x_r, y_r, x_c, y_c, x1, y1, x2, y2
  Real(8) :: rx1, ry1, rx2, ry2, r
  Integer :: root1_flag, root2_flag
  Real(8) :: delta
  Integer :: n_root
  !=========LY modification, 2022-5-30=========
  Real(8) :: Fx1, Fy1, Fx2, Fy2
  !=========LY modification, 2022-5-30=========
  !=========LY modification, 2022-6-1=========
  Real(8) :: slope, theta, thetaA, thetaB, thetaC
  Real(8) :: deltaA, deltaB, deltaC
  !=========LY modification, 2022-6-1=========

  !=========Part 0: Initilization=========
  inters_point = (/0.,0./)
  inters_flag = 0
  
  x_r = Ellipse_radius(1)
  y_r = Ellipse_radius(2)
  
  !=========LY modification, 2022-5-30=========
  Fx1 = Ellipse_Focus(1,1)
  Fy1 = Ellipse_Focus(1,2)
  Fx2 = Ellipse_Focus(2,1)
  Fy2 = Ellipse_Focus(2,2)
  
  x_c = (Fx1+Fx2) / 2.0
  y_c = (Fy1+Fy2) / 2.0
  !=========LY modification, 2022-5-30=========
  
  x1 = line_ends(1,1)
  y1 = line_ends(1,2)

  x2 = line_ends(2,1)
  y2 = line_ends(2,2)
  
  root1 = 0.
  root2 = 0.

  root1_flag = 0
  root2_flag = 0
  !=========Part 0: Initilization=========
  
  !=========Part 1: Calculate intersection node=========
  !=========LY modification, 2022-6-1=========
  If (ABS(Fx1-Fx2)<1.0D-8 .OR. ABS(Fy1-Fy2)<1.0D-8) Then    !vertical or horizontal ellipse.
    If (ABS(x1-x2)<1.0D-8) Then   !vertical edge.
      delta = x_r*x_r -(x1-x_c)*(x1-x_c)
      If (delta < -1.0D-8) Then   !No intersection node.
        n_root = 0
      Elseif (delta > 1.0D-8) Then  !2 intersection nodes.
        n_root = 2
        root1(1) = x1
        root1(2) = y_c + y_r*DSQRT(1-((x1-x_c)/x_r)**2)
        root2(1) = x1
        root2(2) = y_c - y_r*DSQRT(1-((x1-x_c)/x_r)**2)
      Elseif (ABS(delta) < 1.0D-8) Then   !1 intersection node.--tangent
        n_root = 1
        root1(1) = x1
        root1(2) = y_c
        root2(1) = x1
        root2(2) = y_c
      End If
    End If

    If (x1/=x2 .AND. ABS(y1-y2)<1.0D-8) Then  !horizon edge.
      delta = y_r*y_r -(y1-y_c)*(y1-y_c)
      If (delta < -1.0D-8) Then   !No intersection node.
        n_root = 0
      Elseif (delta > 1.0D-8) Then  !2 intersection nodes.
        n_root = 2
        root1(1) = x_c + x_r*DSQRT(1-((y1-y_c)/y_r)**2)
        root1(2) = y1
        root2(1) = x_c - x_r*DSQRT(1-((y1-y_c)/y_r)**2)
        root2(2) = y1
      Elseif (ABS(delta) < 1.0D-8) Then   !1 intersection node.--tangent
        n_root = 1
        root1(1) = x_c
        root1(2) = y1
        root2(1) = x_c
        root2(2) = y1
      End If
    End If
    
  Else    !inclined ellipse
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
    
    If (ABS(x1-x2) < 1.0D-8) Then   !vertical edge.
      deltaA = thetaC
      deltaB = thetaB*x1 - thetaB*x_c - 2*thetaC*y_c
      deltaC = thetaA*x1*x1 - (2*thetaA*x_c+thetaB*y_c)*x1 + (thetaA*x_c*x_c+thetaB*x_c*y_c+thetaC*y_c*y_c-1)
      
      delta = deltaB**2 - 4*deltaA*deltaC
      If (delta < -1.0D-8) Then   !No intersection node.
        n_root = 0
      Elseif (delta > 1.0D-8) Then  !2 intersection nodes.
        n_root = 2
        root1(1) = x1
        root1(2) = (-deltaB+DSQRT(delta)) / (2*deltaA)
        root2(1) = x1
        root2(2) = (-deltaB-DSQRT(delta)) / (2*deltaA)
      Elseif (ABS(delta) < 1.0D-8) Then   !1 intersection node.--tangent
        n_root = 1
        root1(1) = x1
        root1(2) = -deltaB / (2*deltaA)
        root2(1) = x1
        root2(2) = -deltaB / (2*deltaA)
      End If
    End If
    
    If (x1/=x2 .AND. ABS(y1-y2)<1.0D-8) Then    !horizontal edge.
      deltaA = thetaA
      deltaB = thetaB*y1 - thetaB*y_c - 2*thetaA*x_c
      deltaC = thetaC*y1*y1 - (2*thetaC*y_c+thetaB*x_c)*y1 + (thetaA*x_c*x_c+thetaB*x_c*y_c+thetaC*y_c*y_c-1)
      
      delta = deltaB**2 - 4*deltaA*deltaC
      If (delta < -1.0D-8) Then   !No intersection node.
        n_root = 0
      Elseif (delta > 1.0D-8) Then  !2 intersection nodes.
        n_root = 2
        root1(1) = (-deltaB+DSQRT(delta)) / (2*deltaA)
        root1(2) = y1
        root2(1) = (-deltaB-DSQRT(delta)) / (2*deltaA)
        root2(2) = y1
      Elseif (ABS(delta) < 1.0D-8) Then   !1 intersection node.--tangent
        n_root = 1
        root1(1) = -deltaB / (2*deltaA)
        root1(2) = y1
        root2(1) = -deltaB / (2*deltaA)
        root2(2) = y1
      End If
    End If
  End If
  !=========LY modification, 2022-6-1=========
  !=========Part 1: Calculate intersection node=========
  
  !=========Part 2: Judge relationship between intersection node and edge endpoint=========
  If (n_root/=0) Then   !At least 1 intersection node.
    rx1 = x1-root1(1)
	  ry1 = y1-root1(2)
	  rx2 = x2-root1(1)
	  ry2 = y2-root1(2)
    r = rx1*rx2+ry1*ry2
    
    If (r < -1.0D-8) Then
      root1_flag = 1
    Elseif (ABS(r) < 1.0D-8) Then
      If ((rx1*rx1+ry1*ry1) > (rx2*rx2+ry2*ry2)) Then
        root1_flag = -2
      Elseif ((rx1*rx1+ry1*ry1) < (rx2*rx2+ry2*ry2)) Then
        root1_flag = -1
      Else
        Print*,'Check Line_Ellipse_Intersection_BL, STOP.'
        Stop
      End If
    Elseif (r > 1.0D-8) Then
      root1_flag = 0
    End If
    
    rx1=x1-root2(1)
	  ry1=y1-root2(2)
	  rx2=x2-root2(1)
	  ry2=y2-root2(2)
	  r=rx1*rx2+ry1*ry2
    
    If (r < -1.0D-8) Then
      root2_flag = 1
    Elseif (ABS(r) < 1.0D-8) Then
      If ((rx1*rx1+ry1*ry1) > (rx2*rx2+ry2*ry2)) Then
        root2_flag = -2
      Elseif ((rx1*rx1+ry1*ry1) < (rx2*rx2+ry2*ry2)) Then
        root2_flag = -1
      Else
        Print*,'Check Line_Ellipse_Intersection_BL, STOP.'
        Stop
      End If
    Elseif (r > 1.0D-8) Then
      root2_flag = 0
    End If
  
  Elseif (n_root==0) Then   !No intersection node.
    root1_flag = 0
	  root2_flag = 0
  End If
  
  If (root1_flag/=0 .AND. root2_flag==0) Then
    inters_flag 	= root1_flag
    inters_point 	= root1
  Elseif (root1_flag==-1 .AND. root2_flag==-1) Then
    inters_flag 	= root1_flag
    inters_point 	= root1
  Elseif (root1_flag==-2 .AND. root2_flag==-2) Then
    inters_flag 	= root1_flag
    inters_point 	= root1
  Elseif (root1_flag==0 .AND. root2_flag/=0) Then
    inters_flag 	= root2_flag
    inters_point 	= root2
  Elseif ((root1_flag==-1 .AND. root2_flag==-2) .OR. (root1_flag==-2 .AND. root2_flag==-1) .OR. &
          (root1_flag==0 .AND. root2_flag==0)) Then
    inters_flag 	= 0
  Else
    Print*,'n_root=',n_root
    Print*,'root1=',root1, ' root1_flag=',root1_flag
    Print*,'root2=',root2, ' root2_flag=',root2_flag
    Print*,'LineEnd1=',x1,y1
    Print*,'LineEnd2=',x2,y2
    Print*,'Xc=',x_c,' Yc=',y_c,' R=',rx1,ry1,rx2,ry2
    Print*,'Check the Mesh, STOP'
    Stop
  End If
  !=========Part 2: Judge relationship between intersection node and edge endpoint=========
  
End Subroutine Line_Ellipse_Intersection_BL