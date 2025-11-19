Module ModuleParticleBoundary
  Use Constants
  use Numrical
  Use ModuleParticleBundle
  Use ModuleSEE !> ab.ZWZ
  Implicit none
      Type ParticleBoundaryOne
             !Integer(4) :: XStart=0,XEnd=1
             Integer(4) :: ParticleBoundaryModel=15
             Real(8) :: XMin=0.d0,XMax=dble(NxMax-1)!,XLength=99.d0 
             Integer :: CountMinOne=0,CountmaxOne=0
				 Integer :: SEECountMinOne=0,SEECountmaxOne=0
				 Integer :: FEECountMinOne=0,FEECountmaxOne=0
				 Integer :: TEECountMinOne=0,TEECountmaxOne=0
             Integer :: CountMin=0,Countmax=0
             Real(8) :: Gamma=0.01     !ion Secondary Electron emission rate
             !Integer :: SEENparMin,SEENparMax
             Type(ParticleBundle) :: PBLower,PBUpper
             Contains
             Procedure :: AllInit=>InitializationParticleBoundaryOne
		End Type ParticleBoundaryOne
								 
						 
       
    contains
    subroutine InitializationParticleBoundaryOne(PBDO,PB,CF)
        implicit none
        Class(ParticleBoundaryOne),intent(inout) :: PBDO
        Type(ParticleBundle),intent(in) :: PB
        Class(ControlFlow), intent(in) :: CF
		  Type(ParticleBundle) :: TempPB
		  Integer(4) :: i,NPArMax
		  NPArMax = Ceiling(2.d0*PB%NParNormal)
		  TempPB = PB
			If(Allocated( TempPB%PO)) Deallocate(TempPB%PO)
         Allocate(TempPB%PO(NPArMax))
			 Do i=1,NPArMax
             TempPB%PO(i)%X  = 0.d0
				 TempPB%PO(i)%Vx = 0.d0
				 TempPB%PO(i)%Vy = 0.d0
             TempPB%PO(i)%Vz = 0.d0
             TempPB%PO(i)%Ax = 0.d0
             TempPB%PO(i)%Ay = 0.d0
				 TempPB%PO(i)%Az = 0.d0
         End Do
             PBDO%PBLower= TempPB
             PBDO%PBUpper= TempPB
             PBDO%PBLower%NPar=0
             PBDO%PBUpper%NPar=0
             PBDO%XMin=Dble(CF%NxL)
             PBDO%XMax=Dble(CF%NxU)
        return
  end subroutine InitializationParticleBoundaryOne
    
    
   subroutine ParticleAborption(PB,PBDO)
        implicit none
        Type(ParticleBundle),intent(inout) :: PB
        Type(ParticleBoundaryOne),intent(inout) :: PBDO
        Integer :: i
		  PBDO%PBLower%Weight = PB%Weight
		  PBDO%PBUpper%Weight = PB%Weight
        Select case (PBDO%ParticleBoundaryModel)
            case(10)
                do i=PB%NPar,1,-1
                     If (PB%PO(i)%X<=PBDO%XMin.or.PB%PO(i)%X>=PBDO%XMax) then
                        !Write(*,*) "Xmin"
                         Call PB%DelOne(i)  
                     end if
                end do
            case(11)
                PBDO%CountMinOne=0
                PBDO%CountMaxOne=0
                do i=PB%NPar,1,-1
                     If (PB%PO(i)%X<=PBDO%XMin) then
                         PBDO%CountMinOne=PBDO%CountMinOne+1
                         !Write(*,*) "Xmin"
                         Call PB%DelOne(i)  
                     else If(PB%PO(i)%X>=PBDO%XMax) then
                         PBDO%CountMaxOne=PBDO%CountMaxOne+1
                         !Write(*,*) "Xmax"
                        Call PB%DelOne(i)                       
                     end if
                end do
                PBDO%CountMin=PBDO%CountMin+PBDO%CountMinOne
                PBDO%CountMax=PBDO%CountMax+PBDO%CountMaxOne
            case(12)
                PBDO%CountMinOne=0
                PBDO%CountMaxOne=0
                PBDO%PBLower%NPar=0
                PBDO%PBUpper%NPar=0
                do i=PB%NPar,1,-1
                     If (PB%PO(i)%X<=PBDO%XMin) then
                         PBDO%CountMinOne=PBDO%CountMinOne+1
                         Call PBDO%PBLower%AddOne(PB%PO(i))
                         Call PB%DelOne(i)  
                     else If(PB%PO(i)%X>=PBDO%XMax) then
                         PBDO%CountMaxOne=PBDO%CountMaxOne+1
                         Call PBDO%PBUpper%AddOne(PB%PO(i))
                         Call PB%DelOne(i)                          
                     end if
					 end do
                PBDO%CountMin=PBDO%CountMin+PBDO%CountMinOne
                PBDO%CountMax=PBDO%CountMax+PBDO%CountMaxOne
				case(14)
                PBDO%CountMinOne=0
                PBDO%CountMaxOne=0
                PBDO%PBLower%NPar=0
                PBDO%PBUpper%NPar=0
			    	 PBDO%SEECountMaxOne = 0
				    PBDO%SEECountMinOne = 0
				    PBDO%FEECountMaxOne = 0
					 PBDO%FEECountMinOne = 0
					 PBDO%TEECountMaxOne = 0
					 PBDO%TEECountMinOne = 0
                do i=PB%NPar,1,-1
                     If (PB%PO(i)%X<=PBDO%XMin) then
                         PBDO%CountMinOne=PBDO%CountMinOne+1
                         Call PBDO%PBLower%AddOne(PB%PO(i))
                         Call PB%DelOne(i)  
                     else If(PB%PO(i)%X>=PBDO%XMax) then
                         PBDO%CountMaxOne=PBDO%CountMaxOne+1
                         Call PBDO%PBUpper%AddOne(PB%PO(i))
                         Call PB%DelOne(i)                          
                     end if
                end do
                PBDO%CountMin=PBDO%CountMin+PBDO%CountMinOne
                PBDO%CountMax=PBDO%CountMax+PBDO%CountMaxOne 
					 
					case(15)
                PBDO%CountMinOne=0
                PBDO%CountMaxOne=0
                PBDO%PBLower%NPar=0
                PBDO%PBUpper%NPar=0
					 PBDO%SEECountMaxOne = 0
					 PBDO%SEECountMinOne = 0
					 PBDO%FEECountMaxOne = 0
					 PBDO%FEECountMinOne = 0
					 PBDO%TEECountMaxOne = 0
					 PBDO%TEECountMinOne = 0
                do i=PB%NPar,1,-1
                     If (PB%PO(i)%X<=PBDO%XMin) then
                         PBDO%CountMinOne=PBDO%CountMinOne+1
                         Call PBDO%PBLower%AddOne(PB%PO(i))
                         Call PB%DelOne(i)  
                     else If(PB%PO(i)%X>=PBDO%XMax) then
                         PBDO%CountMaxOne=PBDO%CountMaxOne+1
                         Call PBDO%PBUpper%AddOne(PB%PO(i))
                         Call PB%DelOne(i)                          
                     end if
					 end do
                PBDO%CountMin=PBDO%CountMin+PBDO%CountMinOne
                PBDO%CountMax=PBDO%CountMax+PBDO%CountMaxOne 
        End Select 

        return
    end subroutine ParticleAborption
   
    
    
    Subroutine LineIntersection(LineEnds,xyp,InterPoint,InterFlag)
        Implicit none
        Real(8),intent(in) :: LineEnds(2,2),xyp(2,2)
        Real(8),intent(out) :: InterPoint(2)
        Integer(4),intent(out) :: InterFlag
        
        Real(8) :: x1,y1,x2,y2,xi,yi,xp,yp,k1,k2
        Real(8) :: Ferror = 1D-7
        
        x1 = LineEnds(1,1)
        y1 = LineEnds(1,2)
        x2 = LineEnds(2,1)
        y2 = LineEnds(2,2)
        
        
        InterFlag = 0
        InterPoint = (/0.0,0.0/)
        
        If (Abs(x1-x2) < Ferror .and. y1 /= y2) Then ! vertical line
            xi = x1
            k1 = (xyp(2,2) - xyp(1,2)) / (xyp(2,1) - xyp(1,1))
            yi = xyp(1,2) + k1 * ( xi - xyp(1,1) )
            If ( (y1-yi)*(y2-yi)<=0 ) Then ! intersection point lie between the line
                If ((xyp(1,2)-yi)*(xyp(2,2)-yi)<=0 .and. (xyp(1,1)-xi)*(xyp(2,1)-xi)<=0) Then ! Intersection point lie between the trail
                    If (Abs(y1-yi)<Ferror ) Then	    !  Intersection point at first END If			
				        InterFlag		= -1
			        Elseif (Abs(y2-yi)<Ferror ) Then	! Intersection point at second END If
				        InterFlag		= -2
			        Else
				        InterFlag		= 1
			        Endif
			        InterPoint	= (/xi,yi/)
                End If
            End If
            
        Elseif (x1 /= x2 .and. Abs(y1-y2) < Ferror) Then ! horizontal line
            yi=y1
	        k2=(xyp(2,1)-xyp(1,1))/(xyp(2,2)-xyp(1,2))
	        xi=(y1-xyp(1,2))*k2+xyp(1,1)
	        If ((x1-xi)*(x2-xi)<=0 ) Then
		        If ((xyp(1,1)-xi)*(xyp(2,1)-xi)<=0 .and. (xyp(1,2)-yi)*(xyp(2,2)-yi)<=0 ) Then  
			        If (Abs(x1-xi)<Ferror ) Then	
				        InterFlag		= -1
			        Elseif (Abs(x2-xi)<Ferror ) Then	
				        InterFlag		= -2
			        Else
				        InterFlag		= 1
			        Endif
			        InterPoint	= (/xi,yi/)
		        Endif
	        Endif
        Elseif (x1 /= x2 .and. y1 /= y2) Then ! oblique line
            If (xyp(2,1)/=xyp(1,1)) Then
		        k1=(y2-y1)/(x2-x1) 
		        k2=(xyp(2,2)-xyp(1,2))/(xyp(2,1)-xyp(1,1))
		        xi=(y1-xyp(1,2)+xyp(1,1)*k2-x1*k1)/(k2-k1)
		        yi=k1*(xi-x1) + y1
		        If ((x1-xi)*(x2-xi)<=Ferror .and. (y1-yi)*(y2-yi)<=Ferror ) Then
			        If ((xyp(2,1)-xi)*(xyp(1,1)-xi)<=Ferror .and. (xyp(1,2)-yi)*(xyp(2,2)-yi)<=Ferror ) Then
				        If (Abs(y1-yi)<Ferror .and. Abs(x1-xi)<Ferror ) Then	
					        InterFlag		= -1
				        Elseif (Abs(y2-yi)<Ferror .and. Abs(x2-xi)<Ferror ) Then	
					        InterFlag		= -2
				        Else
					        InterFlag		= 1
				        Endif
				        InterPoint	= (/xi,yi/)

			        Endif
		        Endif
	        Else
		        k1=(y2-y1)/(x2-x1)
		        xi=xyp(2,1)
		        yi=k1*(xi-x1) + y1
		        If ((x1-xi)*(x2-xi)<=Ferror .and. (y1-yi)*(y2-yi)<=Ferror ) Then
			        If ((xyp(1,2)-yi)*(xyp(2,2)-yi)<=0 ) Then		
				       If (Abs(y1-yi)<Ferror .and. Abs(x1-xi)<Ferror ) Then				
					        InterFlag		= -1
				        Elseif (Abs(y2-yi)<Ferror .and. Abs(x2-xi)<Ferror ) Then
					        InterFlag		= -2
				        Else
					        InterFlag		= 1
				        Endif
				        InterPoint	= (/xi,yi/)
			        Endif
		        Endif
	        Endif
        Else    ! a point
            xi = x1
            yi = y1
			If ((xyp(2,1)-xi)*(xyp(1,1)-xi)<=Ferror .and. (xyp(1,2)-yi)*(xyp(2,2)-yi)<=Ferror ) Then		
                InterFlag   = -1
				InterPoint	= (/xi,yi/)
			Endif
        End If
    End Subroutine LineIntersection
   
    Subroutine ZCircularPlaneIntersection_3D(Z,R,Line,InterPoint,InterFlag)
        Implicit none
        Real(8),intent(in) :: Z,R
        Real(8),intent(in) :: Line(2,3)
        Real(8),intent(out) :: InterPoint(3)
        Integer(4),intent(out) :: InterFlag
        
        Real(8) :: x0,y0,z0,x1,x2,y1,y2,z1,z2
        
        x1 = Line(1,1)
        y1 = Line(1,2)
        z1 = Line(1,3)
        x2 = Line(2,1)
        y2 = Line(2,2)
        z2 = Line(2,3)
        
        InterFlag = 0
        InterPoint = (/0.,0.,0./)
        If ((Z1-Z)*(Z2-Z)<=0) Then ! line cross the plane
            z0 = Z
            x0 = (x2-x1) * (Z-z1) / (z2-z1) + x1
            y0 = (y2-y1) * (Z-z1) / (z2-z1) + y1
            If ( (x0**2 + y0**2) <= R**2) Then ! intersect in the circle
                InterFlag = 1
                InterPoint = (/x0,y0,z0/)
            End If
        End If
   
    End Subroutine ZCircularPlaneIntersection_3D
    
    ! surface: x^2 + y^2 = R^2
    ! Line: x = t(x2 - x1) + x1
    !       y = t(y2 - y1) + y1
    !       z = t(z2 - z1) + z1
    Subroutine ZCylinderIntersection_3D(R,ZL,ZU,Line,InterPoint,InterFlag)
        Real(8),intent(in) :: R,ZL,ZU
        Real(8),intent(in) :: Line(2,3)
        Real(8),intent(out) :: InterPoint(3)
        Integer(4),intent(out) :: InterFlag
        
        Real(8) :: x0,y0,z0,x1,y1,z1,x2,y2,z2,t1,t2
        Real(8) :: a,b,c,kxz,kyz,kyx,delta
        
        x1 = Line(1,1)
        y1 = Line(1,2)
        z1 = Line(1,3)
        x2 = Line(2,1)
        y2 = Line(2,2)
        z2 = Line(2,3)
        
        InterFlag = 0
        InterPoint = (/0.,0.,0./)
        
        a = (x2-x1)**2 + (y2-y1)**2
        b = 2*x1*(x2-x1)+2*y1*(y2-y1)**2
        c = x1**2+y1**2-R**2
        delta = b**2 - 4*a*c
        If (delta >= 0) Then            
            t1 = (-b+Sqrt(delta))/(2*a)
            x0 = t1*(x2-x1) + x1
            y0 = t1*(y2-y1) + y1
            z0 = t1*(z2-z1) + z1
            If ( (ZL-z0)*(ZU-z0)<=0 ) Then ! lie in the cylinder z section
                If ((x1-x0)*(x2-x0)<=0 .or. (y1-y0)*(y2-y0)<=0 .or. (z1-z0)*(z2-z0)<=0) Then ! lie in the line 
                    InterFlag = 1
                    InterPoint = (/x0,y0,z0/)
                End If
            End If
            t2 = (-b-Sqrt(delta))/(2*a)
            x0 = t2*(x2-x1) + x1
            y0 = t2*(y2-y1) + y1
            z0 = t2*(z2-z1) + z1
            If ( (ZL-z0)*(ZU-z0)<=0) Then ! lie in the cylinder z section
                If ((x1-x0)*(x2-x0)<=0 .or. (y1-y0)*(y2-y0)<=0 .or. (z1-z0)*(z2-z0)<=0) Then ! lie in the line 
                    InterFlag = 1
                    InterPoint = (/x0,y0,z0/)
                End If
            End If
        End If
    End Subroutine ZCylinderIntersection_3D
   

   
End  Module ModuleParticleBoundary