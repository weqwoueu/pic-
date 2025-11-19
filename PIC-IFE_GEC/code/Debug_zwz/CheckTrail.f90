SUBROUTINE CheckTrail(N_Objects,objects)
    USE Constant_Variable_2D
    USE Field_2D
    USE Particle_2D
    USE TimeControl
    USE Object_Data_2D
    USE IFE_INTERFACE, ONLY: AdjustBoundary_2D
    USE PIC_MAIN_PARAM_2D !$ ab.ZWZ
    IMPLICIT NONE
    
    TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
    INTEGER    :: N_Objects
    
    INTEGER :: delta = 1
    
    REAL(8) :: xana, yana, vxana, vyana
    REAL(8)	:: xt
    REAL(8) :: En, Bn
    REAL(8) :: omigc
    INTEGER :: ispe, it, i
    REAL(8) :: vth, vr
    
    ntot = ntot + 1
    !part(1,:) = 0.0
    part(1,2) = 1.0
    part(1,6) = 1.0
    part(1,7) = 1
    ispe  =  INT(part(1,7))
    ns(ispe) = ns(ispe) + 1
    
    
    !ns(ISPI) = 1
    !part_i(1,1) = 0.0
    !part_i(1,2) = 0.0
    !part_i(1,4) = 0.0
    !part_i(1,5) = 0.0
    !part_i(1,6) = 0.0 
    !part_i(1,7) = 2
    
    efx = 0.
    !efy = -1
    efy = 1
    efz = 0.
    
    bfx = 0.
    !bfx = 1.
    bfy = 0.
    bfz = 0.
    !bfz = -1
    !En = 1
    !Bn = 1.
    En = 1
    Bn = 1
    omigc = qs(ispe)*Bn/xm(ispe)
    
    xt = 0.
    dt = 0.01
    
    !delta = 0
    xana = 0
    yana = 0
    vxana = 0
    vyana = 0
    OPEN(1, ACTION = 'WRITE', FILE ='./Trail_e.dat')	      
        WRITE(1,*) 'TITLE = "Trail Plot"'
        WRITE(1,*) 'VARIABLES = "z" "r" "ze" "re" "vx" "vy" "vxe" "vye" "xt"'
        WRITE(1,"(9(E15.6))")xana ,yana , part(1,1:2), vxana, vyana, part(1,4:5), xt
    CLOSE(1)
    
    OPEN(11, ACTION = 'WRITE', FILE ='./Trail_e_circular.dat')	      
        WRITE(11,*) 'TITLE = "Trail Plot"'
        WRITE(11,*) 'VARIABLES = "x" "y" "xana" "yana" "xt"'
        WRITE(11,"(5(E15.6))")part(1,2)*DSIN(part(1,3)),part(1,2)*DCOS(part(1,3)),yana*DSIN(part(1,3)), yana*DSIN(part(1,3)),xt
    CLOSE(11)
    !OPEN(2, ACTION = 'WRITE', FILE ='./OUTPUT/Trail_i.dat')	      
    !    WRITE(2,*) 'TITLE = "Trail Plot"'
    !    WRITE(2,*) 'VARIABLES = "x" "y" "xi" "yi" "xt"'
    !CLOSE(2)
    DO it = ilap+1, nt 
        xt=xt + dt
                                 
        CALL Move_2D(it,dt,delta)
                
        !CALL PrePush
        !CALL PostPush
                
        !CALL AdjustBoundary_2D(xt,dt,delta,N_Objects,objects)	       

            
        omigc = qs(ispe)*Bn/xm(ispe)
                
        xana = En/Bn * xt - xm(ispe)*En/qs(ispe)/Bn**2*sin(OMIGC*xt)
        vxana = En/Bn - xm(ispe)*En/qs(ispe)/Bn**2 * OMIGC*cos(OMIGC*xt)
        yana = xm(ispe)*(-En)/qs(ispe)/Bn**2*(1-cos(OMIGC*xt))
        vyana = xm(ispe)*(-En)/qs(ispe)/Bn**2 * OMIGC*sin(OMIGC*xt)
                
        OPEN(1, ACTION = 'WRITE', FILE ='./Trail_e.dat',POSITION='APPEND')	
            DO i = 1, ns(ispe)
                WRITE(1,"(9(E15.6))")xana ,yana , part(i,1:2), vxana, vyana, part(i,4:5), xt
            END DO
        CLOSE(1)
        
        OPEN(11, ACTION = 'WRITE', FILE ='./Trail_e_circular.dat',POSITION='APPEND')	      
            DO i = 1, ns(ispe)
                WRITE(11,"(5(E15.6))")part(1,2)*DSIN(part(i,3)),part(1,2)*DCOS(part(1,3)),SIN(pi*xt), COS(pi*xt), xt
            END DO
        CLOSE(11)
        
        OPEN(22, ACTION = 'WRITE', FILE ='./sincos.dat',POSITION='APPEND')	      
            DO i = 1, ns(ispe)
                WRITE(22,"(5(E15.6))")xt, sin(xt), cos(xt), dsin(xt), dcos(xt)
            END DO
        CLOSE(22)

        print*,it
        !if(ntot == 0)THEN
        if(part(1,1) > 70. .or. ntot == 0)THEN
            print*,part(1,:)
            print*,ntot
            pause
            exit
        ENDIF
        if (it == 1700) pause
    ENDDO
END SUBROUTINE