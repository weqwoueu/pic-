SUBROUTINE Setup_IFE_Wall_Mesh_2D_NEW(objects)

! Purpose:		Setup IFE Mesh, do objects-mesh intersections and save mesh information to a IFE_Mesh_Filename.
! Last Updated:	4/21/2004 01:21 PM

USE PIC_MAIN_PARAM_2D
USE Wall_2D
USE Particle_2D
USE Domain_2D
USE Object_Data_2D
USE Boundary_Data_2D
USE TimeControl
USE IFE_Interface, ONLY: Curve_tracing_2D, Get_Wave_Direction_2D

IMPLICIT NONE

!INCLUDE 'Curve_tracing_2D.inc'
!INCLUDE 'Delete_Used_Elements_2D.inc'
!INCLUDE 'Get_Wave_Direction_2D.inc'

INTEGER			N_Objects

INTEGER		 num_of_nodes, n_int_elements, i, j,num_of_waindex,ii
REAL(8)		dimensions(2,2),xo,yo,xc,yc,rc,rp,res
INTEGER		count,nnx,nny,n_sect,N_Boundary,nto,nynode,i_obj

INTEGER		nodes(2)

INTEGER, DIMENSION(:,:), POINTER			::	t_c
REAL(8), DIMENSION(:,:), POINTER			::	p_basic
REAL(8), DIMENSION(:), POINTER			::	xp_tmp,yp_tmp
REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp2,p_int_y_tmp2
REAL(8), DIMENSION(:,:), POINTER		::	p_int_x_tmp1,p_int_y_tmp1
REAL(8), DIMENSION(:,:), POINTER			::	information_2
REAL(8), DIMENSION(:,:,:), POINTER		::	wallindex_tmp
REAL(8), DIMENSION(:), ALLOCATABLE		::	ynode
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects

LOGICAL		ATEND

WRITE(6,*) 'Setup_IFE_Wall_Mesh_2D '

N_Boundary = SIZE(boundaries,1)

print *, 'N_Boundary= ',N_Boundary

!!num_of_waindex	= 5*SIZE(p_basic,2)
!!num_of_nodes= SIZE(p_basic,2)
!!n_int_elements= SIZE(information_2,2)
!!ALLOCATE(p_int_x_tmp1(2,n_int_elements),p_int_y_tmp1(2,n_int_elements))
!!ALLOCATE(wallindex_tmp(num_of_waindex,2,N_Boundary))
!!
!!p_int_x_tmp1(1,1:n_int_elements)=information_2(3,1:n_int_elements)
!!p_int_x_tmp1(2,1:n_int_elements)=information_2(5,1:n_int_elements)
!!p_int_y_tmp1(1,1:n_int_elements)=information_2(4,1:n_int_elements)
!!p_int_y_tmp1(2,1:n_int_elements)=information_2(6,1:n_int_elements)
!!wallindex_tmp(:,:,:) = Zero
!!
!!ALLOCATE(nnode(N_Boundary))
!!nnode(:) = Zero

N_Objects	=SIZE(objects,1)
print *, 'N_Objects= ',N_Objects
!IF(irestart==1 .AND. ilap> Erostart )THEN
!
!print *, 'irestart= ',irestart
!
!    OPEN(1, ACTION='READ', FILE='update_wall.msh')
!    DO j = 1, N_Objects
!        IF (objects(j)%Erosion > 0) THEN
!            READ(1,*) nnode(j)
!	    ENDIF
!	    objects(j)%Wall(1)%nxp = nnode(j)
!    END DO
!    DO j=1,N_Objects
!	    ALLOCATE(objects(j)%Wall(1)%node(MAXVAL(nnode),2))
!	    IF (objects(j)%Erosion > 0) THEN
!	        DO i=1,nnode(j)
!		        READ(1,*) objects(j)%Wall(1)%Node(i,1), objects(j)%Wall(1)%Node(i,2)
!	        ENDDO
!	    ENDIF
!    ENDDO 
!    CLOSE(1)
!
!ELSE

print *, 'unrestart '
!!
!!DO j=1,N_Objects
!!
!!    print *,'j,objects(j)%Erosion= ',j,objects(j)%Erosion
!!    IF (objects(j)%Erosion > 0) THEN
!!        objects(j)%Wall(1)%nxp = Zero
!!        
!!        IF	(objects(j)%Shape==3) THEN	!box
!!!			nnx=(objects(j)%Wall(1)%Limits(2,1)-objects(j)%Wall(1)%Limits(1,1))/objects(j)%Wall(1)%stepx+1
!!!			nny=(objects(j)%Wall(1)%Limits(2,2)-objects(j)%Wall(1)%Limits(1,2))/objects(j)%Wall(1)%stepx
!!!			res=objects(j)%Wall(1)%Limits(2,2)-objects(j)%Wall(1)%Limits(1,2)-nny*objects(j)%Wall(1)%stepx
!!!			xo=objects(j)%Wall(1)%Limits(1,1)
!!!			IF (objects(j)%Wall(1)%Channelwall==3) THEN
!!!			    nto=nnx+nny
!!!				nnode(j)=nto
!!!				yo=objects(j)%Wall(1)%Limits(1,2)
!!!				IF(res>0)THEN
!!!				    nnode(j)=nnode(j)+1
!!!				    wallindex_tmp(nnode(j),1:2,j)=objects(j)%Wall(1)%Limits(2,1:2)
!!!				ENDIF
!!!			ELSEIF (objects(j)%Wall(1)%Channelwall==4) THEN
!!!				nto=nnx+nny
!!!				nnode(j)=nto
!!!				yo=objects(j)%Wall(1)%Limits(2,2)
!!!				IF(res>0)THEN
!!!				    nnode(j)=nnode(j)+1
!!!				    wallindex_tmp(nnode(j),2,j)=objects(j)%Wall(1)%Limits(1,2)
!!!				    wallindex_tmp(nnode(j),1,j)=objects(j)%Wall(1)%Limits(2,1)
!!!				ENDIF				
!!!			ELSE
!!!				nnode(j)=2*(nnx+nny-1)
!!!				yo=objects(j)%Wall(1)%Limits(1,2)
!!!			ENDIF
!!!			DO i=1,nto
!!!			IF (objects(j)%Wall(1)%Channelwall==3) THEN
!!!				IF (i<=nnx) THEN
!!!					wallindex_tmp(i,1,j)=xo+(objects(j)%Wall(1)%stepx*(i-1))
!!!					wallindex_tmp(i,2,j)=yo
!!!				ELSE
!!!					wallindex_tmp(i,1,j)=wallindex_tmp(i-1,1,j)
!!!					wallindex_tmp(i,2,j)=yo+(objects(j)%Wall(1)%stepx*(i-1-nnx))
!!!				ENDIF
!!!			ELSEIF (objects(j)%Wall(1)%Channelwall==4) THEN
!!!				IF (i<=nnx) THEN
!!!					wallindex_tmp(i,1,j)=xo+(objects(j)%Wall(1)%stepx*(i-1))
!!!					wallindex_tmp(i,2,j)=yo
!!!				ELSE
!!!					wallindex_tmp(i,1,j)=wallindex_tmp(i-1,1,j)
!!!					wallindex_tmp(i,2,j)=yo-(objects(j)%Wall(1)%stepx*(i-nnx))
!!!				ENDIF
!!!			ELSEIF (objects(j)%Wall(1)%Channelwall==0) THEN
!!!				IF (i<=nny+1) THEN
!!!					wallindex_tmp(i,1,j)=xo
!!!					wallindex_tmp(i,2,j)=yo+(objects(j)%Wall(1)%stepx*(i-1))
!!!				ELSEIF (i<=nnx+nny .and. i>nny+1) THEN
!!!					wallindex_tmp(i,1,j)=xo+(objects(j)%Wall(1)%stepx*(i-1-nny))
!!!					wallindex_tmp(i,2,j)=wallindex_tmp(i-1,2,j)
!!!				ELSEIF (i<=nnx+2*nny .and. i>nnx+nny) THEN
!!!					wallindex_tmp(i,1,j)=wallindex_tmp(i-1,1,j)
!!!					wallindex_tmp(i,2,j)=wallindex_tmp(nnx+nny,2,j)-(objects(j)%Wall(1)%stepx*(i-nnx-nny))
!!!				ELSE
!!!					wallindex_tmp(i,1,j)=wallindex_tmp(nnx+2*nny,1,j)-(objects(j)%Wall(1)%stepx*(i-nnx-2*nny))
!!!					wallindex_tmp(i,2,j)=yo
!!!				ENDIF
!!!			ENDIF
!!!			ENDDO
!!        nnode(j)=n_int_elements+1
!!        print *,nnode(j),j
!!
!!		ALLOCATE(xp_tmp(num_of_nodes),yp_tmp(num_of_nodes))
!!		ALLOCATE(p_int_x_tmp2(2,n_int_elements),p_int_y_tmp2(2,n_int_elements))
!!		xp_tmp=-200
!!		yp_tmp=-200
!!		p_int_x_tmp2(1:2,1:n_int_elements)=p_int_x_tmp1(1:2,1:n_int_elements)
!!		p_int_y_tmp2(1:2,1:n_int_elements)=p_int_y_tmp1(1:2,1:n_int_elements)
!!
!!		xp_tmp(1)=objects(j)%Locations(1,1)
!!		IF(j==1)THEN
!!		    yp_tmp(1)=objects(j)%Locations(2,2)
!!		ELSEIF(j==2)THEN
!!		    yp_tmp(1)=objects(j)%Locations(1,2)
!!		ENDIF
!!
!!		xc=0        !objects(j)%Locations(1,1)
!!		yc=0        !objects(j)%Locations(1,2)
!!
!!		CALL Get_Wave_Direction_2D(n_int_elements, p_int_x_tmp2, 	&
!!						p_int_y_tmp2, xp_tmp, yp_tmp , xc, yc)
!!
!!		count=2
!!141		CONTINUE
!!
!!		CALL Curve_tracing_2D(n_int_elements, xp_tmp(count-1:count+1), yp_tmp(count-1:count+1), p_int_x_tmp2,p_int_y_tmp2)
!!		count=count+1
!!
!!!			IF(xp_tmp(count) == objects(j)%Locations(2,1)) THEN
!!!				GO TO 241
!!!			ENDIF
!!
!!140		IF(xp_tmp(count) /= objects(j)%Wall(1)%Limits(2,1)) GO TO 141
!!
!!		nnode(j)=count
!!!			xp_tmp(count)=objects(j)%Wall(1)%node(n_xo,1)
!!!			yp_tmp(count)=objects(j)%Wall(1)%node(n_xo,2)
!!		p_int_x_tmp1(1:2,1:n_int_elements)=p_int_x_tmp2(1:2,1:n_int_elements)
!!		p_int_y_tmp1(1:2,1:n_int_elements)=p_int_y_tmp2(1:2,1:n_int_elements)
!!!		print *,nnode(j),j
!!		DO i=1,count		!nnode(j)
!!			wallindex_tmp(i,1,j)=xp_tmp(i)
!!			wallindex_tmp(i,2,j)=yp_tmp(i)
!!		ENDDO
!!				
!!		nnode(j)=nnode(j)+2
!!		IF(j==1)THEN
!!		wallindex_tmp(nnode(j)-1,1,j)=wallindex_tmp(nnode(j)-2,1,j)     !objects(j)%Locations(2,1)
!!		wallindex_tmp(nnode(j)-1,2,j)=objects(j)%Locations(1,2)
!!		wallindex_tmp(nnode(j),1,j)=objects(j)%Locations(1,1)
!!		wallindex_tmp(nnode(j),2,j)=objects(j)%Locations(1,2)
!!		ELSEIF(j==2)THEN
!!		wallindex_tmp(nnode(j)-1,1,j)=wallindex_tmp(nnode(j)-2,1,j)     !objects(j)%Locations(2,1)
!!		wallindex_tmp(nnode(j)-1,2,j)=objects(j)%Wall(1)%Limits(2,2)
!!		wallindex_tmp(nnode(j),1,j)=objects(j)%Locations(1,1)
!!		wallindex_tmp(nnode(j),2,j)=objects(j)%Wall(1)%Limits(2,2)
!!		ENDIF
!!	ENDIF
!!objects(j)%Wall(1)%nxp = nnode(j)
!!ENDIF
!!ENDDO			
!!			
!!DO j=1,N_Objects
!!	ALLOCATE(objects(j)%Wall(1)%node(MAXVAL(nnode),2))
!!	IF (objects(j)%Erosion > 0) THEN
!!	    DO i=1,nnode(j)
!!		    objects(j)%Wall(1)%node(i,1) = wallindex_tmp(i,1,j)
!!		    objects(j)%Wall(1)%node(i,2) = wallindex_tmp(i,2,j)
!!	    ENDDO
!!	ENDIF
!!ENDDO
!!
!!ENDIF

DO j=1,N_Objects
    objects(j)%Wall(1)%nxp = 103
	ALLOCATE(objects(j)%Wall(1)%node(103,2))
	IF (objects(j)%Erosion > 0) THEN
	    DO i=1,101
		    objects(j)%Wall(1)%node(i,1) = (i-1)*hx(1)
		    objects(j)%Wall(1)%node(i,2) = 0
	    ENDDO
	    objects(j)%Wall(1)%node(102,1) = 100.0*hx(1)
	    objects(j)%Wall(1)%node(103,1) = 0.0*hx(1)
		objects(j)%Wall(1)%node(102:103,2) = 0
	ENDIF
ENDDO

!IF(irestart==1 .AND. ilap>= Erostart )THEN

print *, 'irestart= ',irestart

!    OPEN(1, ACTION='READ', FILE='update_wall.msh')
!    DO j = 1, 101
!        READ(1,*) xo,yo 
!        objects(1)%Wall(1)%Node(j,2)=yo/2.35E-2
!!        PRINT *,'yo=',yo
!    END DO
!    CLOSE(1)
!ENDIF

!!!!!!!!!!!!!!!!!chj add at 20141215!!!!!!!!!!!!!!!!!!!!
!print*,'hx(1)',hx(1),an_gridt(1)
!DO j = 1, 101
!   objects(1)%Wall(1)%Node(j,2)= 19-0.6*DSIN(TWO_PI*j*hx(1)*4/100-HALF_PI)   !!!!!! 振幅不能为1
!END DO

OPEN(1,ACTION='READ',FILE='Wall13.dat')
    READ(1,*) 
    READ(1,*) 
    DO i=1,101          
    	READ(1,*) xo,yo
    	objects(1)%Wall(1)%Node(i,2)=yo
    END DO

CLOSE(1)
!!!!!!!!!!!!!!!!!chj add at 20141215!!!!!!!!!!!!!!!!!!!!


boundaries(3)%Locations(1,2)=objects(1)%Wall(1)%Node(1,2)
boundaries(4)%Locations(1,2)=objects(1)%Wall(1)%Node(101,2)
print *, 'boundaries(3)%Locations(1,2)= ',boundaries(3)%Locations(1,2),boundaries(4)%Locations(1,2)

OPEN(1, ACTION='WRITE',FILE='Wall_morphologyobj1.dat')
WRITE(1,*) 'TITLE = "Field Plot"'
WRITE(1,*) 'VARIABLES = "x" "y"'
DO j=1,N_Objects
IF (objects(j)%Erosion > 0) THEN
	DO i=1,objects(j)%Wall(1)%nxp
		WRITE(1,51) objects(j)%Wall(1)%Node(i,1), objects(j)%Wall(1)%Node(i,2)
	ENDDO
ENDIF
ENDDO
51 FORMAT (E15.6,' ',E15.6)
CLOSE(1)

!nynode=MAXVAL(nnode)
!ALLOCATE(ynode(nynode))
!ynode=0
!nynode=nynode+1
nynode=103

N_Boundary = SIZE(boundaries,1)

DO j=1,N_Boundary
    boundaries(j)%Wall(1)%nxp = Zero

    IF (boundaries(j)%Erosion > 0) THEN
        ALLOCATE(boundaries(j)%Wall(1)%node(nynode,2))

        i_obj=boundaries(j)%Obound
        IF	(boundaries(j)%Shape==1) THEN	!box

			IF (boundaries(j)%Wall(1)%Channelwall==3 .or. boundaries(j)%Wall(1)%Channelwall==4) THEN

				DO i=1,nynode       !nnode(i_obj)
				        
			        boundaries(j)%Wall(1)%nxp = boundaries(j)%Wall(1)%nxp + 1
				    boundaries(j)%Wall(1)%node(i,1:2)=objects(i_obj)%Wall(1)%node(i,1:2)
!					IF (ABS(objects(i_obj)%Wall(1)%node(i,1)-objects(i_obj)%Locations(2,1))<GTol) THEN
!                        exit
!				    ENDIF
				ENDDO
!				boundaries(j)%Wall(1)%nxp = boundaries(j)%Wall(1)%nxp + 3
!                boundaries(j)%Wall(1)%node(1,2)=boundaries(j)%Wall(1)%node(2,2)
!				IF(i_obj==1)THEN
!					boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp-2,1:2)=boundaries(j)%Locations(2,:)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp-1,1)=boundaries(j)%Locations(2,1)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp-1,2)=objects(i_obj)%Locations(1,2)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp,1)=objects(i_obj)%Locations(1,1)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp,2)=objects(i_obj)%Locations(1,2)
!				ELSEIF(i_obj==2)THEN
!					boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp-2,1:2)=boundaries(j)%Locations(1,:)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp-1,1)=boundaries(j)%Locations(1,1)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp-1,2)=objects(i_obj)%Locations(2,2)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp,1)=objects(i_obj)%Locations(1,1)
!				    boundaries(j)%Wall(1)%node(boundaries(j)%Wall(1)%nxp,2)=objects(i_obj)%Locations(2,2)
!				ENDIF
			ELSEIF (boundaries(j)%Wall(1)%Channelwall==1 .or. boundaries(j)%Wall(1)%Channelwall==2) THEN

                boundaries(j)%Wall(1)%nxp = 2
                boundaries(j)%Wall(1)%node(1:2,1) = boundaries(j)%Locations(1,1)

!				DO i=1,nnode(i_obj)
!				    IF (ABS(objects(i_obj)%Wall(1)%node(i,1)-objects(i_obj)%Locations(2,1))<GTol) THEN
!				        nynode=nynode+1
!					    ynode(nynode)=objects(i_obj)%Wall(1)%node(i,2)
!				    ENDIF
!				ENDDO
				
				IF (boundaries(j)%Wall(1)%Channelwall==1) THEN	
			        boundaries(j)%Wall(1)%node(1,2) = boundaries(j)%Locations(1,2)
				    boundaries(j)%Wall(1)%node(2,2)=boundaries(j)%Locations(2,2)
				ELSEIF(boundaries(j)%Wall(1)%Channelwall==2)THEN
			        boundaries(j)%Wall(1)%node(1,2) = boundaries(j)%Locations(2,2)
				    boundaries(j)%Wall(1)%node(2,2)=boundaries(j)%Locations(1,2)				
				ENDIF
		
			ENDIF
		ENDIF
    ENDIF
ENDDO

OPEN(1, ACTION='WRITE',FILE='Wall morphology11.dat')
WRITE(1,*) 'TITLE = "Field Plot"'
WRITE(1,*) 'VARIABLES = "x" "y"'
DO j=1,N_Boundary
IF (boundaries(j)%Erosion > 0) THEN
	DO i=1,boundaries(j)%Wall(1)%nxp
		WRITE(1,50) boundaries(j)%Wall(1)%Node(i,1), boundaries(j)%Wall(1)%Node(i,2)
	END DO
ENDIF
ENDDO
50 FORMAT (E15.6,' ',E15.6)
CLOSE(1)


!!DEALLOCATE(wallindex_tmp,nnode)
!!DEALLOCATE(ynode)

n_sect= MAXVAL(objects(:)%Wall(1)%nxp)
ALLOCATE(alphax(n_sect,N_Objects,ispe_tot),E_par(n_sect,N_Objects,ispe_tot),nc_par(n_sect,N_Objects,ispe_tot))		!在运行完一次惠根斯子波求解下降距离后需释放或清零这几个数组的内存 
alphax=0
E_par=0
nc_par=0

!!!!!!!!!!!!!!!!!! prepare for recording!!!!!!!!!!!!!!!!!!
WRITE(6,*) 'Setup_IFE_Wall_Mesh_2D  Done!'
END