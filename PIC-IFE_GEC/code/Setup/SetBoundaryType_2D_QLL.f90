SUBROUTINE SetBoundaryType_2D_QLL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: SetBoundaryType_2D.f90                             C
!
!  Purpose: Set outer boundary type for the box
!                                                                      C
!  Reviewer: Yuchuan Chu                              Date: 05-May-12  C
!  Comments: modified for normalization data input in input_2D.f90     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

USE PIC_MAIN_PARAM_2D
USE Domain_2D

IMPLICIT NONE

INTEGER		iface, idim, ifaceindex

WRITE(6,*)
WRITE(6,*) 'SetBoundaryType'

! ------------setup the domain
! default inner face, periodic, 

!	DO iface=1,6
!	outerface(iface)=.FALSE.
!	periodic(iface)=.TRUE.
!	END DO

! field

! check the consistency
DO iface = 1,4
! field
	IF(f_periodic(iface).AND.f_zeroe(iface))	WRITE(6,*) 'Warning: f_zeroe at iface =', iface
! particles
	IF(periodic(iface)) THEN
		IF(pabsorb(iface))			WRITE(6,*) 'Warning: pabsorb at iface=',iface
		IF(preflect(iface))			WRITE(6,*) 'Warning: preflect at iface=',iface
		IF(pemit(iface))			WRITE(6,*) 'Warning: pemit at iface=',iface
	END IF
	IF(pabsorb(iface)) THEN 
		IF(preflect(iface))			WRITE(6,*) 'Warning: pabsorb & preflect at iface =', iface
	END IF
! particle and field
	IF(f_periodic(iface)) THEN
		IF(.NOT.periodic(iface))	WRITE(6,*) 'Warning: periodic & f_periodic at iface =',iface
	END IF 
END DO


DO idim = 1,2
	ifaceindex = 2*idim
	WRITE(6,*) 'i_space = ', idim
	WRITE(6,*) 'f_zeroe   ', f_zeroe(ifaceindex-1), f_zeroe(ifaceindex)
	WRITE(6,*) 'f_periodic', f_periodic(ifaceindex-1), f_periodic(ifaceindex)
END DO

DO idim = 1,3
	ifaceindex = 2*idim
	WRITE(6,*) 'i_space = ', idim
	WRITE(6,*) 'periodic  ', periodic(ifaceindex-1), periodic(ifaceindex)
	WRITE(6,*) 'pabsorb   ', pabsorb(ifaceindex-1), pabsorb(ifaceindex)
	WRITE(6,*) 'preflect  ', preflect(ifaceindex-1), preflect(ifaceindex)
	WRITE(6,*) 'pemit     ', pemit(ifaceindex-1), pemit(ifaceindex)
END DO
!DO idim = 1,2
!	ifaceindex = 2*idim
!	WRITE(6,*) 'i_space = ', idim
!	WRITE(6,*) 'f_zeroe   ', f_zeroe(ifaceindex-1), f_zeroe(ifaceindex)
!	WRITE(6,*) 'f_periodic', f_periodic(ifaceindex-1), f_periodic(ifaceindex)
!	WRITE(6,*) 'periodic  ', periodic(ifaceindex-1), periodic(ifaceindex)
!	WRITE(6,*) 'pabsorb   ', pabsorb(ifaceindex-1), pabsorb(ifaceindex)
!	WRITE(6,*) 'preflect  ', preflect(ifaceindex-1), preflect(ifaceindex)
!	WRITE(6,*) 'pemit     ', pemit(ifaceindex-1), pemit(ifaceindex)
!END DO

END