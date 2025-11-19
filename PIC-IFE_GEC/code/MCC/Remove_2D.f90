SUBROUTINE Remove_2D(i_part, nnp)

! Updated:	4/28/2005 11:15 AM
! Purpose:	Remove particle from the particles array.
!			This subroutine does not output removed particles.

USE Particle_2D
USE MCC_Data_2D

IMPLICIT NONE

INTEGER		i_part, nnp

INTEGER			ispe, i
		

	ispe = part(10,i_part)
	part(:,i_part) = part(:,nnp)
	ENER(i_part) = ENER(nnp)		

	ns(ispe) = ns(ispe) - 1
	nnp = nnp-1

END
