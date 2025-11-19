MODULE FluxLine_2D

IMPLICIT NONE

TYPE FluxLineType
	
	REAL(8) line_start(2)

	REAL(8) line_end(2)

	REAL(8) Q

	INTEGER num_elec, num_ion

	REAL(8) length

END TYPE




END MODULE