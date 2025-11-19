MODULE IFE_MAIN_PARAM

IMPLICIT NONE

REAL(8), PARAMETER	::	pi	= 3.14159265358979D0

! Sparse matrix structure
TYPE SPARSE
	REAL(8)	::	K
	INTEGER	::	JCOL, SROW
END TYPE

! Nonlinear Solver Option
INTEGER		::	NSolver		!	1 Gasuss-Seidel,	2	Newton

! Sparse Linear Solver Option
INTEGER		::	SLSolver	!	1 IMSL PCCG,		2	DLAP PCCG

! Fixed/Floating Nodes Blocking Option
INTEGER		::	Blocking	!	1	YES,	0	NO

! PCG Solver Parameters
INTEGER		::	PCG_MaxIT
REAL(8)		::	PCG_Tol

! Nonlinear Block Gauss-Seidel Solver Parameters
INTEGER		::	BGS_MaxIT
REAL(8)		::	BGS_Tol

! Nonlinear Block Gauss-Seidel Solver Parameters
INTEGER		::	Newton_MaxIT
REAL(8)		::	Newton_Tol

! Mesh Intersection Parameters
REAL(8)		::	Int_El_Frac

! Sparse Structure Parameters
INTEGER, PARAMETER	::	Max_Nodal_Connect = 200

REAL(8), PARAMETER	::	SmallValue = 1.0D-12 

! Machine Numerical Constants
! PC Version
REAL(8), PARAMETER	::	MZero = 2.*EPSILON(1.D0), Zero = 0.D0, One = 1.D0, Half = 0.5D0, Quart = 0.25D0, Three=3.D0, Six = 6.D0

! CRAY Version
!REAL(8), PARAMETER	::	MZero = 2.*EPSILON(1.E0), Zero = 0.E0, One = 1.E0, Half = 0.5E0, Quart = 0.25E0, Six = 6.E0



END MODULE