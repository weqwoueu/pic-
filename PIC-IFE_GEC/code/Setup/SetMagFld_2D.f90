SUBROUTINE SetMagFld_2D(delta)
    USE Constant_Variable_2D
    USE Domain_2D
    USE Field_2D
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: delta
    
    INTEGER :: I, J, K
    INTEGER :: nx_femm_mag, ny_femm_mag
    INTEGER :: temp_nx, temp_ny
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: BX_IMP, BY_IMP, BZ_IMP, B_INPUT
    LOGICAL :: ALIVE

    WRITE(6,*) 'SetMagFld'
    
    nx_femm_mag=141
    ny_femm_mag=191
    
    temp_nx = 141 !877-266+1
    temp_ny = 191 !901-1+1
    
    INQUIRE(FILE='./INPUT/B.DAT', EXIST=ALIVE)
    IF (.NOT.ALIVE) THEN
        WRITE(6,*) 'ERROR: INPUT FILE IS MISSING2 !'
        PAUSE
        STOP
    ELSE
        ALLOCATE(BX_IMP(temp_nx, temp_ny), BY_IMP(temp_nx, temp_ny), BZ_IMP(temp_nx, temp_ny))
        ALLOCATE(B_INPUT(2, temp_nx*temp_ny))
        
        OPEN(1, FILE='./INPUT/B.DAT',STATUS='OLD')
            READ(1,*) B_INPUT
        CLOSE(1)
        
        K = 0
        DO I = 1, temp_nx
            DO J = 1, temp_ny
                K = K + 1
                BX_IMP(I,J) = B_INPUT(1,K)
                BY_IMP(I,J) = B_INPUT(2,K)
                BZ_IMP(I,J) = 0.
            ENDDO
        ENDDO
        DEALLOCATE(B_INPUT)
        
    ENDIF
    
    DO J = 1, ny
        DO I = 1, nx
            BFX(I,J) = BX_IMP(I,J)
            BFY(I,J) = BY_IMP(I,J)
            BFZ(I,J) = 0.0
        ENDDO
    ENDDO
    
    DEALLOCATE(BX_IMP, BY_IMP)
    
    DO j=1,ny
        DO i=1,nx
	        bt(i,j)=sqrt(bfx(i,j)*bfx(i,j)+bfy(i,j)*bfy(i,j)+bfz(i,j)*bfz(i,j))
        END DO
    END DO
    
    DO J = 1, ny
        DO I = 1, nx
            BFX(I,J) = BFX(I,J)/Bfield_ref
            BFY(I,J) = BFY(I,J)/Bfield_ref
            BFZ(I,J) = BFZ(I,J)/Bfield_ref
        ENDDO
    ENDDO
    
    
    OPEN(17, ACTION='WRITE', FILE='./OUTPUT/MagFld.dat')
        WRITE(17,*) 'TITLE = "Field Plot"'
        WRITE(17,*) 'VARIABLES = "x" "y" "bx" "by" "bz" "bt"'
        WRITE(17,"(' ZONE I= ',I6,', J= ',I6)") nx, ny
        DO j = 1, ny
            DO i = 1, nx
                WRITE(17,"(6(E15.6))") VertX(1:2,i,j), bfx(i,j), bfy(i,j), bfz(i,j), bt(i,j)
            ENDDO
        ENDDO
    CLOSE(17)
    
END SUBROUTINE