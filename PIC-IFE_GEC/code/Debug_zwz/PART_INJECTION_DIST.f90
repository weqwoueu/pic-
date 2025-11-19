SUBROUTINE PART_INJECTION_DIST(ILOOPA)
    !USE CONSTANT
    !USE GLOBAL
    USE MDL_DEBUG
    !USE Var_newsetup
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    INTEGER, INTENT(IN) :: ILOOPA
    INTEGER :: I
    
    REAL(4) :: delta_p
    REAL(4) :: min_p
    REAL(4) :: max_p
    
    REAL(4) :: DELTA_V
    REAL(4) :: MIN_V
    REAL(4) :: MAX_V
    REAL(4) :: SUM_DIST_FVA 
    
    IF (DEBUG==1 .AND. PART_STAT_COUNT==PART_STAT_NUM) THEN
        
        !PART_stat_xp = PART_stat_xp
        !PART_stat_yp = PART_stat_yp
        !
        !PART_STAT_VX = PART_STAT_VX * VETH /SQRT(M_FACTOR(ISPA))
        !PART_STAT_VY = PART_STAT_VY * VETH /SQRT(M_FACTOR(ISPA))
        !PART_STAT_VZ = PART_STAT_VZ * VETH /SQRT(M_FACTOR(ISPA))
 
        !$------------ x position -----------
        max_p = MAXVAL(PART_stat_xp)
        min_p = MINVAL(PART_stat_xp)
        delta_p = (max_p-min_p) /real(df_size-1)
        DO I = 1, DF_SIZE-1
            dist_pax(i) = min_p+real(i-1)*delta_p
        END DO
        dist_pax(df_size) = max_p
        CALL DF_STATS(PART_STAT_NUM, PART_STAT_XP, DF_SIZE, DIST_PAX, DIST_FPAX)
        
        !$------------ y position -----------
        max_p = MAXVAL(PART_stat_yp)
        min_p = MINVAL(PART_stat_yp)
        delta_p = (max_p-min_p) /real(df_size-1)
        DO I = 1, DF_SIZE-1
            dist_pay(i) = min_p+real(i-1)*delta_p
        END DO
        dist_pay(df_size) = max_p
        CALL DF_STATS(PART_STAT_NUM, PART_STAT_YP, DF_SIZE, DIST_PAY, DIST_FPAY)
        
        !$------------ x velocity -----------
        MAX_V = MAXVAL(PART_STAT_VX)
        MIN_V = MINVAL(PART_STAT_VX)
        DELTA_V = (MAX_V-MIN_V) /REAL(DF_SIZE-1)
        DO I = 1, DF_SIZE-1
            DIST_VAX(I) = MIN_V+REAL(I-1)*DELTA_V
        END DO
        DIST_VAX(DF_SIZE) = MAX_V
        CALL DF_STATS(PART_STAT_NUM, PART_STAT_VX, DF_SIZE, DIST_VAX, DIST_FVAX)
                
        !$------------ y velocity -----------
        MAX_V = MAXVAL(PART_STAT_VY)
        MIN_V = MINVAL(PART_STAT_VY)
        DELTA_V = (MAX_V-MIN_V) /REAL(DF_SIZE-1)
        DO I = 1, DF_SIZE-1
            DIST_VAY(I) = MIN_V+REAL(I-1)*DELTA_V
        END DO
        DIST_VAY(DF_SIZE) = MAX_V
        CALL DF_STATS(PART_STAT_NUM, PART_STAT_VY, DF_SIZE, DIST_VAY, DIST_FVAY)
        
        !$------------ z velocity -----------
        MAX_V = MAXVAL(PART_STAT_VZ)
        MIN_V = MINVAL(PART_STAT_VZ)
        DELTA_V = (MAX_V-MIN_V) /REAL(DF_SIZE-1)
        DO I = 1, DF_SIZE-1
            DIST_VAZ(I) = MIN_V+REAL(I-1)*DELTA_V
        END DO
        DIST_VAZ(DF_SIZE) = MAX_V
        CALL DF_STATS(PART_STAT_NUM, PART_STAT_VZ, DF_SIZE, DIST_VAZ, DIST_FVAZ)
        
        write(6,*)'writing PART distribution to file...'
        OPEN(40,FILE='./PART_DIST.DAT',STATUS='REPLACE')
            WRITE(40,"(A150)") ' VARIABLES = "XP" "FPX" "YP" "FPY" "VX" "FVX" "VY" "FVY" "VZ" "FVZ"'
            WRITE(40,801) ILOOPA
801			FORMAT('TITLE = "STEP =',I8,'"')
            DO I = 1, DF_SIZE-1
                WRITE(40,"(10(E15.5,1X))")  (DIST_PAX(I)+DIST_PAX(I+1))/2., DIST_FPAX(I+1), &
                                            (DIST_PAY(I)+DIST_PAY(I+1))/2., DIST_FPAY(I+1), &
                                            (DIST_VAX(I)+DIST_VAX(I+1))/2., DIST_FVAX(I+1), &
                                            (DIST_VAY(I)+DIST_VAY(I+1))/2., DIST_FVAY(I+1), &
                                            (DIST_VAZ(I)+DIST_VAZ(I+1))/2., DIST_FVAZ(I+1)
            END DO
		CLOSE(40)
        
        !pause
        
        PART_STAT_COUNT = 0.
        
        !$------------------
        PART_STAT_XP = 0.
        PART_STAT_YP = 0.
        !$------------------
        
        PART_STAT_VX = 0.
        PART_STAT_VY = 0.
        PART_STAT_VZ = 0.
    END IF
    
END SUBROUTINE PART_INJECTION_DIST