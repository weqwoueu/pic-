SUBROUTINE generate_HE_2D(HP, HT, ele_col_number, HE)
! wsy add SIDG 2021/10/27

IMPLICIT NONE

REAL(8), DIMENSION(:,:), POINTER 	    ::	HP
INTEGER, DIMENSION(:,:), POINTER	    ::	HT, HE, HE_temp
INTEGER                               ::  ele_col_number

REAL                                  ::  top, bottom, left, right
INTEGER                               ::  count, i, j

left = MINVAL(HP(1, :))
right = MAXVAL(HP(1, :))
top = MAXVAL(HP(2, :))
bottom = MINVAL(HP(2, :))

count = 0
ALLOCATE(HE_temp(8, SIZE(HP,2)))
HE_temp = 0

DO i = 1, SIZE(HT,2)
    
    IF (HT(5, i) > 0) THEN
        
        DO j = 1, 4
            
            
            IF (j /= 4) THEN
                IF ((ABS(HP(1, HT(j, i))-left)<10D-6 .AND. ABS(HP(1, HT(j + 1, i))-left)<10D-6) .OR. &
                    (ABS(HP(1, HT(j, i))-right)<10D-6 .AND. ABS(HP(1, HT(j + 1, i))-right)<10D-6) .OR. &
                    (ABS(HP(2, HT(j, i))-top)<10D-6 .AND. ABS(HP(2, HT(j + 1, i))-top)<10D-6) .OR. &
                    (ABS(HP(2, HT(j, i))-bottom)<10D-6 .AND. ABS(HP(2, HT(j + 1, i))-bottom)<10D-6)) THEN
                ! IF ( (HP(1, HT(j, i)) == left .AND. HP(1, HT(j + 1, i)) == left) .OR. &
                !     (HP(1, HT(j, i)) == right .AND. HP(1, HT(j + 1, i)) == right) .OR. &
                !     (HP(2, HT(j, i)) == top .AND. HP(2, HT(j + 1, i)) == top) .OR. &
                !     (HP(2, HT(j, i)) == bottom .AND. HP(2, HT(j + 1, i)) == bottom) ) THEN
                    CONTINUE
                ELSE
                    
                    count = count + 1
                    HE_temp(1, count) = HT(j, i)
                    HE_temp(2, count) = HT(j + 1, i)
                    HE_temp(5, count) = i
                    IF (j == 1) THEN
                        HE_temp(3, count) = 0
                        HE_temp(4, count) = -1
                        HE_temp(6, count) = HE_temp(5, count) - 1
                    ELSEIF (j == 2) THEN
                        HE_temp(3, count) = -1
                        HE_temp(4, count) = 0
                        HE_temp(6, count) = HE_temp(5, count) + ele_col_number
                    ELSEIF (j == 3) THEN
                        HE_temp(3, count) = 0
                        HE_temp(4, count) = -1
                        HE_temp(6, count) = HE_temp(5, count) + 1
                    END IF
                    
                END IF
                
            ELSE    
                
                IF ((ABS(HP(1, HT(4, i))-left)<10D-6 .AND. ABS(HP(1, HT(1, i))-left)<10D-6) .OR. &
                    (ABS(HP(1, HT(4, i))-right)<10D-6 .AND. ABS(HP(1, HT(1, i))-right)<10D-6)) THEN
                    CONTINUE
                ELSE
                    
                    count = count + 1
                    HE_temp(1, count) = HT(4, i)
                    HE_temp(2, count) = HT(1, i)
                    HE_temp(5, count) = i
                    HE_temp(3, count) = -1
                    HE_temp(4, count) = 0
                    HE_temp(6, count) = HE_temp(5, count) - ele_col_number
                    
                END IF
                
            END IF
            
        END DO
        
    END IF
    
END DO

AllOCATE(HE(8, count))
DO i = 1, count
    
    HE(:, i) = HE_temp(:, i)
    
END DO

DEALLOCATE(HE_temp)

END

    