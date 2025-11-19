SUBROUTINE Collision_Cross_Section_Table
    
!!! bjw add 2019-03-30 读入截面图上的离散数据  

USE MCC_Data_2D 	
IMPLICIT NONE

INTEGER  :: i, number_ela, number_exc, number_ion

WRITE(*,*) 'Collision_Cross_Section_Table'

OPEN(1,ACTION = 'READ', FILE ='Cross_Section_ela.txt' )
    READ(1,*) number_ela   !!! elastic
    ALLOCATE(section_ela(2,number_ela))
    DO i = 1,number_ela
        READ(1,*) section_ela(1,i), section_ela(2,i)
    END DO
CLOSE(1)

OPEN(1,ACTION = 'READ', FILE ='Cross_Section_exc.txt' )
    READ(1,*) number_exc   !!! excitation
    ALLOCATE(section_exc(2,number_exc))
    DO i = 1,number_exc
        READ(1,*) section_exc(1,i), section_exc(2,i)
    END DO
CLOSE(1)

OPEN(1,ACTION = 'READ', FILE ='Cross_Section_ion.txt' )
    READ(1,*) number_ion  !!! ionization
    ALLOCATE(section_ion(2,number_ion))
    DO i = 1,number_ion
        READ(1,*) section_ion(1,i), section_ion(2,i)
    END DO
CLOSE(1)

END