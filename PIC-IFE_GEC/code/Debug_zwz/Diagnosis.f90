SUBROUTINE DF_STATS(PSIZE, PARA, DFSIZE, XPARA, FPARA)
    ! PSIZE:            样本的个数
    ! PARA(PSIZE):      样本的参数集
    ! DFSIZE-1:         统计参数分布用到的区间个数
    ! XPARA(DFSIZE):    分布区间节点上的参数值
    ! FPARA(2:DFSIZE):  第1到DFSIZE-1个区间上的分布概率；FPARA(1)=0；SUM(FPARA)=1

    IMPLICIT NONE
    INTEGER :: I, J
    INTEGER :: PSIZE, DFSIZE
    REAL(4) :: PARA(PSIZE), XPARA(DFSIZE), FPARA(DFSIZE)
    REAL(4) :: TEMP

    CALL QUICK_SORT(PARA(1:PSIZE), PSIZE, 1, PSIZE)
    
    J = 2
    DO I = 1, PSIZE
200     IF (PARA(I) < XPARA(J)) THEN
            FPARA(J) = FPARA(J) +1.
        ELSE
            J = J+1
            IF (J <= DFSIZE) THEN
                GOTO 200
            ELSE
                EXIT
            END IF
        END IF
    END DO
    
    !FPARA(DFSIZE) = REAL(PSIZE - I)
    FPARA(1) = 0.
    
    TEMP = REAL((PSIZE-I))/REAL(PSIZE)
    IF (TEMP > 0.05) THEN
        WRITE(*,"('POPULATION OUT OF STATISTICS IS ', F10.4, '%')") TEMP*100.
        PAUSE
    END IF
    
    FPARA = FPARA/I
    DO J = 2, DFSIZE
        FPARA(J) = FPARA(J)/(XPARA(J)-XPARA(J-1))
    END DO
    
END SUBROUTINE DF_STATS
    
! 快速排序法的子程序
RECURSIVE SUBROUTINE QUICK_SORT(A,N,S,E)
    IMPLICIT NONE
    
    INTEGER :: N ! 表示类型的大小
    INTEGER :: S ! 传入的参数, 这一组的类型起始位置
    INTEGER :: E ! 传入的参数, 这一组的类型结束位置
    INTEGER :: L,R ! 用来找A(L)>K及A(R)<K时用的
    REAL(4) :: A(N) ! 存放数据的类型
    REAL(4) :: K ! 记录键值A(S)
    REAL(4) :: TEMP ! 交换两个数值时用的
    
    ! 首先要先给定L,R的初值. L要从头开始,E则要从尾开始
    L=S
    R=E
    ! RIGHT值 > LEFT值 时才有必要进行排序 
    IF ( R<=L ) RETURN
    
    K=A(S) ! 设定键值
    DO WHILE(.TRUE.)
        ! 找出A(L)<K的所在
        DO WHILE( .TRUE. )
            L=L+1
            IF ( (A(L) > K) .OR. (L>=E) ) EXIT
        END DO
        ! 找出A(R)>K的所在
        DO WHILE( .TRUE. )
            R=R-1
            IF ( (A(R) < K) .OR. (R<=S) ) EXIT
        END DO
        
        ! 如果RIGHT 跑到 LEFT的左边时, 循环就该结束了
        IF ( R <= L ) EXIT
        ! 交换A(L),A(R)的数值
        TEMP=A(L)
        A(L)=A(R)
        A(R)=TEMP
    END DO
    
    ! 交换A(S),A(R)的数值
    TEMP=A(S)
    A(S)=A(R)
    A(R)=TEMP
    ! 把R之前的数据重新分组,再做排序
    CALL QUICK_SORT(A,N,S,R-1)
    ! 把R之后的数据重新分组,再做排序
    CALL QUICK_SORT(A,N,R+1,E)
    
    RETURN
END SUBROUTINE QUICK_SORT