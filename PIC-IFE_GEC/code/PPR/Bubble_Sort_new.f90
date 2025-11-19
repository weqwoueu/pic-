SUBROUTINE Bubble_Sort_new(A,N)

! wsy revise to shell sort
    
IMPLICIT NONE

Integer :: N
        Integer :: A(N)
        Integer :: i, j
        Integer :: insertVal, insertIndex
        
        i = N/2
        
        Do while(i >= 1)
            
            Do  j = i, N - 1, 1
                
                insertVal = A(j + 1)
                insertIndex = j - i
                Do While (insertIndex >= 0)
                    If (insertVal < A(insertIndex + 1)) Then
                        A(insertIndex + i + 1) = A(insertIndex + 1)
                        insertIndex = insertIndex - i
                    Else
                        exit
                    End If
                End Do
                insertIndex = insertIndex + i
                A(insertIndex + 1) = insertVal
            End Do
            
            i = i/2
        End Do

END SUBROUTINE