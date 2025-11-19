Module ShellSort_2D
    Implicit None
    
    Contains
    
    Subroutine ShellSort(array, index)
    
        Integer :: N, index
        Integer :: i, j
        Integer :: insertIndex
        Real(8), Dimension(:,:) :: array
        Real(8) :: insertVal(SIZE(array,1))
        
        N = SIZE(array, 2)
        
        i = N/2
        
        Do while(i >= 1)
            
            Do  j = i, N - 1, 1
                
                insertVal = array(:, j + 1)
                insertIndex = j - i
                Do While (insertIndex >= 0)
                    If (insertVal(index) < array(index, insertIndex + 1)) Then
                        array(:, insertIndex + i + 1) = array(:, insertIndex + 1)
                        insertIndex = insertIndex - i
                    Else
                        exit
                    End If
                End Do
                insertIndex = insertIndex + i
                array(:, insertIndex + 1) = insertVal
            End Do
            
            i = i/2
        End Do
        
    End Subroutine
    
End Module ShellSort_2D