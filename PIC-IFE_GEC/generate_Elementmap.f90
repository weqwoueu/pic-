Subroutine generate_Elementmap(Element_Map, repeat_refinement)
    Use Domain_2D
    Use IFE_INTERFACE, Only: ParticlePositioning
    
    Implicit None
    
    Integer, Allocatable, Intent(InOut) ::  Element_Map(:,:)
    Integer, Intent(In)                 ::  repeat_refinement
    
    Integer, Allocatable                ::  EdgeArray(:)
    Integer                             ::  n_element, i, j, k, l, traverse_flag
    Real(8)                             ::  dx, dy, x, y
    REAL(8)		                        ::  rxp, ryp
    
    traverse_flag = 0
    n_element = 0
    
    
    Allocate(Element_Map((2**repeat_refinement)*(nx - 1), (2**repeat_refinement)*(ny - 1)))
    
    dx = hx(1) / (2**repeat_refinement)
    dy = hx(2) / (2**repeat_refinement)
    
    Do i = 1, (2**repeat_refinement)*(nx - 1)
        Do j = 1, (2**repeat_refinement)*(ny - 1)
            
            x = dxmin + (dx) * (2*i - 1) / 2
            y = dymin + (dy) * (2*j - 1) / 2
            
            rxp = (x - Vert_o(1)) / hx(1)
            k = INT(rxp)
    
            ryp = (y - Vert_o(2)) / hx(2)
            l = INT(ryp)
    
            If (l == ny) Then
              l = l - 1
            End If
    
            If (k == nx) Then
              k = k - 1
            End If
    
            n_element = (l + (k-1)*(ny-1))
            
            If (n_element < 0) Then
                PRINT*, 'Case4: Check Partition and Code, Stop'
                pause
            End If
            
            Call ParticlePositioning(x, y, n_element, traverse_flag, EdgeArray)
            
            If ( Allocated(EdgeArray) ) Deallocate(EdgeArray)
            
            Element_Map(i, j) = n_element
            
        End Do
    End Do
    
End Subroutine
    