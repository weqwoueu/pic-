Subroutine InjectPositioning(part_x, part_y, n_element)
    Use IFE_Data
    Use Domain_2D
    Implicit None
    
    Real(8), Intent(In)     :: part_x, part_y
    Integer, Intent(InOut)  :: n_element
    
    Real(8)                 :: dx, dy, x, y, rxp, ryp
    Integer                 :: i, j
    
    
    dx = hx(1) / (2**repeat_refinement)
    dy = hx(2) / (2**repeat_refinement)
    
    rxp = (part_x - dxmin) / dx
    i = INT(rxp) + 1
    
    ryp = (part_y - dymin) / dy
    j = INT(ryp) + 1
    
    n_element = Element_Map(i, j)
    
    
End Subroutine