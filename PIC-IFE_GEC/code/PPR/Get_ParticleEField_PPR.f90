Subroutine Get_ParticleEField_PPR(x, y, EleIndex, Ex, Ey)
    ! wsy add for ppr
    Use IFE_Data
    Implicit None
    !---  in out-------
    Real(8), Intent(In) :: x, y
    Integer, Intent(In) :: EleIndex
    Real(8), Intent(Out) :: Ex, Ey
    !------------------
    !--- temp parameter -------
    Integer,parameter :: NodeInEle = 4
    Integer :: i, j, k, node_temp, count
    Real(8) :: x_ref, y_ref
    !--------------------------
    
    Ex = 0
    Ey = 0
    
    If (element_index(EleIndex) > 0) Then
        
        count = 0
        Do i = 1, NodeInEle
            If (node_index(HT(i, EleIndex)) == -1) Then
                node_temp = HT(i, EleIndex)
                count = count + 1
                x_ref = X - HP(1, node_temp)
                y_ref = Y - HP(2, node_temp)
                Ex = Ex + coef_PPR(2,node_temp) + 2.* coef_PPR(4,node_temp)*x_ref + coef_PPR(5,node_temp)*y_ref
                Ey = Ey + coef_PPR(3,node_temp) + 2.* coef_PPR(6,node_temp)*y_ref + coef_PPR(5,node_temp)*x_ref
            End If
        End Do
        Ex = Ex/count
        Ey = Ey/count
        
    Else
        
        Do i = 1, NodeInEle
                node_temp = HT(i, EleIndex)
                x_ref = X - HP(1, node_temp)
                y_ref = Y - HP(2, node_temp)
                Ex = Ex + coef_PPR(2,node_temp) + 2.* coef_PPR(4,node_temp)*x_ref + coef_PPR(5,node_temp)*y_ref
                Ey = Ey + coef_PPR(3,node_temp) + 2.* coef_PPR(6,node_temp)*y_ref + coef_PPR(5,node_temp)*x_ref
        End Do
        Ex = Ex/NodeInEle
        Ey = Ey/NodeInEle
        
    End If
    
    
End Subroutine