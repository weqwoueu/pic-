Subroutine MoveOne(PO,isp,TimeMove, IvelFlag, IposFlag, delta, OriPosi)
    Use Domain_2D
    Use Field_2D
    Use Particle_2D, Only: qm
    Use IMPIC_Data_2D, ONLY: Bfiled_index, Omega
    Use ModuleParticleOne
    !============LY add for SIDG-PIC couple, 2022-7-25============
    Use IFE_Data
    !============LY add for SIDG-PIC couple, 2022-7-25============
    Implicit none
    
    !Class(ParticleOne), intent(inout) :: PO
    Type(ParticleOne), intent(inout) :: PO
    Integer(4), intent(in) :: isp
    Real(8), intent(in) :: TimeMove
    Integer(4), intent(inout) :: IvelFlag, IposFlag
    Integer(4), intent(in) :: delta
    Real(8), intent(out) :: OriPosi(3)
    
    Real(8) :: xp, yp
    Real(8) :: dx, dy, xcellmdx, ycellmdy
    Real(8) :: Rp, R1, R2, den
    Real(8) :: P1, P2, P3, P4
    Integer(4) :: i, j
    
    Real(8) :: f
    Real(8) :: Efield(3), Bfield(3)
    Real(8) :: SVelocity(1:3)=0.d0
    Real(8) :: deta_vx, deta_vy, deta_vz
    REAL(8) :: X, Y, Z
    REAL(8) :: zefield, refield, tefield
    REAL(8) :: xefield, yefield
    REAL(8) :: zbfield, rbfield, tbfield
    REAL(8) :: xbfield, ybfield
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    Integer :: n_element_old
    Real(8) :: hx_partition, hy_partition
    Real(8) :: W1, W2, W3, W4
    Real(8) :: ex_part, ey_part
    Integer :: trial_basis_type
    Real(8) :: vertices(2,4)
    Integer :: index_ele, piece_flag, n
    Real(8) :: beta1, beta2
    Integer :: information_vector_1(18)
    Real(8) :: information_vector_2(8)
    Real(8) :: Ex_basis, Ey_basis
    Real(8) :: Efield_old(3), Efield_bro(3)
    Real(8) :: Bfield_old(3), Bfield_bro(3)
    Real(8), Parameter :: PI_LY	= 3.14159265358979D0
    
    !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    
    
    OriPosi = (/PO%X,PO%Y,PO%Z/)
    
    If(IvelFlag == 1) Then
      
      !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
      
      !==========Part 1 : Positioning Particle and Moving==========
      
      !xp = (PO%X - Vert_o(1))*hxi(1)
      !i = INT(xp)
      !
      !yp = (PO%Y - Vert_o(2))*hxi(2)
      !j = INT(yp)
      !
      !If (j == ny) Then
      !  j = j - 1
      !End If
      !
      !If (i == nx) Then
      !  i = i -1
      !End If
      
        !n_element_old = (j + (i-1)*(ny-1))
        n_element_old = PO%Location

        !There denote the particle locate inside the element.
        hx_partition = HP(1, HT(2,n_element_old)) - HP(1, HT(1,n_element_old))
        hy_partition = HP(2, HT(4,n_element_old)) - HP(2, HT(1,n_element_old))

        dx = (PO%X - HP(1, HT(1,n_element_old))) / hx_partition
        dy = (PO%Y - HP(2, HT(1,n_element_old))) / hy_partition

        xcellmdx = 1.0 - dx
        ycellmdy = 1.0 - dy

        If (delta == 0) Then  !2D Cartesian coordinates.
            W1 = xcellmdx * ycellmdy
            W2 = dx       * ycellmdy
            W3 = dx       * dy
            W4 = xcellmdx * dy
        Elseif (delta == 1) Then  !2D axis-symmetric coordinates.
            R1 = dymin + HP(2,HT(1,n_element_old))
            R2 = dymin + HP(2,HT(4,n_element_old))
            Rp = PO%Y
            den = R2*R2 - R1*R1

            !Here, we not consider the unequal weights for the moment, 2022-4-8.
            W1 = xcellmdx * (R2*R2-Rp*Rp) / den
            W2 = dx       * (R2*R2-Rp*Rp) / den
            W3 = dx       * (Rp*Rp-R1*R1) / den
            W4 = xcellmdx * (Rp*Rp-R1*R1) / den
        End If

        !Firstly, we need to calculate the Efield between interface element and non-interface element.
        If (element_index(n_element_old) >0) Then  !There denote the particle locate on interface element. We need use basis function interpolation.
            ex_part = 0.0
            ey_part = 0.0
          
            Call Get_ParticleEField_PPR(PO%X, PO%Y, n_element_old, ex_part, ey_part)
            
            !trial_basis_type = 1
            !vertices = HP(1:2, HT(1:4, n_element_old))
            !index_ele = element_index(n_element_old)
            !beta1 = information_2(1, index_ele)
            !beta2 = information_2(2, index_ele)
            !information_vector_1 = information_1(1:18, index_ele)
            !information_vector_2 = information_2(1:8, index_ele)
            !
            !If (ABS(beta1 - out_object_beta) < SmallValue) Then
            !  piece_flag = 1
            !Elseif (ABS(beta2 - out_object_beta) < SmallValue) Then
            !  piece_flag = 2
            !End If
            !
            !Do n = 1, 4
            !  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
            !                                  piece_flag, trial_basis_type, n, 1, 0, Ex_basis, n_element_old)
            !  ex_part = ex_part + Phi(HT(n, n_element_old), 1) * Ex_basis
            !
            !  Call Retangular_local_basis_IFE(PO%X, PO%Y, vertices, information_vector_1, information_vector_2,&
            !                                  piece_flag, trial_basis_type, n, 0, 1, Ey_basis, n_element_old)
            !  ey_part = ey_part + Phi(HT(n, n_element_old), 1) * Ey_basis
            !End Do

            Efield(1) = -ex_part
            Efield(2) = -ey_part

        Else  !There denote the particle locate on non-interface element. We need use linear interpolation.
            Efield(1) = W1 * efx(HT(1,n_element_old),1) + W2 * efx(HT(2,n_element_old),1) + &
                        W3 * efx(HT(3,n_element_old),1) + W4 * efx(HT(4,n_element_old),1)
            Efield(2) = W1 * efy(HT(1,n_element_old),1) + W2 * efy(HT(2,n_element_old),1) + &
                        W3 * efy(HT(3,n_element_old),1) + W4 * efy(HT(4,n_element_old),1)
        End If
        !=========LY modification for 2D3V PIC model, 2022-5-27=========
        Efield(3) = W1 * efz(HT(1,n_element_old),1) + W2 * efz(HT(2,n_element_old),1) + &
                    W3 * efz(HT(3,n_element_old),1) + W4 * efz(HT(4,n_element_old),1)
        !=========LY modification for 2D3V PIC model, 2022-5-27=========

        !Secondly, we need to calculate the Efield force and Bfield force.
        If (Bfiled_index) Then  !with magnetic field.
            If (delta_global == 0) Then  !2D Cartesian coordinates.
            deta_vx = 0.5 * qm(isp+1) * Efield(1) * TimeMove
            deta_vy = 0.5 * qm(isp+1) * Efield(2) * TimeMove
            deta_vz = 0.5 * qm(isp+1) * Efield(3) * TimeMove

            Bfield(1) = W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                        W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)
            Bfield(2) = W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                        W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)
            Bfield(3) = W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                        W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)

            Omega(1) = Bfield(1) * qm(isp+1) * 0.5 * TimeMove
            Omega(2) = Bfield(2) * qm(isp+1) * 0.5 * TimeMove
            Omega(3) = Bfield(3) * qm(isp+1) * 0.5 * TimeMove

            Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
            !$ axisymmetric e field
            zefield = Efield(1)
            refield = Efield(2)
            tefield = Efield(3)
            !$ Convert axisymmetric efield to Cartesian efield
            xefield = refield*DCOS(PO%Z) - tefield*DSIN(PO%Z)
            yefield = refield*DSIN(PO%Z) + tefield*DCOS(PO%Z)
          
            deta_vx = 0.5 * qm(isp+1) * zefield * TimeMove
            deta_vy = 0.5 * qm(isp+1) * xefield * TimeMove
            deta_vz = 0.5 * qm(isp+1) * yefield * TimeMove

            !$ axisymmetric b field
            zbfield = W1 * bfx(HT(1,n_element_old),1) + W2 * bfx(HT(2,n_element_old),1) + &
                        W3 * bfx(HT(3,n_element_old),1) + W4 * bfx(HT(4,n_element_old),1)
            rbfield = W1 * bfy(HT(1,n_element_old),1) + W2 * bfy(HT(2,n_element_old),1) + &
                        W3 * bfy(HT(3,n_element_old),1) + W4 * bfy(HT(4,n_element_old),1)
            tbfield = W1 * bfz(HT(1,n_element_old),1) + W2 * bfz(HT(2,n_element_old),1) + &
                        W3 * bfz(HT(3,n_element_old),1) + W4 * bfz(HT(4,n_element_old),1)
            !$ Convert axisymmetric bfield to Cartesian bfield
            xbfield = rbfield*DCOS(PO%Z) - tbfield*DSIN(PO%Z)
            ybfield = rbfield*DSIN(PO%Z) + tbfield*DCOS(PO%Z)
          
            Omega(1) = zbfield * qm(isp+1) * 0.5 * TimeMove
            Omega(2) = xbfield * qm(isp+1) * 0.5 * TimeMove
            Omega(3) = ybfield * qm(isp+1) * 0.5 * TimeMove
            End If

            PO%Vx = PO%Vx + deta_vx
            PO%Vy = PO%Vy + deta_vy
            PO%Vz = PO%Vz + deta_vz

            f = 2.0/(1.0 + Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3))

            SVelocity(1) = (PO%Vx + PO%Vy*Omega(3) - PO%Vz*Omega(2)) * f
            SVelocity(2) = (PO%Vy + PO%Vz*Omega(1) - PO%Vx*Omega(3)) * f
            SVelocity(3) = (PO%Vz + PO%Vx*Omega(2) - PO%Vy*Omega(1)) * f
        
            PO%Vx = PO%Vx + SVelocity(2)*Omega(3) - SVelocity(3)*Omega(2) + deta_vx
            PO%Vy = PO%Vy + SVelocity(3)*Omega(1) - SVelocity(1)*Omega(3) + deta_vy
            PO%Vz = PO%Vz + SVelocity(1)*Omega(2) - SVelocity(2)*Omega(1) + deta_vz

        Else  !without magnetic field.
            If (delta_global == 0) Then  !2D Cartesian coordinates.
            deta_vx = qm(isp+1) * Efield(1) * TimeMove
            deta_vy = qm(isp+1) * Efield(2) * TimeMove
          
            PO%Vx = PO%Vx + deta_vx
            PO%Vy = PO%Vy + deta_vy
            Elseif (delta_global == 1) Then  !2D axis-symmetric coordinates.
            !$ axisymmetric e field
            zefield = Efield(1)
            refield = Efield(2)
            tefield = Efield(3)
            !$ Convert axisymmetric efield to Cartesian efield
            xefield = refield*DCOS(PO%Z) - tefield*DSIN(PO%Z)
            yefield = refield*DSIN(PO%Z) + tefield*DCOS(PO%Z)
          
            deta_vx = qm(isp+1) * zefield * TimeMove
            deta_vy = qm(isp+1) * xefield * TimeMove
            deta_vz = qm(isp+1) * yefield * TimeMove
          
            PO%Vx = PO%Vx + deta_vx
            PO%Vy = PO%Vy + deta_vy
            PO%Vz = PO%Vz + deta_vz
            End If
        End If

        IvelFlag = 0
    End If
    
    IF (delta == 0) THEN
          PO%X = PO%X + TimeMove*PO%Vx
		  PO%Y = PO%Y + TimeMove*PO%Vy
		  PO%Z = PO%Z + TimeMove*PO%Vz
    ELSEIF (delta == 1) THEN
      !>transform positions to cartesian coordinates
		  X=PO%Y*DCOS(PO%Z)
		  Y=PO%Y*DSIN(PO%Z)
		  Z=PO%X
      !>update cartesian positions
		  X=X+TimeMove*PO%Vy
		  IF(X == 0) X=hx(2)*1.0E-5
      Y=Y+TimeMove*PO%Vz
      PO%X = PO%X + TimeMove*PO%Vx
      !>update polar positions
      PO%Y = DSQRT(X*X+Y*Y)
      PO%Z = DATAN(Y/X)
      !>place particle in proper quadrant
      IF (X <= 0.0) THEN
        PO%Z=PO%Z+PI_LY
      ENDIF
    ENDIF
    IposFlag = 0
    
End Subroutine MoveOne
