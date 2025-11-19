SUBROUTINE OneAndChi_2D
   
USE Domain_2D 
USE Particle_2D
USE Field_2D
USE TimeControl
USE IMPIC_Data_2D
!=========LY modification, 2022-7-25=========
Use IFE_Data, Only:Field_Size
!=========LY modification, 2022-7-25=========
IMPLICIT NONE
INTEGER   :: ispe, i, j


!$ modified by ZWZ
!IF (.NOT.ALLOCATED(Chi)) THEN
!    ALLOCATE(Chi(0:nx+1,0:ny+1))
!ENDIF
!ALLOCATE(Chi(0:nx+1,0:ny+1))

!=========LY modification, 2022-7-25=========
Chi = 0.0

IF (Bfiled_index) THEN    !with Bfield

  TransChi = 0.0

  DO ispe = 1, ispe_tot
    Do i = 1, Field_Size
      Chi(i,1) = 0.5 * qm(ispe) * dt * dt * rho_s(i,1,ispe)   !0.5*q_i*rho_i*dt*dt/m_i
      
      Omega(1) = 0.5 * bfx(i,1) * qm(ispe) * dt
      Omega(2) = 0.5 * bfy(i,1) * qm(ispe) * dt
      Omega(3) = 0.5 * bfz(i,1) * qm(ispe) * dt
      OmegaT = 1.d0+Omega(1)*Omega(1)+Omega(2)*Omega(2)+Omega(3)*Omega(3)
      
      TransB(1,1)=1.d0+Omega(1)*Omega(1)
      TransB(1,2)=Omega(1)*Omega(2)+Omega(3)
      TransB(1,3)=Omega(1)*Omega(3)-Omega(2)
      TransB(2,1)=Omega(1)*Omega(2)-Omega(3)
      TransB(2,2)=1.d0+Omega(2)*Omega(2)
      TransB(2,3)=Omega(2)*Omega(3)+Omega(1)
      TransB(3,1)=Omega(1)*Omega(3)+Omega(2)
      TransB(3,2)=Omega(2)*Omega(3)-Omega(1)
      TransB(3,3)=1.d0+Omega(3)*Omega(3)
      
      TransB(1:3,1:3)=TransB(1:3,1:3)/OmegaT
      
      TransChi(1,i,1) = TransChi(1,i,1) + TransB(1,1) * Chi(i,1)
      TransChi(2,i,1) = TransChi(2,i,1) + TransB(1,2) * Chi(i,1)
      TransChi(3,i,1) = TransChi(3,i,1) + TransB(1,3) * Chi(i,1)
      TransChi(4,i,1) = TransChi(4,i,1) + TransB(2,1) * Chi(i,1)
      TransChi(5,i,1) = TransChi(5,i,1) + TransB(2,2) * Chi(i,1)
      TransChi(6,i,1) = TransChi(6,i,1) + TransB(2,3) * Chi(i,1)
      TransChi(7,i,1) = TransChi(7,i,1) + TransB(3,1) * Chi(i,1)
      TransChi(8,i,1) = TransChi(8,i,1) + TransB(3,2) * Chi(i,1)
      TransChi(9,i,1) = TransChi(9,i,1) + TransB(3,3) * Chi(i,1)
    End Do
  END DO

  !!!! 对于二维情况，实际上 TransChi 只有1，2，4，5 有用
  !!!! 注意：外面加的时单位矩阵，不是全1阵

  !!! ***** 2019-6-21 bjw 修改，认为在 Function_coefficient_2D.f90 里加Epsilon更合理 ******
  !TransChi(1,:,:) = TransChi(1,:,:) + 1.0
  !TransChi(2,:,:) = TransChi(2,:,:) + 0.0
  !TransChi(3,:,:) = TransChi(3,:,:) + 0.0
  !TransChi(4,:,:) = TransChi(4,:,:) + 0.0
  !TransChi(5,:,:) = TransChi(5,:,:) + 1.0
  !TransChi(6,:,:) = TransChi(6,:,:) + 0.0
  !TransChi(7,:,:) = TransChi(7,:,:) + 0.0
  !TransChi(8,:,:) = TransChi(8,:,:) + 0.0
  !TransChi(9,:,:) = TransChi(9,:,:) + 1.0
  !!! ***** 2019-6-21 bjw 修改，认为在 Function_coefficient_2D.f90 里加Epsilon更合理 ******


ELSE  !without Bfield

  DO ispe = 1, ispe_tot
    Do i = 1, Field_Size
      Chi(i,1) = Chi(i,1) + 0.5 * qm(ispe) * dt * dt * rho_s(i,1,ispe)  ! SUM(0.5*q_i*rho_i*dt*dt/m_i）
    End Do
  END DO

  !IF (.NOT.ALLOCATED(OneAndChi)) ALLOCATE(OneAndChi(0:nx+1,0:ny+1))
  !OneAndChi = 0.0
  !DO i = 1, nx
  !    DO j = 1, ny
  !        !OneAndChi(i,j) = 1.0 + Chi(i,j)
  !        !!! ***** 2019-6-21 bjw 修改，认为在 Function_coefficient_2D.f90 里加Epsilon更合理 ******
  !        OneAndChi(i,j) = 0.0 + Chi(i,j)
  !    END DO
  !END DO

END IF

!DEALLOCATE(Chi)
!WRITE(6,*) 'OneAndChi_2D END'


END