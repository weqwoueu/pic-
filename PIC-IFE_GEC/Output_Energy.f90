SUBROUTINE Output_Energy(it)
    
use Particle_2D
Use ModuleMCCInterface,ONLY:ControlFlowGlobal, ParticleGlobal,JtoeV
Use Constant_Variable_2D
IMPLICIT NONE

    INTEGER    :: N  ! 能量区间数量
    INTEGER    :: i, it, isp, i_part,n_ini
    REAL(8)    :: Ek_max,k_energy,R,R0,Vx,Vy,Vz
    REAL(8)    :: Ek_one,Vk_one
    REAL(8)    :: Vk_max(3),k_v(3)
    REAL(8)    :: X_Emin, X_Emax, Y_Emin, Y_Emax, k_area
    REAL(8), PARAMETER	::	pii	= 3.14159265358979D0
    CHARACTER(LEN=50) :: fname, filename
    REAL(8)    :: dN_dE(ControlFlowGlobal%Ns+1, 100), dN_dV(ControlFlowGlobal%Ns+1, 100)
    REAL(8)    :: Ek_per(ControlFlowGlobal%Ns+1, 101), Vk_per(ControlFlowGlobal%Ns+1, 101)
    REAL(8)    :: N_per(ControlFlowGlobal%Ns+1, 101), NV_per(ControlFlowGlobal%Ns+1, 101)

    Ek_max = 5.0
    Vk_max(1) = 10
    Vk_max(2) = 0.1
    Vk_max(3) = 0.01
    N = 100
    R0=8.0
    !X_Emin=4.0      !目标区域
    !X_Emax=8.0
    !Y_Emin=4.0
    !Y_Emax=8.0
    k_energy=10.0
    k_v(1)=20.0
    k_v(2)=1000
    k_v(3)=10000
    !k_area=4.0
    n_ini=4000*4*128   !MCC中粒子给的多少 初始密度




    ! 初始化能量区间
    DO isp = 0, ControlFlowGlobal%Ns
        DO i = 0, N
            Ek_per(isp+1, i+1) = Ek_max * REAL(i) / N
            Vk_per(isp+1, i+1) = Vk_max(isp+1) * REAL(i-(0.5*N)) / N
            N_per(isp+1, i+1) = 0._8
            NV_per(isp+1, i+1) = 0._8
        END DO
    END DO
    R=0._8
    
       
    ! 统计每个能量区间中的粒子数

    If (Mod(it,2000).eq.0 .Or. it == 1) Then
        DO isp = 0, ControlFlowGlobal%Ns
            DO i = 1, ParticleGlobal(isp)%Npar
                Ek_one = ParticleGlobal(isp)%PO(i)%Energy(ParticleGlobal(isp)%Mass, ParticleGlobal(isp)%VFactor)/ (JtoeV)
                Vx=ParticleGlobal(isp)%PO(i)%Vx
                Vy=ParticleGlobal(isp)%PO(i)%Vy
                Vz=ParticleGlobal(isp)%PO(i)%Vz
                Vk_one=Vx
                    DO i_part = 1, N
                        IF (Ek_one <= Ek_per(isp+1, i_part+1) .AND. Ek_one > Ek_per(isp+1, i_part)) THEN
                            N_per(isp+1, i_part) = N_per(isp+1, i_part) + 1
                        END IF
                        IF (Vk_one<= Vk_per(isp+1, i_part+1) .AND. Vk_one > Vk_per(isp+1, i_part)) THEN
                            NV_per(isp+1, i_part) = NV_per(isp+1, i_part) + 1
                        END IF
                    END DO
            END DO
        END DO


    !! 计算能量区间的分布函数 f(E)
    DO isp = 0, ControlFlowGlobal%Ns


        ! 计算分布函数的导数 dN/dE
        DO i = 1, N
            dN_dE(isp+1,i)=0._8
            dN_dE(isp+1,i) = k_energy*N_per(isp+1,i)
            dN_dE(isp+1,i)=dN_dE(isp+1,i)/n_ini
            dN_dV(isp+1,i)=0._8
            dN_dV(isp+1,i) = k_V(isp+1)*NV_per(isp+1,i)
            dN_dV(isp+1,i)=dN_dV(isp+1,i)/n_ini
        END DO

    END DO
    
    ! 准备文件名并打开文件
    WRITE(fname, '(I6.6)') it
    filename = './OUTPUT/Energy/energy_IJ2_' // TRIM(fname) // '.dat'
    OPEN(543, ACTION = 'WRITE', FILE = TRIM(filename))
    
    ! 写入文件头
    WRITE(543, *) 'TITLE = "EnergySpectrum Plot"'
    WRITE(543, '(A150)') 'VARIABLES = "dN/dE_ele" "dN/dE_ion" "E" "V1" "V2" "dN/dV_ele" "dN/dV_ion" '
    
    ! 写入数据
    DO i = 1, N
        WRITE(543, 17) &
            dN_dE(1, i), dN_dE(2, i),dN_dE(3, i), Ek_per(1, i), Vk_per(1, i), Vk_per(2, i), Vk_per(3, i), dN_dV(1, i), dN_dV(2, i),dN_dV(3, i)
        WRITE(543, *)
    END DO
    
    17 FORMAT ((5X, F9.3), (5X, F9.3), (5X, F9.3),(5X, F8.3),(5X, F9.5),(5X, F9.5),(5X, F8.4), (5X, F8.4))
    
    CLOSE(543)

    
    endif
    


END SUBROUTINE Output_Energy





    

        