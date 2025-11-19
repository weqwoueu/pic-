SUBROUTINE OUTPUT_velocity(it)
    
use Particle_2D
Use ModuleMCCInterface,ONLY:ControlFlowGlobal, ParticleGlobal,JtoeV
Use Constant_Variable_2D
IMPLICIT NONE

        integer :: i, isp,it,N,num
        REAL    :: Vt0(2)
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: totalVelocity,driftVelocity
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: variance,thermalVelocity
        REAL, DIMENSION(:,:), ALLOCATABLE :: x_per,x_num
        CHARACTER(LEN=50) :: fname, filename

    IF (.NOT.ALLOCATED(totalVelocity)) THEN
        ALLOCATE(totalVelocity(ControlFlowGlobal%Ns+1,1024, 3))
        totalVelocity = 0.0
    END IF

    IF (.NOT.ALLOCATED(driftVelocity)) THEN
        ALLOCATE(driftVelocity(ControlFlowGlobal%Ns+1,1024, 3))
        driftVelocity = 0.0
    END IF

    IF (.NOT.ALLOCATED(variance)) THEN
        ALLOCATE(variance(ControlFlowGlobal%Ns+1,1024, 3))
        variance = 0.0
    END IF

    IF (.NOT.ALLOCATED(thermalVelocity)) THEN
        ALLOCATE(thermalVelocity(ControlFlowGlobal%Ns+1,1024, 3))
        thermalVelocity = 0.0
    END IF
    
    IF (.NOT.ALLOCATED(x_per)) THEN
        ALLOCATE(x_per(ControlFlowGlobal%Ns+1, 1025))
        x_per = 0.0
    END IF
    
    IF (.NOT.ALLOCATED(x_num)) THEN
        ALLOCATE(x_num(ControlFlowGlobal%Ns+1, 1024))
        x_num = 0.0
    END IF
    
         Vt0(1)=1.0
         Vt0(2)=0.01
         N=1024
    If (Mod(it,1000).eq.0 .Or. it == 1) Then
        DO isp = 0, ControlFlowGlobal%Ns
            DO i = 0, N
                x_per(isp+1, i+1) = REAL(i)
            END DO
        END DO
         
        DO isp = 0, ControlFlowGlobal%Ns
                DO i = 1,ParticleGlobal(isp)%Npar
                    num=MIN(MAX(1, FLOOR(ParticleGlobal(isp)%PO(i)%X)+1), N)
                            totalVelocity(isp+1,num,1) = totalVelocity(isp+1,num,1) + ParticleGlobal(isp)%PO(i)%Vx
                            totalVelocity(isp+1,num,2) = totalVelocity(isp+1,num,2) + ParticleGlobal(isp)%PO(i)%Vy
                            totalVelocity(isp+1,num,3) = totalVelocity(isp+1,num,3) + ParticleGlobal(isp)%PO(i)%Vz
                            x_num(isp+1,num)=x_num(isp+1,num)+1
                end do
        end do
        do isp=0,ControlFlowGlobal%Ns
            do num=1,N
                if(x_num(isp+1,num)==0)then 
                    driftVelocity(isp+1,num,1)=0
                    driftVelocity(isp+1,num,2)=0
                    driftVelocity(isp+1,num,3)=0
                else 
                    driftVelocity(isp+1,num,1)=totalVelocity(isp+1,num,1)/x_num(isp+1,num)
                    driftVelocity(isp+1,num,2)=totalVelocity(isp+1,num,2)/x_num(isp+1,num)
                    driftVelocity(isp+1,num,3)=totalVelocity(isp+1,num,3)/x_num(isp+1,num)
                end if
            end do
        end do
        
        

        DO isp = 0, ControlFlowGlobal%Ns
                DO i = 1, ParticleGlobal(isp)%Npar
                    num=MIN(MAX(1, FLOOR(ParticleGlobal(isp)%PO(i)%X)+1), N)
                    variance(isp+1,num,1) = variance(isp+1,num,1) + (ParticleGlobal(isp)%PO(i)%Vx - driftVelocity(isp+1,num,1))**2
                    variance(isp+1,num,2) = variance(isp+1,num,2) + (ParticleGlobal(isp)%PO(i)%Vy - driftVelocity(isp+1,num,2))**2
                    variance(isp+1,num,3) = variance(isp+1,num,3) + (ParticleGlobal(isp)%PO(i)%Vz - driftVelocity(isp+1,num,3))**2
                end do
        end do
        do isp=0,ControlFlowGlobal%Ns
          do num=1,N
              if(x_num(isp+1,num)==0)then
                  thermalVelocity(isp+1,num,1) = 0
                  thermalVelocity(isp+1,num,2) = 0
                  thermalVelocity(isp+1,num,3) = 0
              else
                  variance(isp+1,num,1)=variance(isp+1,num,1)/x_num(isp+1,num)
                  variance(isp+1,num,2)=variance(isp+1,num,2)/x_num(isp+1,num)
                  variance(isp+1,num,3)=variance(isp+1,num,3)/x_num(isp+1,num)
                  thermalVelocity(isp+1,num,1) = sqrt(variance(isp+1,num,1))
                  thermalVelocity(isp+1,num,2) = sqrt(variance(isp+1,num,2))
                  thermalVelocity(isp+1,num,3) = sqrt(variance(isp+1,num,3))
              end if
          end do
    end do
 
                !DO i = 1, ParticleGlobal(0)%Npar
                !    ParticleGlobal(0)%PO(i)%Vx=((ParticleGlobal(0)%PO(i)%Vx-driftVelocity(1,1))*(Vt0(1)/thermalVelocity(1,1))+driftVelocity(1,1))
                !    ParticleGlobal(0)%PO(i)%Vy=((ParticleGlobal(0)%PO(i)%Vy-driftVelocity(1,2))*(Vt0(1)/thermalVelocity(1,2))+driftVelocity(1,2))
                !end do

        
        !variance = 0.0
        !thermalVelocity=0.0
        !
        !DO isp = 0, ControlFlowGlobal%Ns
        !        DO i = 1, ParticleGlobal(isp)%Npar
        !            variance(isp+1,num,1) = variance(isp+1,num,1) + (ParticleGlobal(isp)%PO(i)%Vx - driftVelocity(isp+1,num,1))**2
        !            variance(isp+1,num,2) = variance(isp+1,num,2) + (ParticleGlobal(isp)%PO(i)%Vy - driftVelocity(isp+1,num,2))**2
        !            variance(isp+1,num,3) = variance(isp+1,num,3) + (ParticleGlobal(isp)%PO(i)%Vz - driftVelocity(isp+1,num,3))**2
        !        end do
        !end do
        !do isp=0,ControlFlowGlobal%Ns
        !  variance(isp+1,num,1)=variance(isp+1,num,1)/x_num(isp+1,num)
        !  variance(isp+1,num,2)=variance(isp+1,num,2)/x_num(isp+1,num)
        !  variance(isp+1,num,3)=variance(isp+1,num,3)/x_num(isp+1,num)
        !  thermalVelocity(isp+1,1) = sqrt(variance(isp+1,num,1))
        !  thermalVelocity(isp+1,2) = sqrt(variance(isp+1,num,2))
        !  thermalVelocity(isp+1,3) = sqrt(variance(isp+1,num,3))
        !end do
        !            
        
        WRITE(fname, '(I6.6)') it
        filename = './OUTPUT/Velocity/velocity_IJ_3' // TRIM(fname) // '.dat'
        OPEN(544, ACTION = 'WRITE', FILE = TRIM(filename)) 
        do num=1,N
            !WRITE(544, '(A150)') 'VARIABLES = "driftVelocity_ele" "thermalVelocity_ion" "thermalVelocity_ele"  "thermalVelocity_ion" '
            WRITE(544, 18) driftVelocity(1,num, 1), driftVelocity(2,num, 1), driftVelocity(3,num, 1), thermalVelocity(1,num, 1), thermalVelocity(2,num, 1), thermalVelocity(3,num, 1),x_num(1,num),x_num(2,num),x_num(3,num),num
        end do
    
        
18      FORMAT ((5X, F7.4), (5X, F7.4), (5X, F7.4),(5X, F8.4),,(5X, F8.4)(5X, F8.5),(4X,F9.1),(4X,F9.1),(4X,F9.1),(4X,I5))
        close(544)
    end if
    
        DEALLOCATE(totalVelocity, driftVelocity, variance, thermalVelocity,x_per,x_num)
    
end subroutine OUTPUT_velocity