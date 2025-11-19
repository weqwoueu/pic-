! This version is finished by 2016/12/09, the code can run well and it can be a produciable code, but still improvements are possible, please refer next version.
! Note Due to the crosssections are different from the ones used in our code, there are some differences in the results.    
    
Program PIC_MCC_for_CCP
   Use ModuleOneStep
   Use ModuleDiagOneStep
   Implicit none
   Integer(4) :: i,j,k,M
   real(8) Cpu1,Cpu2 
   
    !  Integer(4) :: NRun=1,NDiagShort=20,NDiagLong=0
   !Integer(4) :: NRun=10,NDiagShort=0,NDiagLong=0
   !Integer(4) :: NRun=10000,NDiagShort=200,NDiagLong=200
      
   Call AllInitilalization()
   Call OneStepRestart() 
   !Call CPU_TIME(CPU1)
   DO j=1,ControlFlowGlobal%NRun
       write(*,*) j
      do i=1,ControlFlowGlobal%Period    !t = j * i * (T/dt)     ControlFlowGlobal%Period = T/dt   T is the The time period of the change of the RF field
        Call OneStep()
		  ControlFlowGlobal%Timer = ControlFlowGlobal%Timer + 1
		  DO m=0,ControlFlowGlobal%Ns
          If (ParticleGlobal(m)%Npar>int(1.5*ParticleGlobal(m)%NParNormal)) then
            do k=0,ControlFlowGlobal%Ns     !if the particle is too many , use this to contract quantity of particle, 
                !Delect hulf of the particle and add it as weight to Remaining particle. N(real particle) = N * weight          
                Write(*,*) ParticleGlobal(k)%Npar,k,"before"
                 Call ParticleBundleNormalization(ParticleGlobal(k),INT(0.75*ParticleGlobal(k)%Npar))
                 Write(*,*) ParticleGlobal(k)%Npar,k,"after" 
           end do
			End If  
		 end do
     ENd DO
     open (10,position='append',file='ParticleNumber.dat')
     write (10,*) ParticleGlobal%NPar,ParticleGlobal%weight
     close (10)
   
     Call OneStepRestart()  !Temporary save ,one period or more
     Write(*,*) 'Period=',j,'   particle quanity=',ParticleGlobal%NPar,'   Total Timer:',ControlFlowGlobal%Timer
   ENd Do
   !Call CPU_TIME(CPU2)
   !Write(*,*) 'Period ',CPU2-CPU1,ParticleGlobal%NPar
   !Write(*,*) 'Period ', j,' Complete!'  
   Call DiagInitilalization(ControlFlowGlobal)
     do j=1,ControlFlowGlobal%NDiagShort
		   Write(*,*) 'Diag Period', j
         do i=1,ControlFlowGlobal%Period
 !            call Diag_Tracking_Particle(i,j,ControlFlowGlobal)
             Call OneStep()
             Call DiagOneStep(i,j)
         End do
     ENd Do
     Call DiagOneStepFinal()
   !

pause
stop
	end  Program



	
	
