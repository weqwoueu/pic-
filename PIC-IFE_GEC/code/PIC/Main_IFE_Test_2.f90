!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name:Main_IFE_Test_2.f90                                 C
!
!  Purpose: Main Program
!                                                                      C
!  Reviewer: Yuchuan Chu                              Date: 03-May-12  C
!  Comments: Normalization finished in code                            C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

USE IMPIC_Data_2D
USE PIC_Main_Param_2D
USE Domain_2D
USE Field_2D
USE Particle_2D
USE TimeControl
USE Object_Data_2D
USE IFE_Boundary
USE IFE_Data
USE Wall_2D         !? this module is not used, delete?
USE IFE_Interface_cell_volume, ONLY: Setup_Cell_Volume_2D   !? the Setup_Cell_Volume_2D is not used, delete the mudule?
USE IFE_INTERFACE, ONLY: IFE_START_2D, Input_2D, Setup_IFE_Mesh_2D, &
                        Setup_IFE_Wall_Mesh_2D_NEW, AdjustBoundary_Sputter_2D, &                !? the Setup_IFE_Wall_Mesh_2D_NEW and AdjustBoundary_Sputter_2D are not used, delete?
                        Huygens_Wave_2D,Output_To_Tecplot_Wall_2D, Update_IFE_Start_2D, &       !? the Huygens_Wave_2D and Update_IFE_Start_2D are not used, delete?
                        AdjustBoundary_2D, CheckTrail, InjectBeams_2D   !$ ab.ZWZ 2021/7/9
USE ModuleMCCInterface   !$ ab.ZWZ 2021/12/19 for JW's MCC
Use ModuleDiagOneStep
IMPLICIT NONE

!? ------------------- what are they used for? -------------------------
INTEGER, PARAMETER	::	nbkgd = 1	    
REAL(8)					IFE_phi_bkgd(nbkgd), IFE_Te_bkgd(nbkgd), IFE_Rho_bkgd(nbkgd)
!? ------------------- what are they used for? --------------------------

INTEGER		it, jj,i,j,i_part   !? jj and i_part arer not used in Main, delete ?
REAL(8)		xt

REAL(8)		R_1     !? R_1 is not used in Main, delete ?

! delta  0: 2-D Cartesian coordinates, 1: cylindrical coordinate, axisymmetric situation
INTEGER		::  delta 
INTEGER		::	Ndisf,Update        !? The Ndisf is used in SetupPartInject_QLL, delete? ; What is "Update" used for?

!? the follow varibles are not used in Main, delete?
!	Output
INTEGER ni,nj, nindex
REAL(8), DIMENSION(:), POINTER ::	R_full, U_full
!? the above varibles are not used in Main, delete?

!   Input
REAL(8)		xmin, xmax, ymin, ymax, zmin, zmax   !, PIC_Phi_BC(6)
INTEGER		nnx, nny, N_Objects, N_Boundary
!INTEGER		nnx, nny, N_Boundary
INTEGER     nnz
TYPE(ObjectType), DIMENSION(:), POINTER		::	objects

!? the follow varibles are not used in Main, delete?
INTEGER   N_theta, ii
!REAL(8)  :: R, rxp, ryp, dx, dy, f, g
REAL(8)  ::  rxp, ryp, dx, dy, f, g !$ mb.ZWZ 2021/12/19 for JW's MCC
REAL(8), PARAMETER	::	pii	= 3.14159265358979D0
!REAL(8),DIMENSION(:,:),ALLOCATABLE     :: P1
!? the above varibles are not used in Main, delete?

REAL(8) :: xp, yp, resu, error  !$ ab.ZWZ to check convergence 2021/7/9     !? used in convergence checking, delete ?
Logical :: DumpFlag

Integer(4) :: isp 


CALL My_Label(  '2 Dimension I F E - P I C', 'Immersed Finite Element Particle-In-Cell Code', &
				'Y. Cao and R. Kafafy', 'Dr. J. Wang and Dr. T. Lin')
				
CALL Input_2D(xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects, bc_point_1, bc_point_2, bc_index, bc_value)


!---- Setup IFE Wall Mesh
IF(N_Objects/=0)THEN
    !CALL Setup_IFE_Wall_Mesh_2D_NEW(objects)   !? this subroutine is not used, delete ?
ENDIF

!---- Setup IFE Mesh
CALL Setup_IFE_Mesh_2D(delta, xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, objects)

!$ ab.ZWZ 2021/7/9 for passing value to Domain module
delta_global = delta 
region_type=0      !0--初始分布为矩形，1--扇形
dxmin = xmin
dxminmin = 0        !zzj for initialization
dxmax = xmax
dxmaxmax = 200.0
dymin = ymin
dyminmin = 0
dymax = ymax
dymaxmax = 4.0
Radius_L = 0
Radius_U = 8
Theta_L = 0
Theta_U = pii/2

!---- Setup Cell Volume
!CALL Setup_Cell_Volume_2D(delta, xmin, xmax, ymin, ymax, nnx, nny, N_Objects, N_Boundary, Objects) !? this subroutine is not used, delete ?

!---- Initialization and Start 
CALL Start_2D

!	SETUP GEOMETRY AND BOUNDARY CONDITIONS
!---- Setup Domain
CALL SetupGrids_2D_QLL(delta, xmin, xmax, ymin, ymax, zmin, zmax, nnx, nny, nnz)

!---- Setup Outer Boundary
CALL SetBoundaryType_2D_QLL

!---- Setup Outer Boundary Phi BC 
!CALL SetPhiOuterBC_2D           ! useless in SIDG

!---- Setup PIC Boundary
!CALL Set_PIC_Boundary_2D(nnx, nny,N_Objects)
!CALL Set_PIC_Boundary_2D(nnx, nny)    ! useless in SIDG

!---- Setup Objects
IF(N_Objects/=0)THEN
    CALL SetObjects_2D ! should revise when test
ENDIF
    
    
!---- Set the Electron Density & Temperature Reference
!CALL SetEleRef_2D   ! useless in SIDG

!!!! bjw add 2019.10.23
!CALL SetupBfield_circle(xmin, xmax, ymin, ymax, zmin, zmax, nnx, nny, nnz, VertX)  !? can this subroutine be delete?

!$ ab.ZWZ 2021/9/8 for LH's magnetic field
IF (Bfiled_Index) Then
    !CALL SetMagFld_2D(delta) ! should revise when use SIDG
Endif
    
!---- Set the flux boundary
!CALL FluxBoundary_2D       !? cannot find this subroutine, delete?

!---- Setup Load Particles ****
CALL SetupPartInject_QLL(Ndisf)     !? used in DSMC, retain?



!IF(ispe_tot>2 )THEN
!    CALL Init_DSMC_2D  !? DSMC, retain?
!ENDIF

!CALL Init_MCC_2D       !? the original MCC(not JW's) init, delete?

!$ ========= ab.ZWZ 2021/12/19 for JW's MCC =========== \\
CALL AllInitialization() 
Call DiagInitilalization(ControlFlowGlobal)
!$ ========= ab.ZWZ 2021/12/19 for JW's MCC =========== //

!---- Comment this IF you wanna start from scratch
CALL Restart_2D  !$ Should be revise When use in SIDG wsy 2022 7 22


!---- Initialize Phi
IF (irestart/=1) THEN
	CALL InitialPhi_2D  !? what is this subroutine used for ?
END IF

xt = 0.

!---------------
!	Initialize IFE Solver
CALL IFE_Start_2D(	IFE_phi_bkgd, IFE_Te_bkgd, IFE_Rho_bkgd, nbkgd,		&
				    Phi, nx, ny, delta, xt		)


 CALL Get_PointPatch(nnx, nny, xmin, xmax, ymin, ymax)
 
 
!? the follow code is used for convergence checking and particle trail checking, delete?
!$ ====================== ab.ZWZ to check convergence ========================== \\
!CALL IFE_Solve_2D(Rho, Phi, nx, ny, delta, xt)
!
!CALL Get_PointPatch(nnx, nny, xmin, xmax, ymin, ymax)
!CALL GetEField_SIDG_PPR
!
!CALL Output_To_Tecplot_IJK_2D(	0, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
!                                VertX, phi, rho, rho_s, efx, efy, part, 1, P_average, HP, HT, SIZE(HP,2), SIZE(HT,2), SIZE(P_average,2))

!=========IFE=========
!CALL Error_Analysis_2D(delta, nx)
!=========IFE=========

!OPEN(1, ACTION = 'WRITE', FILE = './error_field.dat')
!WRITE(1,*) 'TITLE = "Field Plot"'
!WRITE(1,"(A50)") 'VARIABLES = "x" "y" "Phi" "Phi_Real" "error"'!$ modified by ZWZ  2020/10/28
!WRITE(1,"(' ZONE I = ',I6,', J= ',I6)") nnx, nny
!	DO j=1,nny
!		DO i=1,nnx
!            nindex = j+(i-1)*nny
!            xp = p_basic(1, nindex)
!            yp = p_basic(2, nindex)
!            CALL Function_True(delta, 10._8, xp, yp, 0, 0, resu)
!            error = resu - phi(i,j)
!			WRITE(1,"(E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6)") VertX(1:2,i,j), Phi(i,j), resu, error
!		END DO
!	END DO
!CLOSE(1)
!$ ====================== ab.ZWZ to check convergence ========================== //
!$ ====================== ab.ZWZ to check trail of particle ==================== \\
!CALL CheckTrail(N_Objects,objects)
!$ ====================== ab.ZWZ to check trail of particle ==================== //


!!******Solve the Initial Field *********** 
xt = 0.

IF (irestart/=1) THEN
    
    !CALL Load_Charged_Particle(delta, xmin, xmax, ymin, ymax)
    CALL AdjustBoundary_2D(xt,dt,delta,N_Objects,objects)	

!>------ Get particle charge   -------------------------------- 
	CALL GetParChg_2D(delta) ! shoule be revise in SIDG

!>------- Solve the PoissonEq using IFE
	CALL IFE_Solve_2D(Rho, Phi, nx, ny, delta, xt)
    
    
    !Do i = 1, Size(HP,2)
    !Phi(i, 1) = 773*(1 - ((HP(1,i))/(maxval(HP(1,:))))**6)
    !!Phi(i, 1) = 0
    !End Do
    
	!CALL GetEfield_2D
   Call GetEField_SIDG_PPR
    
    CALL Output_To_Tecplot_IJK_2D(	0, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
								    VertX, phi, rho, rho_s, efx, efy, part, 1, P_average, HP, HT, SIZE(HP,2), SIZE(HT,2), SIZE(P_average,2))
    !pause
    CALL Dump_2D(0) !$ ab.ZWZ
    
    do isp=0,ControlFlowGlobal%Ns
        Call ParticleGlobal(isp)%Dump(0)
    End do
    
   !$ ================================= mb.ZWZ for writing title for output file =========================================\\
    OPEN(540,ACTION='WRITE', FILE='./OUTPUT/PartcountReal.dat',POSITION='APPEND', STATUS='REPLACE')
        WRITE(540,"(A100)")'VARIABLES = "it" "tot" "num_e" "num_i"'
    CLOSE(540)
    
    Open(541,Action='Write', File='./OUTPUT/ElectronChange.dat',Position='Append',Status='Replace')
        Write(541,"(A400)")'Variables = "it" "Change_ele"&
                                        "Gen_tot_ele" "Gen_BL_ele" "Gen_BR_ele" "Gen_MCC_ele" &
                                        "Loss_tot_ele" "Loss_B1_ele" "Loss_B2_ele" "Loss_B3_ele" "Loss_B4_ele" "Loss_MCC_ele"'
    Close(541)
    
    Open(542,Action='Write', File='./OUTPUT/IonChange.dat',Position='Append',Status='Replace')
        Write(542,"(A400)")'Variables = "it" "Change_ion"&
                                        "Gen_tot_ion" "Gen_BL_ion" "Gen_BR_ion" "Gen_MCC_ion" &
                                        "Loss_tot_ion" "Loss_B1_ion" "Loss_B2_ion" "Loss_B3_ion" "Loss_B4_ion" "Loss_MCC_ion"'
    Close(542)
    
    OPEN(540,FILE='./OUTPUT/PartcountReal.dat',POSITION='APPEND')
        WRITE(540,*) 0, ntot * ParticleGlobal(0)%Weight, (ns(isp+1)* ParticleGlobal(isp)%Weight,isp=0,ControlFlowGlobal%Ns)
    CLOSE(540)
    !$ ============================================ mb.ZWZ ================================================================//     
    
ELSE

!------ get efield --------------------
  Call GetEField_SIDG_PPR   !LY modification, 2022-7-25
  
!***start
	xt=ilap*dt
    
END IF

Call SYSTEM_CLOCK(Time1)
!!!**********begin the time loop******************************* 
DO it = ilap+1, nt
	xt=xt + dt
	WRITE(6,*) 
    ntot = ntot/ParticleGlobal(0)%Weight
	WRITE(6,*) ' ***** it=',it, '  ntot=', ntot, '   nPeriod=',Real(it)/Real(PStep),'*****'
    
    Do isp=0,ControlFlowGlobal%Ns
        ntot = ntot + ParticleGlobal(isp)%Npar 
        ns(isp+1) = ParticleGlobal(isp)%Npar 
    End Do
    
    call OUTPUT_velocity(it)
    call Output_Energy(it)

!------ Inject new slice of Particles from upstream
    !> inject particles from boundary
 !   write(*,*) 'inject'
	!CALL InjectBeams_2D(it, delta, ParticleGlobal)
  
    !> supplement electron form boundary to maintain quasi state
    !CALL Quasi_2D(delta,it)
    !Call Quasi_cdp(ParticleGlobal(0), ParticleGlobal(1), delta, it)
    
    ! =============== ab.ZWZ for adding JW's MCC code ============== \\
    !Call MCCEntrance
    ! =============== ab.ZWZ for adding JW's MCC code ============== //
    !CALL GetParChg_2D(delta)
    !CALL Output_To_Tecplot_IJK_2D(it, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
    !                        VertX, phi, rho, rho_s, efx, efy, part, 1, P_average, HP, HT, SIZE(HP,2), SIZE(HT,2), SIZE(P_average,2))
    !pause
    
    IF (IMPIC_index) THEN   !!!! 隐式
       
      Do i=0,ControlFlowGlobal%Ns
        Call MoveAll(ParticleGlobal(i), ControlFlowGlobal,N_objects,objects,delta,1)
      End Do

      CALL GetParChg_2D(delta)

      CALL OneAndChi_2D
      CALL IFE_Start_2D(	IFE_phi_bkgd, IFE_Te_bkgd, IFE_Rho_bkgd, nbkgd,		&
                          Phi, nx, ny, delta, xt	)

      CALL IFE_Solve_2D(Rho, Phi, nx, ny, delta, xt)
      
      Call GetEField_SIDG_PPR

      Do i=0,ControlFlowGlobal%Ns
        Call MoveAll(ParticleGlobal(i), ControlFlowGlobal,N_objects,objects,delta,2)
      End Do
        
    ELSE  !!!! 显式
         !write(*,*) 'move'
      Do i=0,ControlFlowGlobal%Ns
        Call MoveAll(ParticleGlobal(i), ControlFlowGlobal,N_objects,objects,delta,0)
      End Do
         !write(*,*) 'gpt'
      CALL GetParChg_2D(delta)
      
      CALL IFE_Solve_2D(Rho, Phi, nx, ny, delta, xt)
            
      Call GetEField_SIDG_PPR
        
    END IF    
    
    !write(*,*) 'inject'
	!CALL InjectBeams_2D(it, delta, ParticleGlobal)
    
    ! ========================OUTPUT================================
    
    
    ! === Average value accumulating === \\
    !Do isp = 1,ispe_tot
    !  !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    !  Do i = 1, Field_Size
    !    If (isp == 1) Then
    !      Aver_phi(i,1) = Aver_phi(i,1) + Phi(i,1)
    !      Aver_rho(i,1) = Aver_rho(i,1) + Rho(i,1)
    !    End If
    !    Aver_rho_s(i,1,isp) = Aver_rho_s(i,1,isp) + rho_s(i,1,isp)
    !    Aver_Ek_s(i,1,isp) = Aver_Ek_s(i,1,isp) + Ek_s(i,1,isp)
    !    Aver_Ek_tot(i,1,isp) = Aver_Ek_tot(i,1,isp) + Ek_tot(i,1,isp)
    !  End Do
    !  !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
    !End Do
    
    !Call DiagOneStep(ControlFlowGlobal) !> get the EPF and accumulate it
    !$ === Average value accumulating === //
    
    !> output average value every period and zero out them
    IF(it - iStep_previous == Pstep .And. .False. ) THEN
        nPeriod = nPeriod + 1
        
        !> Averaging Grid Value data
        AverStep = AverStep + Pstep
        If (Mod(nPeriod,32).eq.0) Then
            !> Output EPF
            !Call DiagOneStepFinal(ControlFlowGlobal)    !> output the averaged EPF
            
            !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
            Do i = 1, Field_Size
              Aver_phi(i,1) = Aver_phi(i,1)/AverStep
              Aver_rho(i,1) = Aver_rho(i,1)/AverStep
              Do isp = 1,ispe_tot
                Aver_rho_s(i,1,isp) = Aver_rho_s(i,1,isp)/AverStep
                Aver_Ek_s(i,1,isp) = Aver_Ek_s(i,1,isp)/AverStep
                Aver_Ek_tot(i,1,isp) = Aver_Ek_tot(i,1,isp)/AverStep
              End Do
            End Do
            !=========LY modification for Multi-Layer-Grid, 2022-7-25=========        
        
            !> Output real time data
            CALL Output_To_Tecplot_IJK_2D(it, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
                            VertX, phi, rho, rho_s, efx, efy, part, 3, P_average, HP, HT, SIZE(HP,2), SIZE(HT,2), SIZE(P_average,2))
            !> Output average field data
            !CALL Output_To_Tecplot_IJK_2D(it, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
            !                VertX, phi, rho, rho_s, efx, efy, part, 3, P_average, HP, HT, SIZE(HP,2), SIZE(HT,2), SIZE(P_average,2))
        
            !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
            Do i = 1, Field_Size
              Aver_phi(i,1) = 0.0
              Aver_rho(i,1) = 0.0
              Do isp = 1,ispe_tot
                Aver_rho_s(i,1,isp) = 0.0
                Aver_Ek_s(i,1,isp) = 0.0
                Aver_Ek_tot(i,1,isp) = 0.0
              End Do
            End Do
            !=========LY modification for Multi-Layer-Grid, 2022-7-25=========
            AverStep = 0
        Endif
        
        iStep_previous = it
        Pstep = INT(Period + Period_res)
        Period_res = Period + Period_res - Pstep
        
        !If (Mod(nPeriod,32) .eq. 0) Then
        !    DumpFlag = .True.
        !EndIf
        
        
        
        
        If (Mod(nPeriod,1) .eq. 0) Then
            Do i = 0,ControlFlowGlobal%Ns
                If ( ParticleGlobal(i)%UnequalWeightFlag) Then
                    !ParticleGlobal(i)%WQOne = 2 * ParticleGlobal(i)%WQOne
                    Write(*,*)  ParticleGlobal(i)%NPar,'NParBefore'    
                    Call WeightRenormalization(ParticleGlobal(i)) 
                    Write(*,*)  ParticleGlobal(i)%NPar,'NParAfter'
                    Call Coalescing(ParticleGlobal(i),nx-1,ny-1)
                    Write(*,*)  ParticleGlobal(i)%NPar,'NParAfter2'
                Endif
            Enddo
        !        
            !> for checking effect of particle coalescing
            !CALL GetParChg_2D(delta)
            !!> Output average field data
            !CALL Output_To_Tecplot_IJK_2D(	it+1, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
								    !VertX, phi, rho, rho_s, efx, efy, part, 1)
        Endif
        
        
    END IF
    
    !> If the number of simulate particle is too much, 
    do isp=0,ControlFlowGlobal%Ns
        If (ParticleGlobal(isp)%Npar>ParticleGlobal(isp)%NParNormal) then
            Do i = 0,ControlFlowGlobal%Ns
                Write(*,*) ParticleGlobal(i)%Npar,i,"before"
                If (ParticleGlobal(i)%UnequalWeightFlag) Then   
                    !Call ParticleBundleNormalizationUnequalWeight(ParticleGlobal(i),ParticleGlobal(i)%Npar/2,ny-1)
                    ParticleGlobal(i)%WQOne = 2 * ParticleGlobal(i)%WQOne
                    OPEN(1,ACTION='WRITE', FILE='./OUTPUT/Coalescing.dat',POSITION='APPEND')
                        Write(1,*)  ParticleGlobal(i)%NPar,'NParBefore' ,it
                        Call Coalescing(ParticleGlobal(i),nx-1,ny-1)
                        Write(1,*)  ParticleGlobal(i)%NPar,'NParAfter',it
                        Call WeightRenormalization(ParticleGlobal(i)) 
                        Write(1,*)  ParticleGlobal(i)%NPar,'NParAfter2',it
                    Close(1)
                Else    !> Increase weight of particles and randomly delete certain fraction (EqualWeight)
                    Call ParticleBundleNormalization(ParticleGlobal(i),ParticleGlobal(i)%Npar/2)
                Endif
                Write(*,*) ParticleGlobal(i)%Npar,i,"after" 
            End Do
        End If   
    end do

    If (Mod(it,1000).eq.0 .Or. it == 1) Then
    !If (it - iStep_previous == Int(Pstep/4) ) Then
        CALL Output_To_Tecplot_IJK_2D(it, nx, ny, ispe_tot, ntot, SIZE(part,1), SIZE(part,2), n_stride, &
                            VertX, phi, rho, rho_s, efx, efy, part, 1, P_average, HP, HT, SIZE(HP,2), SIZE(HT,2), SIZE(P_average,2))
    Endif
	
    !$ ============================ mb.ZWZ ==================\\
    ntot = 0
    ns = 0
    Do isp=0,ControlFlowGlobal%Ns
        If (ParticleGlobal(isp)%UnequalWeightFlag) Then
            Do i = 1,ParticleGlobal(isp)%Npar
                ntot = ntot + ParticleGlobal(isp)%PO(i)%WQ
                ns(isp+1) = ns(isp+1) + ParticleGlobal(isp)%PO(i)%WQ
            End Do
        Else
            ntot = ntot + ParticleGlobal(isp)%Npar * ParticleGlobal(isp)%Weight
            ns(isp+1) = ParticleGlobal(isp)%Npar *  ParticleGlobal(isp)%Weight
        Endif
    End Do
    
    ! !$ ============================ zzj ==================\\
    !IF (.NOT.ALLOCATED(Ek_one))	ALLOCATE(Ek_one(0:ControlFlowGlobal%Ns,1:ParticleGlobal(isp)%Npar))
    !Do isp=0,ControlFlowGlobal%Ns
    !    Do i=1,ParticleGlobal(isp)%Npar
    !        Ek_one(isp,i) = ParticleGlobal(isp)%PO(i)%Energy(ParticleGlobal(isp)%Mass,ParticleGlobal(isp)%VFactor)
    !    end do
    !end do
    !
    
    
    
    OPEN(540,FILE='./OUTPUT/PartcountReal.dat',POSITION='APPEND')
        WRITE(540,*) it, ntot, (ns(isp+1),isp=0,ControlFlowGlobal%Ns)
    CLOSE(540)
    !If (ParticleGlobal(0)%UnequalWeightFlag) Then
    !    Do isp=0,ControlFlowGlobal%Ns
    !        Do i = 1,ParticleGlobal(isp)%Npar
    !            ntot = ntot + ParticleGlobal(isp)%PO(i)%WQ
    !            ns(isp+1) = ns(isp+1) + ParticleGlobal(isp)%PO(i)%WQ
    !        End Do
    !    End Do
    !    OPEN(540,FILE='./OUTPUT/PartcountReal.dat',POSITION='APPEND')
    !        WRITE(540,*) it, ntot, (ns(isp+1),isp=0,ControlFlowGlobal%Ns)
    !    CLOSE(540)
    !Else
    !    Do isp=0,ControlFlowGlobal%Ns
    !        ntot = ntot + ParticleGlobal(isp)%Npar
    !        ns(isp+1) = ParticleGlobal(isp)%Npar
    !    End Do
    !    OPEN(540,FILE='./OUTPUT/PartcountReal.dat',POSITION='APPEND')
    !        WRITE(540,*) it, ntot * ParticleGlobal(0)%Weight, (ns(isp+1)* ParticleGlobal(isp)%Weight,isp=0,ControlFlowGlobal%Ns)
    !    CLOSE(540)
    !Endif
    
     
    
    !>> ab.ZWZ for checking particle generation and loss for debuging
    !Do isp=0,ControlFlowGlobal%Ns
    !    ParticleGlobal(isp)%nGen(0)=Sum(ParticleGlobal(isp)%nGen(1:3))
    !    ParticleGlobal(isp)%nLoss(0)=Sum(ParticleGlobal(isp)%nLoss(1:5))
    !End Do
    !
    !OPEN(541,FILE='./OUTPUT/ElectronChange.dat',POSITION='APPEND')
    !    WRITE(541,*) it, (ParticleGlobal(0)%nGen(0)-ParticleGlobal(0)%nLoss(0))*ParticleGlobal(0)%Weight,&
    !                 ParticleGlobal(0)%nGen(0:3)*ParticleGlobal(0)%Weight, ParticleGlobal(0)%nLoss(0:5)*ParticleGlobal(0)%Weight
    !CLOSE(541)
    !
    !OPEN(542,FILE='./OUTPUT/IonChange.dat',POSITION='APPEND')
    !    WRITE(542,*) it, (ParticleGlobal(1)%nGen(0)-ParticleGlobal(1)%nLoss(0))*ParticleGlobal(1)%Weight,&
    !                 ParticleGlobal(1)%nGen(0:3)* ParticleGlobal(1)%Weight, ParticleGlobal(1)%nLoss(0:5)* ParticleGlobal(1)%Weight
    !CLOSE(542)
    !
    !Do isp=0,ControlFlowGlobal%Ns
    !    ParticleGlobal(isp)%nGen=0
    !    ParticleGlobal(isp)%nLoss=0
    !End Do
    !$ ============================ mb.ZWZ ==================//
    
    
    !IF(MOD(it,n_dump).eq.0) CALL Dump_2D(it)
    If (Mod(it,1000) .eq. 0) Then
            DumpFlag = .True.
        EndIf
	IF(DumpFlag) Then
        CALL Dump_2D(it)
        do isp=0,ControlFlowGlobal%Ns
            Call ParticleGlobal(isp)%Dump(0)
        End do
        DumpFlag = .False.
    End If
END DO

WRITE(6,*) 
WRITE(6,*) 'run finish at after it=', it-1 

CLOSE(540)  !$ close partcount file

     
END

