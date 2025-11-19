SUBROUTINE Output_To_Tecplot_IJK_2D(	it, nnx, nny, ispe_max, ntot, size_part1, size_part2, n_stride,		&
										VertX, Phi, Rho, Rho_s, efx, efy, part, Output_Option, P_average, HP, HT, Node_Size, Element_Size, Aver_Size)

!USE Wall_2D
USE Field_2D, only : Aver_phi, Aver_rho, Aver_rho_s, Aver_Ek_s, Aver_Ek_tot, Ek_s, Ek_tot
USE IFE_Data, only : node_map, node_type
USE Constant_Variable_2D
USE TimeControl, Only : nPeriod
IMPLICIT NONE

INTEGER :: it, nnx, nny, size_part1, size_part2,n_sect
INTEGER :: ispe_max
INTEGER :: n_stride(1:ispe_max)	

REAL(8) :: VertX(2,0:nnx+1,0:nny+1), ntot
REAL(8) :: part(1:size_part1,1:size_part2)

INTEGER :: Output_Option

INTEGER :: i, j, isp, temp_count, NumOfNodeX,i_part,j_part
CHARACTER*50	fname, sname, filename, num_n, num_e

REAL(8),PARAMETER :: SmallValue = 1.0D-5 
INTEGER :: count, n_node, n_element, count_grid, Node_Size, Element_Size, Aver_Size,n_node_part
REAL(8) :: h_partition, h_flag, radius_max

!LY REVISE, 2022-1-15
!REAL(8),DIMENSION(:,:),POINTER :: HP, P_average
!INTEGER,DIMENSION(:,:),POINTER :: HT
REAL(8) :: HP(2,Node_Size), P_average(6,Aver_Size)
INTEGER :: HT(5,Element_Size)
REAL(8),DIMENSION(:),ALLOCATABLE :: Phi_aver_axis, Rho_aver_axis, Rho_ele_aver_axis, Rho_ion_aver_axis
REAL(8),DIMENSION(:),ALLOCATABLE :: Ek_s_ele_aver_axis, Ek_s_ion_aver_axis, Ek_tot_ele_aver_axis, Ek_tot_ion_aver_axis
REAL(8),DIMENSION(:),ALLOCATABLE :: radius
REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: Rho_temp

REAL(8) :: Rho(Node_Size, 1), Rho_s(Node_Size, 1, 1:ispe_max)
REAL(8) :: Phi(Node_Size, 1)
REAL(8) :: efx(Node_Size, 1), efy(Node_Size, 1)

ALLOCATE(radius(SIZE(HP,2)))
ALLOCATE(Rho_temp(SIZE(HP,2), 1, 2))
!PRINT*, "Saving Data in Tecplot IJ Format ...."
!
!WRITE(fname,999) it
!999 FORMAT(I6.6)
!
!Check_Option: SELECT CASE (Output_Option)
!CASE (1)
!! Field Phi, Rho (IJK format), Species Vj (IJK format), Phase (IJK format)
!


!! Phase in TecPlot format: IJK Ordered
!filename = 'Phase_'//TRIM(fname)//'.dat'
!PRINT*, 'Writing to file: ', TRIM(filename)
!OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))	
!WRITE(1,*) 'TITLE = "Phase Plot"'
!WRITE(1,*) 'VARIABLES = "x" "y" "z" "Vx" "Vy" "Vz" "isp"'
!DO i = 1, ntot
!	IF (MOD(i, n_stride(1))==0) THEN
!		WRITE(1,80) part(i,1:3), part(i,4:6), INT(part(i,7))
!	END IF
!END DO
!80 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',I2)
!CLOSE(1)
!
!CASE (2)
!! Position (IJK format)
!
!! Trace in TecPlot format: IJK Ordered
!!filename = 'Trace_'//TRIM(fname)//'.dat'
!!PRINT*, 'Writing to file: ', TRIM(filename)
!!OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))	
!!WRITE(1,*) 'TITLE = "Trace Plot"'
!!!WRITE(1,*) 'VARIABLES = "x" "y" "isp"'
!!WRITE(1,*) 'VARIABLES = "x" "y" "z" "Vx" "Vy" "Vz" "isp"'
!!DO i = 1, ntot
!!	IF (MOD(i, n_stride(1))==0) THEN
!!!		WRITE(1,90) part(i,1:2), INT(part(i,7))
!!		WRITE(1,80) part(i,1:3), part(i,4:6), INT(part(i,7))
!!	END IF
!!END DO
!!90 FORMAT (E15.6,' ',E15.6,' ',I2)
!!CLOSE(1)
!
!
!END SELECT Check_Option
!
!PRINT*, "Saving Data in Tecplot IJK Format Done."

WRITE(6,*) 'Output_To_Tecplot_IJK_2D'


    n_node = SIZE(HP,2)
    n_element = SIZE(HT,2)
    WRITE(num_n,"(I10)")  n_node
    WRITE(num_e,"(I10)")  n_element
    n_node_part = SIZE(HP,2)
    radius=0.0
    Rho_temp=0.0
        do j=1, SIZE(HP,2)
            radius(j)=SQRT (HP(1,j)*HP(1,j)+HP(2,j)*HP(2,j))
        end do 
    
    !radius_max=MAXVAL (radius)
    !j_part=1.0
    
    !do i=1,SIZE (HP,2)
    !    Rho_temp(i, 1, 2)=Rho_s(i, 1, 2)
    !    Rho_temp(i, 1, 1)=Rho_s(i, 1, 1)
    !end do
    !

    !do i=1,2
    !    do j=1,SIZE(HP,2)
    !        do i_part=i,SIZE(HP,2)
    !            if(radius(j)==radius(i_part))then 
    !            Rho_temp(j,1,i)=Rho_temp(j,1,i)+Rho_temp(i_part,1,i)
    !            Rho_temp(i_part,1,i)=0
    !            j_part=j_part+1.0
    !            radius(i_part)=radius_max+j_part
    !            end if 
    !        end do
    !    end do
    !end do
           
    
    
        

    PRINT*, "Saving Data in Tecplot IJ Format ...."
    WRITE(fname,999) it
999 FORMAT(I6.6)
    
!    filename = './OUTPUT/rho/radius_IJ_'//TRIM(fname)//'.dat'
!    OPEN(317, ACTION = 'WRITE', FILE = TRIM(filename))
!    WRITE(317,('(A150)')) 'VARIABLES = "radius" "Rho_ele" "Rho_ion"'
!    DO i = 1, SIZE(HP,2)
!          WRITE(317,17) radius(i), Rho_temp(i,1,1), Rho_temp(i,1,2)
!    END DO
!17 format(E15.6,' ',E15.6,' ',E15.6)    
        
    Check_Option: SELECT CASE (Output_Option)
        CASE (1)
        !==========the value of Phi, Rho, Rho_ele, Rho_ion==========
        filename = './OUTPUT/Field/field_IJ_'//TRIM(fname)//'.dat'
        PRINT*, 'Writing to file: ', TRIM(filename)
        OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))
        WRITE(1,*) 'TITLE = "Field Plot"'
        WRITE(1,('(A150)')) 'VARIABLES = "x" "y" "R" "Phi" "Rho" "Rho_ele" "Rho_H" "Rho_C" "efx" "efy" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion" "node_type1"  "node_type2"'
        WRITE(1,*) 'ZONE N='//num_n//', E='//num_e//',F=FEPOINT,ET=QUADRILATERAL'
        DO i = 1, SIZE(HP,2)
          WRITE(1,50) HP(1:2,i), radius(i),Phi(i,1), Rho(i,1), Rho_s(i,1,1), Rho_s(i,1,2),Rho_s(i,1,3), efx(i, 1), efy(i, 1), Ek_s(i,1,1), Ek_s(i,1,2), Ek_tot(i,1,1), Ek_tot(i,1,2), node_type(1, i), node_type(2, i)
        END DO
        DO i = 1, SIZE(HT,2) 
          WRITE(1,44) HT(1, i),HT(2, i),HT(3, i),HT(4, i)
        END DO
        44 FORMAT (I15, I15, I15, I15)
        50 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,'',E15.6,' ',E15.6,' ',E15.6, E15.6,' ',E15.6,' ',E15.6,' ',E15.6, ' ', I15, ' ', I15)
        CLOSE(1)
        !==========the value of Phi, Rho, Rho_ele, Rho_ion==========

 
        NumOfNodeX = size(node_map,2)
        
        !---real time output---- 
        ALLOCATE(Phi_aver_axis(1:NumOfNodeX), Rho_aver_axis(1:NumOfNodeX), Rho_ele_aver_axis(1:NumOfNodeX), Rho_ion_aver_axis(1:NumOfNodeX))
        Allocate(Ek_s_ele_aver_axis(1:NumOfNodeX), Ek_s_ion_aver_axis(1:NumOfNodeX), Ek_tot_ele_aver_axis(1:NumOfNodeX), Ek_tot_ion_aver_axis(1:NumOfNodeX))
         Phi_aver_axis(1:NumOfNodeX) = 0.0
         Rho_aver_axis(1:NumOfNodeX) = 0.0
         Rho_ele_aver_axis(1:NumOfNodeX) = 0.0
         Rho_ion_aver_axis(1:NumOfNodeX) = 0.0
         Ek_s_ele_aver_axis(1:NumOfNodeX) = 0.0
         Ek_s_ion_aver_axis(1:NumOfNodeX) = 0.0
         Ek_tot_ele_aver_axis(1:NumOfNodeX) = 0.0
         Ek_tot_ion_aver_axis(1:NumOfNodeX) = 0.0
         
         !h_partition = VertX(1,2,1)-VertX(1,1,1)
         !WRITE(6,*) 'Output_To_Tecplot: h_partition = ', h_partition
         h_flag = 0.0
         DO j = 1, NumOfNodeX
  
           !h_flag = VertX(1,1,1) + (j-1)*h_partition
            h_flag = node_map(2,j)
             
           count_grid = 0
           DO i = 1, SIZE(HP,2)
    
             IF (ABS(h_flag-HP(1,i))<SmallValue) THEN
               Phi_aver_axis(j) = Phi_aver_axis(j) + phi(i, 1)
               Rho_aver_axis(j) = Rho_aver_axis(j) + Rho(i, 1)
               Rho_ele_aver_axis(j) = Rho_ele_aver_axis(j) + rho_s(i, 1, 1)
               Rho_ion_aver_axis(j) = Rho_ion_aver_axis(j) + rho_s(i, 1, 2)
               Ek_s_ele_aver_axis(j) = Ek_s_ele_aver_axis(j) + Ek_s(i, 1, 1)
               Ek_s_ion_aver_axis(j) = Ek_s_ion_aver_axis(j) + Ek_s(i, 1, 2)
               Ek_tot_ele_aver_axis(j) = Ek_tot_ele_aver_axis(j) + Ek_tot(i, 1, 1)
               Ek_tot_ion_aver_axis(j) = Ek_tot_ion_aver_axis(j) + Ek_tot(i, 1, 2)
               count_grid = count_grid + 1
             END IF
           END DO
  
           IF (count_grid==0) THEN
             WRITE(6,*) 'Output_To_Tecplot: count_grid = ', count_grid
             WRITE(6,*) 'Output_To_Tecplot: h_flag = ', h_flag
             STOP
           END IF
  
           Phi_aver_axis(j) = Phi_aver_axis(j) / count_grid
           Rho_aver_axis(j) = Rho_aver_axis(j) / count_grid
           Rho_ele_aver_axis(j) = Rho_ele_aver_axis(j) / count_grid
           Rho_ion_aver_axis(j) = Rho_ion_aver_axis(j) / count_grid
           Ek_s_ele_aver_axis(j) = Ek_s_ele_aver_axis(j) / count_grid
           Ek_s_ion_aver_axis(j) = Ek_s_ion_aver_axis(j) / count_grid
           Ek_tot_ele_aver_axis(j) = Ek_tot_ele_aver_axis(j) / count_grid
           Ek_tot_ion_aver_axis(j) = Ek_tot_ion_aver_axis(j) / count_grid
         END DO

         filename = './OUTPUT/Field/Average_x_'//TRIM(fname)//'.dat'
         PRINT*, 'Writing to file: ', TRIM(filename)
         OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))
         WRITE(1,*) 'TITLE = "Field Plot"'
         WRITE(1,*) 'VARIABLES = "x" "Phi" "Rho" "Rho_ele" "Rho_ion" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion"'
         WRITE(1,60) NumOfNodeX
         DO j=1,NumOfNodeX
           WRITE(1,70) node_map(2,j)*L_ref, Phi_aver_axis(j)*Phi_ref, Rho_aver_axis(j)*n_ref, Rho_ele_aver_axis(j)*n_ref, Rho_ion_aver_axis(j)*n_ref, &
                    Ek_s_ele_aver_axis(j), Ek_s_ion_aver_axis(j), Ek_tot_ele_aver_axis(j), Ek_tot_ion_aver_axis(j)
         END DO
         !60  FORMAT (' ZONE I = ',I6)
         !70 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ')
         CLOSE(1)
        !==========the average value of Field on structured grid(No SIDG mesh)==========
        DEALLOCATE(Phi_aver_axis, Rho_aver_axis, Rho_ele_aver_axis, Rho_ion_aver_axis, &
                Ek_s_ele_aver_axis, Ek_s_ion_aver_axis, Ek_tot_ele_aver_axis, Ek_tot_ion_aver_axis)
    
        CASE (3)
        NumOfNodeX = size(node_map,2)
        !---avg output 32T---- 
         ALLOCATE(Phi_aver_axis(1:NumOfNodeX), Rho_aver_axis(1:NumOfNodeX), Rho_ele_aver_axis(1:NumOfNodeX), Rho_ion_aver_axis(1:NumOfNodeX))
         Allocate(Ek_s_ele_aver_axis(1:NumOfNodeX), Ek_s_ion_aver_axis(1:NumOfNodeX), Ek_tot_ele_aver_axis(1:NumOfNodeX), Ek_tot_ion_aver_axis(1:NumOfNodeX))
         Phi_aver_axis(1:NumOfNodeX) = 0.0
         Rho_aver_axis(1:NumOfNodeX) = 0.0
         Rho_ele_aver_axis(1:NumOfNodeX) = 0.0
         Rho_ion_aver_axis(1:NumOfNodeX) = 0.0
         Ek_s_ele_aver_axis(1:NumOfNodeX) = 0.0
         Ek_s_ion_aver_axis(1:NumOfNodeX) = 0.0
         Ek_tot_ele_aver_axis(1:NumOfNodeX) = 0.0
         Ek_tot_ion_aver_axis(1:NumOfNodeX) = 0.0
         
         !h_partition = VertX(1,2,1)-VertX(1,1,1)
         !WRITE(6,*) 'Output_To_Tecplot: h_partition = ', h_partition
         h_flag = 0.0
         DO j = 1, NumOfNodeX
  
           !h_flag = VertX(1,1,1) + (j-1)*h_partition
           h_flag = node_map(2,j)
           
           count_grid = 0
           DO i = 1, SIZE(HP,2)
    
             IF (ABS(h_flag-HP(1,i))<SmallValue) THEN
               Phi_aver_axis(j) = Phi_aver_axis(j) + Aver_phi(i, 1)
               Rho_aver_axis(j) = Rho_aver_axis(j) + Aver_Rho(i, 1)
               Rho_ele_aver_axis(j) = Rho_ele_aver_axis(j) + Aver_rho_s(i, 1, 1)
               Rho_ion_aver_axis(j) = Rho_ion_aver_axis(j) + Aver_rho_s(i, 1, 2)
               Ek_s_ele_aver_axis(j) = Ek_s_ele_aver_axis(j) + Aver_Ek_s(i, 1, 1)
               Ek_s_ion_aver_axis(j) = Ek_s_ion_aver_axis(j) + Aver_Ek_s(i, 1, 2)
               Ek_tot_ele_aver_axis(j) = Ek_tot_ele_aver_axis(j) + Aver_Ek_tot(i, 1, 1)
               Ek_tot_ion_aver_axis(j) = Ek_tot_ion_aver_axis(j) + Aver_Ek_tot(i, 1, 2)
               count_grid = count_grid + 1
             END IF
           END DO
  
           IF (count_grid==0) THEN
             WRITE(6,*) 'Output_To_Tecplot: count_grid = ', count_grid
             WRITE(6,*) 'Output_To_Tecplot: h_flag = ', h_flag
             STOP
           END IF
  
           Phi_aver_axis(j) = Phi_aver_axis(j) / count_grid
           Rho_aver_axis(j) = Rho_aver_axis(j) / count_grid
           Rho_ele_aver_axis(j) = Rho_ele_aver_axis(j) / count_grid
           Rho_ion_aver_axis(j) = Rho_ion_aver_axis(j) / count_grid
           Ek_s_ele_aver_axis(j) = Ek_s_ele_aver_axis(j) / count_grid
           Ek_s_ion_aver_axis(j) = Ek_s_ion_aver_axis(j) / count_grid
           Ek_tot_ele_aver_axis(j) = Ek_tot_ele_aver_axis(j) / count_grid
           Ek_tot_ion_aver_axis(j) = Ek_tot_ion_aver_axis(j) / count_grid
         END DO

         filename = './OUTPUT/Average/Average_J_'//TRIM(fname)//'.dat'
         PRINT*, 'Writing to file: ', TRIM(filename)
         OPEN(1, ACTION = 'WRITE', FILE = TRIM(filename))
         WRITE(1,*) 'TITLE = "Field Plot"'
         WRITE(1,*) 'VARIABLES = "x" "Phi" "Rho" "Rho_ele" "Rho_ion" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion"'
         WRITE(1,60) NumOfNodeX
         DO j=1,NumOfNodeX
           WRITE(1,70) node_map(2,j)*L_ref, Phi_aver_axis(j)*Phi_ref, Rho_aver_axis(j)*n_ref, Rho_ele_aver_axis(j)*n_ref, Rho_ion_aver_axis(j)*n_ref, &
                Ek_s_ele_aver_axis(j), Ek_s_ion_aver_axis(j), Ek_tot_ele_aver_axis(j), Ek_tot_ion_aver_axis(j)
         END DO
         60  FORMAT (' ZONE I = ',I6)
         70 FORMAT (E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6,' ',E15.6)
         CLOSE(1)
        !==========the average value of Field on structured grid(No SIDG mesh)==========

         DEALLOCATE(Phi_aver_axis, Rho_aver_axis, Rho_ele_aver_axis, Rho_ion_aver_axis, &
         Ek_s_ele_aver_axis, Ek_s_ion_aver_axis, Ek_tot_ele_aver_axis, Ek_tot_ion_aver_axis)
 
 
END SELECT Check_Option

PRINT*, "Saving Data in Tecplot IJK Format Done."
!----------LY REVISE for SIDG-PIC, 2022-1-12----------
END SUBROUTINE

