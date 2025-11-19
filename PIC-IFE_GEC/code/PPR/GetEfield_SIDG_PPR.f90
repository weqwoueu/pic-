SUBROUTINE GetEField_SIDG_PPR
    
  USE IFE_Data
  USE Field_2D
  USE Domain_2D
  USE Object_Data_2D
  USE Gauss_Data
  USE IFE_INTERFACE, ONLY :Generate_Gauss_local_ppr,Generate_Gauss_local_triangle,&
                          Linear_IFE_Basis_Coeff_2D,&
                          Linear_FE_Basis_Coeff_2D , Linear_FE_Basis_Eval_2D, Inter_Tetra_Partition_2D,&
                          BRINV, Bubble_sort_new
  !=========LY modification for Periodic boundary and Neumann boundary, 2022-7-8=========
  Use IFE_Boundary
  !=========LY modification for Periodic boundary and Neumann boundary, 2022-7-8=========
  
  IMPLICIT NONE
  
  ! ============wsy add=================
  INTEGER,PARAMETER :: In_index = -2, Out_index = -1
  !=====================================

  INTEGER nnx,nny
  !TYPE(ObjectType), DIMENSION(:), POINTER :: objects
  INTEGER  i, j, m, n, num, k, st
  INTEGER  n_node_sidg, n_elem_sidg
  INTEGER  row, col
  
  
  !find adjacent element for all nodes
  INTEGER num_node_in_elem
  
  
  !build point patch information
  INTEGER, DIMENSION(:,:), POINTER ::  ab_point_patch
  !INTEGER, DIMENSION(:,:), POINTER :: point_patch
  
  INTEGER  n_node_in_one_patch
  
  !recover gradient
  INTEGER   NI, ith, jth, aa, x, y, L
  INTEGER ::IS(6), JS(6)
  REAL(8)   h, h_temp, xx, yy, temp1, temp2, error_inf, r1
  REAL(8), DIMENSION(:,:), ALLOCATABLE::	phi_vec
  
  !REAL(8), DIMENSION(:,:), POINTER :: coef_ppr
  !REAL(8), DIMENSION(:,:), POINTER :: re_grad
  REAL(8), DIMENSION(:,:), POINTER :: pnts
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: phi_1!,efx_1,efy_1
  REAL(8), DIMENSION(:,:), POINTER :: AT,A,A_MUL
  REAL(8), DIMENSION(:,:), POINTER :: BX,BY
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: temp_coef
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: temp_a,temp_b
  REAL(8)                          ::efx_true,efy_true,x_cor,y_cor
  REAL(8)                          :: H1_semi_norm_ppr, ab_errorx, ab_errory
  REAL(8)                          :: h1, h2
  REAL(8), DIMENSION(:,:),POINTER  :: efx_temp_x, efy_temp_y
  
  
  TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
  REAL(8)                               :: x_loc,y_loc,   x_loc1,y_loc1
  INTEGER, DIMENSION(:), POINTER        :: ab_ppr_patch_index_temp, ab_ppr_patch_index
  REAL(8), DIMENSION(:,:), POINTER      :: ab_ppr_patch
  REAL(8), DIMENSION(:), POINTER        :: ab_ppr_point_x, ab_ppr_point_y, ab_ppr_point
  
  INTEGER                                         n_sub_pt, sub_region_ind, el_type
  INTEGER, DIMENSION(:,:), POINTER		  ::	    t_sub_pt
  REAL(8), DIMENSION(:,:), POINTER		  ::	    p_sub_pt
  INTEGER                                         dind(2),pointer_reference_to_local(4)
  INTEGER, DIMENSION(:), POINTER		    ::      tri_sign
  REAL(8)                                         eint(9),u_val_x(9), u_val_y(9),&
                                                u_val_ppr_x(9),u_val_ppr_y(9),err_val(9),&
                                                x_cord(9),y_cord(9)
  REAL(8), DIMENSION(:,:),  POINTER     ::      Intrs_pts
  REAL(8), DIMENSION(:), ALLOCATABLE		::	gwght, gwght_IFE
  REAL(8), DIMENSION(:,:), ALLOCATABLE	::	gnodes, gnodes_IFE
  REAL(8)                                     error_sum
  INTEGER                                         pts_index(2,4)
  REAL(8)                                     vert(2,4), vert_st(2,3)
  REAL(8)                                         coef, r0
  REAL(8)                                         dis, dis_x,dis_y
  INTEGER                                     it, res, ind, ie, o_ind
  REAL(8)                                         node_ave
  REAL(8)                                     beta_minus,beta_plus
  INTEGER                                      part_num, part_index
  INTEGER, DIMENSION(:), ALLOCATABLE		    :: repeat_times_ppr
  !INTEGER, DIMENSION(:,:), ALLOCATABLE  	  :: 
  REAL(8)                                   :: part_x, part_y
  INTEGER                                   :: part_ex_1, part_ex_2, part_ex_3, part_ex_4
  INTEGER, DIMENSION(:), ALLOCATABLE		    :: ab_sample_node
  REAL(8), DIMENSION(:,:), ALLOCATABLE      :: ab_sample_cor
  REAL(8), DIMENSION(:), ALLOCATABLE        :: ab_point_x, ab_point_y
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: ab_point
  REAL(8), DIMENSION(:), ALLOCATABLE        :: ab_efx
  REAL(8), DIMENSION(:), ALLOCATABLE        :: ab_efy
  REAL(8)                                   :: x_cordi, y_cordi
  INTEGER                                   :: num_intrs_elem
  ! zyz new add for ?
    INTEGER, DIMENSION(:), POINTER            :: temp
  INTEGER                                   :: length,node_ab
  REAL(8)                                   :: x_re_cor, y_re_cor
  REAL(8)                                   ::  error_ini ,error_temp,phi_temp
  REAL(8)     index,r
  !REAL(8), DIMENSION(:), ALLOCATABLE        ::
  !INTEGER, DIMENSION(:), ALLOCATABLE      :: re_times
  
  ! wsy add for linear interpolation
  Integer                                   :: Node1, Node2, Node3, Node4
  Real(8)                                      :: dx(9), xcellmdx(9), dy(9), ycellmdy(9), P1(9), P2(9), P3(9), P4(9)
  
  !LY modification for Periodic boundary and Neumann boundary, 2022-7-8
  Real(8) :: left_boundary
  Real(8) :: right_boundary
  Real(8) :: bottom_boundary
  Real(8) :: top_boundary
  
  
  !r0 = 0.500253607
  !r0 = 3.1415926
  !beta_minus = Global_Beta(2)  ! 内
  !beta_plus  = Global_Beta(1)  ! 外

  
  n_elem_sidg  = size(HT,2)
  n_node_sidg  = size(HP,2)
  
  ALLOCATE(phi_1(n_node_sidg,1))!把phi(i,j)转化为向量phi_1
  !HP(1:2,i), Phi(i,1), i所以在这里phi和hp里面的点是一一对应的关系。
  DO i=1,n_node_sidg
      phi_1(i,1)=phi(i,1)
  ENDDO
  
  num_node_in_elem = 4
  num_intrs_elem =MAXVAL(element_index)
  

  WRITE(6,*) 
  WRITE(6,*) 'GetEfield from PPR method'

!===========================开始重构===================================================
!=====================================================================================
  WRITE(*,*) 'ready to recover'
  ! IF (.NOT.ALLOCATED(efx_1)) ALLOCATE(efx_1(n_node_sidg,1), efy_1(n_node_sidg,1))
  IF (.NOT.ALLOCATED(coef_ppr)) ALLOCATE(coef_ppr(6,n_node_sidg))
  IF (.NOT.ALLOCATED(temp_coef)) ALLOCATE(temp_coef(6,1))
  IF (.NOT.ALLOCATED(temp_a)) ALLOCATE(temp_a(1,1),temp_b(1,1))
  !ALLOCATE(re_grad(2,n_node))
  efx_1=0.D0; efy_1=0.D0
  coef_ppr=0.D0
  temp_coef=0.D0
  !re_grad=0.D0
  !BX=0.0 ; BY=0.0
  
  DO i=1,n_node_sidg
    n_node_in_one_patch=point_patch(i,1)
    IF(n_node_in_one_patch/=0)THEN  
      IF(PPRPointIndex(i) == 1 .Or. .True.)THEN
        IF(re_times(i)==0) THEN
          ALLOCATE(BX(1,n_node_in_one_patch),BY(1,n_node_in_one_patch))
          ALLOCATE(pnts(2,n_node_in_one_patch))
          BX=0.0 ; BY=0.0
          pnts=0.0
          h=0
          h_temp=0
          ALLOCATE(phi_vec(n_node_in_one_patch,1))
          DO j=1,n_node_in_one_patch
            k=point_patch(i,j+1)
            pnts(:,j)=HP(:,k)
            phi_vec(j,1)=phi_1(k,1) !bh,方程右侧的解向量,只有六个元素
        
            h_temp=MAX(ABS(pnts(1,j)-HP(1,i)),ABS(pnts(2,j)-HP(2,i)))
            IF(h_temp>h) THEN 
              h=h_temp
            ENDIF
          ENDDO
          !from local coordinate to reference coordinate
          pnts(1,:)=(pnts(1,:)-HP(1,i))/h  !ζ=(x-xi)/h
          pnts(2,:)=(pnts(2,:)-HP(2,i))/h  !η=(y-yi)/h
      
          ALLOCATE(AT(6,n_node_in_one_patch))
          ALLOCATE(A(n_node_in_one_patch,6))
          ALLOCATE(A_MUL(6,6))
          AT(1,:)=1.0                         !1
          AT(2,:)=pnts(1,:)                   !ζ
          AT(3,:)=pnts(2,:)                   !η
          DO num=1,n_node_in_one_patch
            AT(4,num)=pnts(1,num)*pnts(1,num) !ζ**2
            AT(5,num)=pnts(1,num)*pnts(2,num) !ζ*η
            AT(6,num)=pnts(2,num)*pnts(2,num) !η**2
          ENDDO
      
          A=TRANSPOSE(AT)
          A_MUL=MATMUL(AT,A)
          CALL BRINV(A_MUL,6,L,IS,JS)!(AT*A)-1 6*6
          AT=MATMUL(A_MUL,AT)!AT=((AT*A)-1)*AT   6*n_node_in_one_patch
          DO num=1,n_node_in_one_patch
            BX(1,num)=AT(2,num)!每执行一个i的循环得到一个1*6的行向量，总共n_node次循环，BXi只有n_node_in_one_patch个数据不为零，其他为零
            BY(1,num)=AT(3,num)
          ENDDO
          temp_a = (MATMUL(BX,phi_vec))/h
          temp_b = (MATMUL(BY,phi_vec))/h
          efx_1(i,1) = temp_a(1,1)
          efy_1(i,1) = temp_b(1,1)
          !DO num=1,n_node_in_one_patch
          !  AT2(num)=AT(2,num)!/h
          !  AT3(num)=AT(3,num)!/h
          !  AT4(num)=AT(4,num)!/(h*h)
          !  AT5(num)=AT(5,num)!/(h*h)
          !  AT6(num)=AT(6,num)!/(h*h)
          !ENDDO
          !
          !a1_hat=MATMUL(AT2,phi_vec)
          !a2_hat=MATMUL(AT3,phi_vec)
          !a3_hat=MATMUL(AT4,phi_vec)
          !a4_hat=MATMUL(AT5,phi_vec)
          !a5_hat=MATMUL(AT6,phi_vec)
          !
          !a1=a1_hat/h
          !a2=a2_hat/h
          !a3=a3_hat/(h*h)
          !a4=a4_hat/(h*h)
          !a5=a5_hat/(h*h)

          !efx=a1+2*a3*x+a4*y
          !efy=a2+2*a5*y+a4*x!这里的xy是相对坐标还是绝对坐标？
          !efx_1(i,1)=DOT_PRODUCT(BX,phi_vec)!重构的梯度x
          !efy_1(i,1)=DOT_PRODUCT(BY,phi_vec)!重构的梯度y

          temp_coef=MATMUL(AT,phi_vec) !高阶近似的系数向量a帽,6个元素的列向量
          coef_ppr(1:6,i)=temp_coef(:,1)
          coef_ppr(2:3,i)=temp_coef(2:3,1)/h
          coef_ppr(4:6,i)=temp_coef(4:6,1)/(h*h)!除以h之后得到的是a，不是a帽了coef的第i列是第i个节点的a向量，从a0到a5
          DEALLOCATE(BX,BY)
          DEALLOCATE(A_MUL,A,AT,pnts,phi_vec)
        ELSE
          efx_1(i,1) = efx_1(re_times(i),1)
          efy_1(i,1) = efy_1(re_times(i),1)
          coef_ppr(1:6,i)=coef_ppr(1:6,re_times(i))
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  WRITE(*,*) 'finish recovering gradient'
  
  

!===========================完成PPR重构===================================================
  
  
!===================特殊点的处理:公式法=============================================  
  !ALLOCATE(ab_efx(n_node_sidg), ab_efy(n_node_sidg))
  !ab_efx=0.D0
  !ab_efy=0.D0
  !
  !DO i =1,n_node_sidg
  !  IF(point_patch(i,1)== 0)THEN
  !    IF(re_times(1,i)==0) THEN
  !      DO j=1,16
  !        IF(ab_point(i,j)/=0) THEN
  !          x_cordi = HP(1,i)-HP(1,ab_point(i,j))
  !          y_cordi = HP(2,i)-HP(2,ab_point(i,j))
  !          ab_efx(i) = coef_ppr(2,ab_point(i,j)) + 2*coef_ppr(4,ab_point(i,j))*x_cordi + coef_ppr(5,ab_point(i,j))*y_cordi
  !          ab_efy(i) = coef_ppr(3,ab_point(i,j)) + 2*coef_ppr(6,ab_point(i,j))*y_cordi + coef_ppr(5,ab_point(i,j))*x_cordi
  !          k=k+1
  !        ENDIF
  !      ENDDO
  !      ab_efx(i) = ab_efx(i)/k
  !      ab_efy(i) = ab_efy(i)/k
  !      k=0
  !      efx_1(i,1) = ab_efx(i)
  !      efy_1(i,1) = ab_efy(i)
  !    ELSE
  !      efx_1(i,1) = efx_1(re_times(1,i),1)
  !      efy_1(i,1) = efy_1(re_times(1,i),1)
  !      coef_ppr(1:6,i)=coef_ppr(1:6,re_times(1,i))
  !    ENDIF
  !  ENDIF
  !ENDDO
  !WRITE(*,*) 'abnormal nodes treatment, done'
   
!=========================特殊点用其他点来平均=============================================
!========================================================================================  

  !DO i =1,n_node_sidg
  !  IF(point_patch(i,1)== 0)THEN
  !    x_loc = HP(1,i)
  !    y_loc = HP(2,i)
  !    DO j=1,ab_point_patch(i,1)
  !      
  !      n = ab_point_patch(i,j+1)
  !      x_re_cor=x_loc - HP(1,ab_point_patch(i,j+1))
  !      y_re_cor=y_loc - HP(2,ab_point_patch(i,j+1))
  !      efx_1(i,1) = efx_1(i,1)+ coef_PPR(2,n) + 2.* coef_PPR(4,n)*x_re_cor + coef_PPR(5,k)*y_re_cor
  !      efy_1(i,1) = efy_1(i,1)+ coef_PPR(3,n) + 2.* coef_PPR(6,n)*y_re_cor + coef_PPR(5,k)*x_re_cor
  !    ENDDO
  !    efx_1(i,1) = efx_1(i,1)/ab_point_patch(i,1)
  !    efy_1(i,1) = efy_1(i,1)/ab_point_patch(i,1)
  !  ENDIF
  !ENDDO
  
  !=========Electric Field final value, it will moving particle=========
  DO i=1,n_node_sidg
    efx(i,1) = -efx_1(i,1)
    efy(i,1) = -efy_1(i,1)
  End Do
  
  !======LY modification for Periodic boundary and Neumann boundary condition, 2022-7-8======
  !left_boundary = MINVAL(HP(1,1:Size(HP,2)))
  !right_boundary = MAXVAL(HP(1,1:Size(HP,2)))
  !bottom_boundary = MINVAL(HP(2,1:Size(HP,2)))
  !top_boundary = MAXVAL(HP(2,1:Size(HP,2)))
  
  !======Periodic boundary condition======
  !If (bc_index(1)==-1 .OR. bc_index(2)==-1 .OR. bc_index(3)==-1 .OR. bc_index(4)==-1) Then
  !  If (bc_index(1) == -1) Then !Left and Right periodic boundary
  !    Do i = 1, Size(PairNodes,2)
  !      efx(PairNodes(1,i),1) = (efx(PairNodes(1,i),1) + efx(PairNodes(2,i),1)) / 2.0
  !      efx(PairNodes(2,i),1) = efx(PairNodes(1,i),1)
  !    End Do
  !  End If
  !  
  !  If (bc_index(3) == -1) Then !Bottom and Top periodic boundary
  !    Do i = 1, Size(PairNodes,2)
  !      efy(PairNodes(1,i),1) = (efy(PairNodes(1,i),1) + efy(PairNodes(2,i),1)) / 2.0
  !      efy(PairNodes(2,i),1) = efy(PairNodes(1,i),1)
  !    End Do
  !  End If
  !End If
  !======Periodic boundary condition======
  
  !======Zero-Neumann boundary condition======
  !If (bc_index(1)==0 .OR. bc_index(2)==0 .OR. bc_index(3)==0 .OR. bc_index(4)==0) Then
  !  Do i = 1, Size(HP,2)
  !    If (bc_index(1) == 0) Then      !Left field boundary
  !      If (ABS(HP(1,i)-left_boundary) < SmallValue) Then
  !        efx(i,1) = 0._8
  !      End If
  !    Elseif (bc_index(2) == 0) Then  !Right field boundary
  !      If (ABS(HP(1,i)-right_boundary) < SmallValue) Then
  !        efx(i,1) = 0._8
  !      End If
  !    Elseif (bc_index(3) == 0) Then  !Bottom field boundary
  !      If (ABS(HP(2,i)-bottom_boundary) < SmallValue) Then
  !        efy(i,1) = 0._8
  !      End If
  !    Elseif (bc_index(4) == 0) Then  !Top field boundary
  !      If (ABS(HP(2,i)-top_boundary) < SmallValue) Then
  !        efy(i,1) = 0._8
  !      End If
  !    End If
  !  End Do
  !End If
  !======Zero-Neumann boundary condition======
  !======LY modification for Periodic boundary and Neumann boundary condition, 2022-7-8======
  !=========Electric Field final value, it will moving particle========= 
    
    
    !ab_errorx = 0
    !ab_errory = 0
    !Do i = 1 , n_node_sidg
    !    xx = HP(1,i)
    !    yy = HP(2,i)
    !    if ((xx)**2 + (yy)**2 > r0**2) then
    !        !efx(i,1) = -3*(xx-5)*((xx - 5)**2 + (yy - 5)**2)**(0.5)/10
    !        !efy(i,1) = -3*(yy-5)*((xx - 5)**2 + (yy - 5)**2)**(0.5)/10
    !        
    !        temp1 = abs( 3*(xx)*((xx)**2 + (yy)**2)**(0.5)/10 +  efx(i,1))
    !        temp2 = abs( 3*(yy)*((xx )**2 + (yy )**2)**(0.5)/10 +  efy(i,1))
    !        if (temp2 > ab_errory) then
    !            aa = i
    !        End if
    !        ab_errorx = max(ab_errorx, temp1)
    !        ab_errory = max(ab_errory, temp2)
    !        
    !        
    !    else 
    !        temp1 = abs( 3*(xx)*((xx)**2 + (yy)**2)**(0.5) +  efx(i,1))
    !        temp2 = abs( 3*(yy)*((xx)**2 + (yy)**2)**(0.5) +  efy(i,1))
    !        if (temp2 > ab_errory) then
    !            aa = i
    !        End if
    !        ab_errorx = max(ab_errorx, temp1)
    !        ab_errory = max(ab_errory, temp2)
    !    end if
    !End DO
    !write(*,*) 'ab_errorx = ',ab_errorx
    !write(*,*) 'ab_errory = ',ab_errory
!error_ini=0
!!ALLOCATE(error_temp(1,1))
!!error_temp = 0.0
!index=1.1
!DO i = 1, SIZE(HP,2)
!  CALL Function_True(0, index, HP(1, i), HP(2, i), 0, 0, r)
!  !error_temp(1,1) = ABS(Phi(i, 1) - temp)
!  phi_temp =Phi(i, 1)
!  error_temp = ABS(phi_temp - r)
!  
!  if(error_temp>error_ini) then
!    error_ini = error_temp
!  endif
!  
!END DO
!write(*,*) error_ini
  
  

!=====================================完成平均================================================= 

  !!cg的时候可以这样，但是dg的时候不能这样了，行列不再是这样排布
  !DO i=1,n_node_cg
  !  row=MOD(i,nny)
  !  IF(row==0) THEN
  !    row=i/nny
  !    col=nny
  !  ELSE
  !    col=MOD(i,nny)
  !    row=(i-col)/nny+1
  !  ENDIF
  !  efx(row,col)=efx_1(i,1)
  !  efy(row,col)=efy_1(i,1)
  !ENDDO
  
  
  !=========Only Using in calculate IFE-Efield error norm, not using PIC========= 
  !=========Only Using in calculate IFE-Efield error norm, not using PIC========= 
  !=========Only Using in calculate IFE-Efield error norm, not using PIC========= 
  !!===================================jyj误差分析================================================
  !!=============================================================================================
  ! ALLOCATE(Intrs_pts(2,2))
  ! ALLOCATE(tri_sign(2))
  ! ALLOCATE(gwght(9), gnodes(2,9)) 
  ! err_val=0.D0
  ! error_sum = 0.D0
  ! pts_index = 0               ! 4个节点所在区域标记
  ! 
  ! DO i=1,n_elem_sidg
  !   vert = HP(:,HT(1:4,i)) 
  !   IF(element_index(i)<0)THEN      ! 非界面网格的四个顶点有patch
  !     hx = ABS(vert(:,3)-vert(:,1))!zyz ================================================================================sidg
  !       CALL Generate_Gauss_local_ppr(Gauss_coefficient_reference_Nine,Gauss_point_reference_Nine, &
  !                                   vert(1:2,1), hx, Gauss_coefficient_local_Nine, Gauss_point_local_Nine) !这里hx需要改动
  !       gwght = Gauss_coefficient_local_Nine
  !       gnodes= TRANSPOSE(Gauss_point_local_Nine)  ! 2*9
  !      
  !       IF(element_index(i)==-1)THEN        ! out
  !           coef =  beta_plus  
  !         
  !       ELSEIF(element_index(i)==-2)THEN    ! in
  !           coef =  beta_minus   
  !       ENDIF
  !     
  !       !u_val_x = 3.0D0*(gnodes(1,:))*( (gnodes(1,:))*(gnodes(1,:))+(gnodes(2,:))*(gnodes(2,:)) )**(0.5)/ coef
  !       !u_val_y = 3.0D0*(gnodes(2,:))*( (gnodes(1,:))*(gnodes(1,:))+(gnodes(2,:))*(gnodes(2,:)) )**(0.5)/ coef
  !       
  !       u_val_x = 2.0D0*gnodes(1,:) * EXP( gnodes(1,:)*gnodes(1,:)+gnodes(2,:)*gnodes(2,:) ) / coef
  !       u_val_y = 2.0D0*gnodes(2,:) * EXP( gnodes(1,:)*gnodes(1,:)+gnodes(2,:)*gnodes(2,:) ) / coef
  !       !u_val_x = 2.0D0*gnodes(1,:) / coef
  !       !u_val_y = 2.0D0*gnodes(2,:) / coef!真值
         !u_val_ppr_x = 0.D0      ! 每次遍历非界面单元都要清零
         !u_val_ppr_y = 0.D0
         !! ------  ppr part ------ 
         !DO j= 1,4
         !    k=HT(j,i)
         !    x_cord(:) = gnodes(1,:) - HP(1,k)!相对坐标
         !    y_cord(:) = gnodes(2,:) - HP(2,k)
         !  
         !    u_val_ppr_x = u_val_ppr_x + coef_PPR(2,k) + 2.* coef_PPR(4,k)*x_cord(:) + coef_PPR(5,k)*y_cord(:)
         !    u_val_ppr_y = u_val_ppr_y + coef_PPR(3,k) + 2.* coef_PPR(6,k)*y_cord(:) + coef_PPR(5,k)*x_cord(:)
         !  
         !ENDDO
         !u_val_ppr_x = u_val_ppr_x / 4.      ! 取平均
         !u_val_ppr_y = u_val_ppr_y / 4.
         !! ---------------
         !
         !!--------------linear interpolation------------------
         !!Node1 = HT(1,i)
         !!Node2 = HT(2,i)
         !!Node3 = HT(3,i)
         !!Node4 = HT(4,i)
         !!
         !!dx = (gnodes(1,:) - HP(1, Node1)) / hx(1)
         !!xcellmdx = ( HP(1, Node2) - gnodes(1,:) ) / hx(1)
         !!dy = (gnodes(2,:) - HP(2, Node1)) / hx(2)
         !!ycellmdy = (HP(2, Node4) - gnodes(2 ,:)) / hx(2)
         !!
         !!P1 = xcellmdx *ycellmdy
         !!P2 = dx       *ycellmdy
         !!P3 = xcellmdx *dy
         !!P4 = dx       *dy
         !!
         !!u_val_ppr_x = -(efx(Node1, 1)*P1 + efx(Node2, 1)*P2 + efx(Node3, 1)*P3 + efx(Node4, 1)*P4)
         !!u_val_ppr_y = -(efy(Node1, 1)*P1 + efy(Node2, 1)*P2 + efy(Node3, 1)*P3 + efy(Node4, 1)*P4)
         !!----------------------------------------------
         !
         !
         !! ------------- test -----------------
         !!!Node1 = HT(1,i)
         !!!Node2 = HT(2,i)
         !!!Node3 = HT(3,i)
         !!!Node4 = HT(4,i)
     !    !!
     !    !!h_partition(1) = vert(1,2) - vert(1,1)
     !    !!h_partition(2) = vert(2,4) - vert(2,1)
     !    !!
     !    !!u_val_ppr_x = 0
     !    !!u_val_ppr_y = 0
     !    !!Do n = 1, 9
     !    !!    DO m = 1, 4
     !    !!      CALL Retangular_local_basis(gnodes(1 ,n), gnodes(2 ,n), vert(1:2,1), h_partition, 1, m, 0, 0, r1)
     !    !!      u_val_ppr_x(n) = u_val_ppr_x(n) - efx(HT(m, i), 1) * r1
     !    !!      u_val_ppr_y(n) = u_val_ppr_y(n) - efy(HT(m, i), 1) * r1
     !    !!    END DO
     !    !!End Do
     !    ! ------------------------------------
     !    
     !    err_val = (DABS(u_val_ppr_x - u_val_x))**2. + (DABS(u_val_ppr_y - u_val_y))**2.
     !    error_sum = error_sum + SUM(gwght*err_val)
     !ELSEIF(element_index(i)>0)THEN  ! 界面网格 的四个顶点中 可能有没有patch的顶点
     !
     !    ie   = element_index(i)
     !    pointer_reference_to_local(1:4) = information_1(11:14,ie) 
     !    tri_sign(1)=information_1(15,ie)
	 	  !   tri_sign(2)=information_1(16,ie) 
        !
        ! Intrs_pts(1,1) =information_2(3,ie)     ! Dx
        ! Intrs_pts(1,2) =information_2(4,ie)     ! Dy
        ! Intrs_pts(2,1) =information_2(5,ie)     ! Ex
	 	     !Intrs_pts(2,2) =information_2(6,ie)     ! Ey
        ! el_type = information_1(6,ie)
        !
        ! !print *,el_type
        ! CALL Inter_Tetra_Partition_2D(tri_sign, pointer_reference_to_local, el_type, vert, Intrs_pts, t_sub_pt, p_sub_pt)
        ! n_sub_pt = SIZE(t_sub_pt,2)             ! 4
        !
        ! DO j =1,4
        !    k = HT(j,i)
        !    x_loc = HP(1,k)
        !    y_loc = HP(2,k)    ! 第一行为节点所在区域
        !    pts_index(2,j) = k      ! 第二行为界面网格每个顶点的全局索引
        !    IF (node_index(k)==In_index) THEN
        !        pts_index(1,j) = -2   ! in
        !    ELSE    
        !        pts_index(1,j) = -1   ! out
        !    ENDIF
        ! ENDDO     
        ! !print *,pts_index(1,:)
        ! !print *,pts_index(2,:)
         !!print *,n_sub_pt
         !  
         !DO st = 1,n_sub_pt
         !   vert_st = p_sub_pt(:,t_sub_pt(1:3,st))
         !   CALL Generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle_Nine,Gauss_point_reference_triangle_Nine, &
         !                                     vert_st,Gauss_coefficient_local_triangle_Nine,Gauss_point_local_triangle_Nine)
         !   !print*,vert_st(1,:)
         !   !print*,vert_st(2,:)
         !   !=========================================zyz====================================
         !   IF(el_type == 1) THEN
         !     m=HT(pointer_reference_to_local(st),i)
         !     
         !   ELSEIF(el_type == 2) THEN
         !     IF (st == 1 .OR. st==2) THEN 
         !       n = t_sub_pt(1,st)
         !     ELSEIF(st == 3)THEN
         !       n = t_sub_pt(2,st)
         !     ELSE
         !       n = t_sub_pt(3,st)
         !     ENDIF
            !  m = HT(pointer_reference_to_local(n),i)
            !ENDIF
            !
            !IF(node_index(m) == In_index)THEN
            !      coef =  beta_minus ! in                  
            !ELSE
            !      coef =  beta_plus  ! out
            !ENDIF
            !!=========================================zyz====================================
            !gnodes = TRANSPOSE(Gauss_point_local_triangle_Nine) 
            !gwght  = Gauss_coefficient_local_triangle_Nine
            !
            !!dis_x = ( p_sub_pt(1,t_sub_pt(1, st)) + p_sub_pt(1,t_sub_pt(2, st)) + p_sub_pt(1,t_sub_pt(3, st)) ) /3.
            !!dis_y = ( p_sub_pt(2,t_sub_pt(1, st)) + p_sub_pt(2,t_sub_pt(2, st)) + p_sub_pt(2,t_sub_pt(3, st)) ) /3.
            !!dis   = dsqrt( (dis_x)**2. + (dis_y)**2.)
            !!IF(dis <= r0)THEN
            !!    coef =  beta_minus ! in                  
            !!ELSE
            !!    coef =  beta_plus  ! out
            !!ENDIF
            !
            !!u_val_x = 3.0D0*(gnodes(1,:))*( (gnodes(1,:))*(gnodes(1,:))+(gnodes(2,:))*(gnodes(2,:)) )**(0.5)/ coef
            !!u_val_y = 3.0D0*(gnodes(2,:))*( (gnodes(1,:))*(gnodes(1,:))+(gnodes(2,:))*(gnodes(2,:)) )**(0.5)/ coef
            !u_val_x = 2*gnodes(1,:) * EXP( gnodes(1,:)*gnodes(1,:)+gnodes(2,:)*gnodes(2,:) ) / coef
            !u_val_y = 2*gnodes(2,:) * EXP( gnodes(1,:)*gnodes(1,:)+gnodes(2,:)*gnodes(2,:) ) / coef
            !!u_val_x = 2*gnodes(1,:)  / coef
            !!u_val_y = 2*gnodes(2,:)  / coef
            !
            !node_ave = 0.D0 
            !u_val_ppr_x = 0.D0
            !u_val_ppr_y = 0.D0
            !DO j= 1,4
            !    k = pts_index(2,j)
            !    IF( (pts_index(1,j) == -1 .AND.  coef == beta_plus) .OR. & 
            !        (pts_index(1,j) == -2 .AND.  coef == beta_minus )   )THEN
            !      IF(point_patch(k,1) /= 0) THEN
            !        
            !        x_cord(:) = gnodes(1,:) - HP(1,k)
            !        y_cord(:) = gnodes(2,:) - HP(2,k)
            !        u_val_ppr_x = u_val_ppr_x + coef_PPR(2,k) + 2.* coef_PPR(4,k)*x_cord(:) + coef_PPR(5,k)*y_cord(:)
            !        u_val_ppr_y = u_val_ppr_y + coef_PPR(3,k) + 2.* coef_PPR(6,k)*y_cord(:) + coef_PPR(5,k)*x_cord(:)
            !        node_ave = node_ave + 1
            !      ELSE
            !        IF(st == 1) THEN
            !          !print *,k
            !          node_ab = 0
            !          DO m = 1,ab_point_patch(k,1)
                        !n = ab_point_patch(k,m+1)
    !                    x_cord(:) = gnodes(1,:) - HP(1,n)
    !                    y_cord(:) = gnodes(2,:) - HP(2,n)
    !                    u_val_ppr_x = u_val_ppr_x + coef_PPR(2,n) + 2.* coef_PPR(4,n)*x_cord(:) + coef_PPR(5,n)*y_cord(:)
    !                    u_val_ppr_y = u_val_ppr_y + coef_PPR(3,n) + 2.* coef_PPR(6,n)*y_cord(:) + coef_PPR(5,n)*x_cord(:)
    !                    node_ab = node_ab + 1
    !                  ENDDO
    !                  u_val_ppr_x = u_val_ppr_x / node_ab
    !                  u_val_ppr_y = u_val_ppr_y / node_ab
    !                  !node_ave = node_ab
    !                  node_ave = node_ave + 1
    !                ENDIF
    !              ENDIF   
    !            ENDIF      
    !        ENDDO           
    !        u_val_ppr_x = u_val_ppr_x / node_ave
    !        u_val_ppr_y = u_val_ppr_y / node_ave
    !        err_val = (DABS(u_val_ppr_x-u_val_x))**2. + (DABS(u_val_ppr_y - u_val_y))**2.
    !        error_sum = error_sum + SUM(gwght*err_val)      
    !      ENDDO
    !  ENDIF 
    !ENDDO
    !WRITE(*,*)'H1_error_ppr_new =',DSQRT(error_sum)
    !WRITE(*,*) '---'
  !=========Only Using in calculate IFE-Efield error norm, not using PIC========= 
  !=========Only Using in calculate IFE-Efield error norm, not using PIC========= 
  !=========Only Using in calculate IFE-Efield error norm, not using PIC=========
  
  
  
    !stop
  !DO i=1,n_node_cg
  !  ith= MOD(i,nny)
  !  IF(ith==0)THEN
  !      ith = i/nny
  !      jth = nny
  !  ELSE
  !      jth = MOD(i,nny)
  !      ith = (i-jth)/nny +1 
  !  ENDIF
  !  !----------------------------------------------------------
  !  IF(ith==1)THEN          ! efx_2  x方向 向右为正
  !      efx_2(i,1)= (-3.*phi(1,jth)+4.*phi(2,jth)-phi(3,jth))/(2*hx(1))
  !  ELSEIF(ith ==nnx)THEN
  !      efx_2(i,1)= (3.*phi(nnx,jth)-4.*phi(nnx-1,jth)+phi(nnx-2,jth))/(2*hx(1))
  !  ELSE   
  !      efx_2(i,1)= (phi(ith+1,jth)-phi(ith-1,jth)) /(2*hx(1))
  !  ENDIF
  !  efx_temp_x(ith,jth) = efx_2(i,1)
  !  
  !  
  !  !----------------------------------------------------------
  !  IF(jth==1)THEN          ! efy_2  y方向 向上为正
  !      efy_2(i,1)= (-3.*phi(ith,1)+4.*phi(ith,2) -phi(ith,3))/(2*hx(2))
  !  ELSEIF(jth ==nnx)THEN
  !      efy_2(i,1)= (3.*phi(ith,nny)-4.*phi(ith,nny-1)+phi(ith,nny-2))/(2*hx(2))
  !  ELSE   
  !      efy_2(i,1)= (phi(ith,jth+1)-phi(ith,jth-1)) /(2*hx(2))
  !  ENDIF
  !  efy_temp_y(ith,jth) = efy_2(i,1)    
  !ENDDO
  !write(*,*)'interpolation'
  
  
  
END SUBROUTINE
    
