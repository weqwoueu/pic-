SUBROUTINE Get_PointPatch(nnx, nny, xmin, xmax, ymin, ymax)
    
  USE IFE_Data
  USE Field_2D
  USE Domain_2D
  USE Object_Data_2D
  USE Gauss_Data
  USE IFE_INTERFACE, ONLY :Generate_Gauss_local_ppr,Generate_Gauss_local_triangle,&
                        Linear_IFE_Basis_Coeff_2D,&
                        Linear_FE_Basis_Coeff_2D , Linear_FE_Basis_Eval_2D, Inter_Tetra_Partition_2D,&
                        BRINV, Bubble_sort_new
  IMPLICIT NONE
  
  ! =========wsy add =========
  INTEGER, PARAMETER :: In_index = -2, Out_index = -1
  !=============================

  INTEGER nnx,nny
  !TYPE(ObjectType), DIMENSION(:), POINTER :: objects
  INTEGER  i, j, m, n, num, k, st, xx, yy
  INTEGER  n_node_sidg, n_elem_sidg, n_node_cg, n_elem_cg, elem_index, n_n
  INTEGER  row, col
  
  !find index for boundary points
  INTEGER, DIMENSION(:), POINTER :: boundary_node
  INTEGER, DIMENSION(:), POINTER :: interface_node, temp_node, all_bd_node!边界相关
  INTEGER num_of_bdnode, num_intrs_node
  INTEGER num_intrs_elem, num_all_bd_node
  
  !find adjacent element for all nodes
  INTEGER, DIMENSION(:,:), POINTER :: point_patch_elem, point_patch_elem_dg, point_patch_dg, point_patch_element
  INTEGER, DIMENSION(:,:), POINTER :: point_patch_element_debug! wsy add for test
  !patch单元，第一列是邻接单元数目，后面是临界单元的单元编号
  INTEGER num_node_in_elem
  
  
  !build point patch information
  INTEGER, DIMENSION(:,:), POINTER :: tmp_big_patch,  tmp_big_patch_cg!patch内部所有的节点
  INTEGER, DIMENSION(:,:), POINTER :: temp_patch, point_patch_cg, point_patch_temp, ab_point_patch
  !INTEGER, DIMENSION(:,:), POINTER :: point_patch
  INTEGER, DIMENSION(:), POINTER   :: temp, temp_intrs, ppe_temp !temp是非界面非边界所用的临时数组，temp_intrs是去除界面单元的重复节点编号所用数组
  INTEGER  length   !temp,temp_intrs数组长度
  INTEGER  n_node_in_one_patch, n_node_in_one_patch_cg
  
  
  TYPE(ObjectType), DIMENSION(:), POINTER		::	objects
  REAL(8)                               :: x_loc,y_loc,   x_loc1,y_loc1
  INTEGER, DIMENSION(:), POINTER        :: ab_ppr_patch_index_temp, ab_ppr_patch_index
  REAL(8), DIMENSION(:,:), POINTER      :: ab_ppr_patch
  REAL(8), DIMENSION(:), POINTER        :: ab_ppr_point_x, ab_ppr_point_y, ab_ppr_point
  
  INTEGER                                         x
  
  
 
  !boundary情况
  INTEGER  n_boundary_node_re
  INTEGER, DIMENSION(:), ALLOCATABLE		    :: boundary_node_re, boundary_node_re_temp
  INTEGER                                   :: n_node_cg_re
  REAL(8), DIMENSION(:,:), ALLOCATABLE      :: P_order
  INTEGER                                      part_num, part_index
  INTEGER, DIMENSION(:), ALLOCATABLE		    :: repeat_times_ppr
  !INTEGER, DIMENSION(:,:), ALLOCATABLE  	  :: 
  REAL(8)                                   :: part_x, part_y
  INTEGER                                   :: part_ex_1, part_ex_2, part_ex_3, part_ex_4
  REAL(8)                                   :: xmin, xmax, ymin, ymax, distance
  INTEGER, DIMENSION(:), ALLOCATABLE		    :: ab_sample_node
  REAL(8), DIMENSION(:,:), ALLOCATABLE      :: ab_sample_cor
  REAL(8), DIMENSION(:), ALLOCATABLE        :: ab_point_x, ab_point_y
  INTEGER, DIMENSION(:,:), ALLOCATABLE      :: ab_point
  REAL(8), DIMENSION(:), ALLOCATABLE        :: ab_efx
  REAL(8), DIMENSION(:), ALLOCATABLE        :: ab_efy
  REAL(8)                                   :: x_cordi, y_cordi
  !INTEGER, DIMENSION(:), ALLOCATABLE      :: re_times
  INTEGER                                   :: re_nodes
  REAL(8)                                       r0
  

  
  n_node_sidg  = size(HP,2)
  num_node_in_elem = 4

  ALLOCATE(re_times(n_node_sidg))
  re_times=0
  re_nodes=0
  DO i=1,n_node_sidg
    
    DO j=i+1,n_node_sidg
      IF(re_times(j)==0) THEN
        IF(ABS(HP(1,i)-HP(1,j))<SmallValue .AND. ABS(HP(2,i)-HP(2,j))<SmallValue) THEN
          re_nodes=re_nodes+1
          re_times(j)=i
          !re_times(2,re_nodes)=j 
        ENDIF
      ENDIF
    ENDDO
      
  ENDDO
  WRITE(*,*) 'find adjacent element for all nodes,done'

  
!==========================construct point patch element =======================
!=========================================================================  
  ALLOCATE(point_patch_element(n_node_sidg,6),ppe_temp(4))
  point_patch_element(:,:) = 0
  ppe_temp(:) = 0
  DO i=1,n_node_sidg
    IF(re_times(i)==0)THEN 
      part_x = HP(1,i)
      part_y = HP(2,i)
      CALL ParPosition_ppr_elem(part_x, part_y, part_ex_1, part_ex_2, part_ex_3, part_ex_4)
      IF(DABS(part_x-xmin)< SmallValue .AND. DABS(part_y -ymin)< SmallValue) THEN!lower left
        
        point_patch_element(i,1) = 1
        point_patch_element(i,2) = 1
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_x-xmin)< SmallValue .AND. DABS(part_y -ymax)< SmallValue) THEN!upper left
        
        point_patch_element(i,1) = 1
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_x-xmax)< SmallValue .AND. DABS(part_y -ymin)< SmallValue) THEN!lower right
        
        point_patch_element(i,1) = 1
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_x-xmax)< SmallValue .AND. DABS(part_y -ymax)< SmallValue) THEN!upper right
        
        point_patch_element(i,1) = 1
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_x-xmin)< SmallValue) THEN!left boundary
        
        point_patch_element(i,1) = 2
        point_patch_element(i,2) = part_ex_2
        point_patch_element(i,3) = part_ex_3
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_y -ymax)< SmallValue) THEN!top boundary
        
        point_patch_element(i,1) = 2
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,3) = part_ex_2
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_x-xmax)< SmallValue) THEN!right boundary
        
        point_patch_element(i,1) = 2
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,3) = part_ex_4
        point_patch_element(i,6) = 2
      ELSEIF (DABS(part_y -ymin)< SmallValue) THEN!bottom boundary
        
        point_patch_element(i,1) = 2
        point_patch_element(i,2) = part_ex_4
        point_patch_element(i,3) = part_ex_3
        point_patch_element(i,6) = 2
    
      ELSEIF(element_index(part_ex_1)== In_index .AND.element_index(part_ex_2)== In_index.AND.&
                element_index(part_ex_3)== In_index .AND.element_index(part_ex_4)== In_index) THEN
      !inside 
        point_patch_element(i,1) = 4
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,3) = part_ex_2
        point_patch_element(i,4) = part_ex_3
        point_patch_element(i,5) = part_ex_4
        point_patch_element(i,6) = 1
      ELSEIF((element_index(part_ex_1)== Out_index .AND.element_index(part_ex_2)==Out_index.AND.&
                element_index(part_ex_3)==Out_index.AND.element_index(part_ex_4)==Out_index)) THEN
      !outside 
        point_patch_element(i,1) = 4
        point_patch_element(i,2) = part_ex_1
        point_patch_element(i,3) = part_ex_2
        point_patch_element(i,4) = part_ex_3
        point_patch_element(i,5) = part_ex_4
        point_patch_element(i,6) = 1
      
      ELSE!interface meshs
        point_patch_element(i,6) = 3
        IF(element_index(part_ex_1) < 0) THEN
          point_patch_element(i,1) = point_patch_element(i,1) + 1
          point_patch_element(i,2) = part_ex_1
        ENDIF
        IF(element_index(part_ex_2) < 0) THEN
          point_patch_element(i,1) = point_patch_element(i,1) + 1
          point_patch_element(i,3) = part_ex_2
        ENDIF
        IF(element_index(part_ex_3) < 0) THEN
          point_patch_element(i,1) = point_patch_element(i,1) + 1
          point_patch_element(i,4) = part_ex_3
        ENDIF
        IF(element_index(part_ex_4) < 0) THEN
          point_patch_element(i,1) = point_patch_element(i,1) + 1
          point_patch_element(i,5) = part_ex_4
        ENDIF
      ENDIF
     
       !delete repeated elements
      length=0
      DO j=1,4
        IF(point_patch_element(i,j+1)/=0)THEN
          length=length+1
          ppe_temp(length)=point_patch_element(i,j+1)
        ENDIF
      ENDDO
    
      ALLOCATE(temp(length))
      temp(1:length)=ppe_temp(1:length)
      CALL Bubble_Sort_new(temp,length)  
      point_patch_element(i,1:5)=0
      IF(length == 1)THEN
        point_patch_element(i,1)=1
        point_patch_element(i,2)=temp(1)
      ELSE
        DO x=1,length-1
          IF (temp(x)/=temp(x+1)) THEN
            point_patch_element(i,1)=point_patch_element(i,1)+1
            point_patch_element(i,point_patch_element(i,1)+1)=temp(x)
            IF(x==length-1) THEN
              point_patch_element(i,1)=point_patch_element(i,1)+1
              point_patch_element(i,point_patch_element(i,1)+1)=temp(x+1)
            ENDIF
          ELSEIF(x==length-1 .AND. temp(x)==temp(x+1) ) THEN
            point_patch_element(i,1)=point_patch_element(i,1)+1
            point_patch_element(i,point_patch_element(i,1)+1)=temp(x)
          ENDIF
        ENDDO
      ENDIF
    ELSE
      point_patch_element(i,:) = point_patch_element(re_times(i),:)
    ENDIF
    
  ENDDO
  
  ! wsy add for debug
  !allocate(point_patch_element_debug(6, size(point_patch_element, 1)))
  !Do i = 1, size(point_patch_element, 1)
  !      point_patch_element_debug(:, i) = point_patch_element(i, :)
  !End Do
  

  WRITE(*,*) 'find adjacent element for all nodes,done'
!=========================== point patch element======================== 
! write(*,*) point_patch_element(1560, :)
! Write(*,*) HP(:, 1560)
! write(*,*) element_index(1475), element_index(1500)
! write(*,*) HP(:,HT(1:4, 1500))
! stop
  
!===========================initial point patch================================ 
!============================================================================  
  n_node_in_one_patch = 50
  ALLOCATE(point_patch(n_node_sidg,n_node_in_one_patch))
  ALLOCATE(tmp_big_patch(n_node_sidg,n_node_in_one_patch))
  point_patch(:,:)=0
  tmp_big_patch(:,:)=0
  DO i=1,n_node_sidg
    k= point_patch_element(i,1)
    IF(k/=0) THEN!abnormal nodes with k=0
      IF(re_times(i)==0)THEN 
        DO j=2,k+1
          IF(point_patch(i,1)==0) THEN
            point_patch(i,1) = 4
            point_patch(i,2:5)=HT(1:4,point_patch_element(i,j))
          ELSE
            length = num_node_in_elem+point_patch(i,1)
            ALLOCATE(temp(length))
            temp(1:num_node_in_elem)=HT(1:4,point_patch_element(i,j))
            n=point_patch(i,1)
            temp((num_node_in_elem+1):length)=point_patch(i,2:(n+1))
            CALL Bubble_Sort_new(temp,length)
        
            point_patch(i,:)=0
            !delete repeated nodes
            DO x=1,length-1
              IF (temp(x)/=temp(x+1)) THEN
                point_patch(i,1)=point_patch(i,1)+1
                point_patch(i,point_patch(i,1)+1)=temp(x)
                IF(x==length-1) THEN
                  point_patch(i,1)=point_patch(i,1)+1
                  point_patch(i,point_patch(i,1)+1)=temp(x+1)
                ENDIF
              ELSEIF(x==length-1 .AND. temp(x)==temp(x+1) ) THEN
                point_patch(i,1)=point_patch(i,1)+1
                point_patch(i,point_patch(i,1)+1)=temp(x)
              ENDIF
            ENDDO
            DEALLOCATE(temp)
          ENDIF
        ENDDO
      ELSE
       point_patch(i,:) = point_patch(re_times(i),:)
      ENDIF
    ENDIF
  ENDDO
  
  WRITE(*,*) 'point patch for all nodes without extending, done'    
!===========================finish point patch================================                
 
  
!=================== abnormal nodes extend patch,result similiar to the relative coordinate method  ==================
  !ALLOCATE(ab_sample_node(16),ab_sample_cor(16,2),ab_point_x(n_node_sidg), ab_point_y(n_node_sidg), ab_point(n_node_sidg,16))
  !ALLOCATE(ab_point_patch(n_node_sidg,16))
  !ab_sample_node=0.D0
  !ab_sample_cor=0.D0
  !ab_point_x=0.D0
  !ab_point_y=0.D0
  !ab_point = 0.D0
  !ab_point_patch=0.D0
  !
  !DO i =1,n_node_sidg
  !  IF(point_patch(i,1)== 0)THEN
  !    x_loc = HP(1,i)!x_location
  !    y_loc = HP(2,i)
  !    CALL ParPosition_ppr_elem(x_loc, y_loc, part_ex_1, part_ex_2, part_ex_3, part_ex_4)
  !    ab_sample_node(1:4)   = HT(:,part_ex_1)
  !    ab_sample_node(5:8)   = HT(:,part_ex_2)
  !    ab_sample_node(9:12)  = HT(:,part_ex_3)
  !    ab_sample_node(13:16) = HT(:,part_ex_4)
  !    k=0
  !    IF (DSQRT(x_loc*x_loc+y_loc*y_loc) >= r0+SmallValue) THEN!外部
  !      DO j=1,16
  !        ab_sample_cor(j,:)=HP(:,ab_sample_node(j))
  !        IF(point_patch(ab_sample_node(j),1) /= 0)THEN
  !          IF(DSQRT(ab_sample_cor(j,1)*ab_sample_cor(j,1) + ab_sample_cor(j,2)*ab_sample_cor(j,2)) >= r0+SmallValue) THEN
  !            ab_point(i,k+1) = ab_sample_node(j)
  !            
  !            !ab_point_x(i) = ab_point_x(i) + efx_1(ab_sample_node(j),1)
  !            !ab_point_y(i) = ab_point_y(i) + efy_1(ab_sample_node(j),1)
  !            k=k+1
  !          ENDIF
  !        ENDIF
  !      ENDDO
  !    ELSE!内部
  !      DO j=1,16
  !        ab_sample_cor(j,:)=HP(:,ab_sample_node(j))
  !        IF(point_patch(ab_sample_node(j),1) /= 0)THEN
  !          IF(DSQRT(ab_sample_cor(j,1)*ab_sample_cor(j,1) + ab_sample_cor(j,2)*ab_sample_cor(j,2)) < r0+SmallValue) THEN!同为内部
  !            ab_point(i,k+1) = ab_sample_node(j)
  !            !ab_point_x(i) = ab_point_x(i) + efx_1(ab_sample_node(j),1)
  !            !ab_point_y(i) = ab_point_y(i) + efy_1(ab_sample_node(j),1)
  !            k=k+1
  !          ENDIF
  !        ENDIF
  !      ENDDO
  !    ENDIF
  !    !ab_point_x(i)=ab_point_x(i)/k
  !    !ab_point_y(i)=ab_point_y(i)/k
  !    !efx_1(i,1)=ab_point_x(i)
  !    !efy_1(i,1)=ab_point_y(i)
  !    DO j=1,16
  !      IF(ab_point(i,j)/=0)THEN
  !        ab_point_patch(i,1)=ab_point_patch(i,1)+1
  !        ab_point_patch(i,j+1)=ab_point(i,j)
  !      ENDIF
  !    ENDDO
  !    
  !  ENDIF
  !  
  !ENDDO
  !
  !DO i =1,n_node_sidg
  !  IF(ab_point_patch(i,1)/=0 )THEN
  !  DO j=1,16
  !    point_patch(i,j)=ab_point_patch(i,j)
  !  ENDDO
  !  ENDIF
  !ENDDO
  !
  !WRITE(*,*) 'abnormal nodes treatment, done'
  
!===========================================================================  
   DO i=1,n_node_in_one_patch
    DO j=1,n_node_sidg
        tmp_big_patch(j,i) = point_patch(j,i)
    ENDDO
  ENDDO
  
!========================sampling nodes <6 ,not enough to calculate ==============================
!===========================================================================
  xx = 0; yy = 0
  ALLOCATE(temp_patch(n_node_sidg,n_node_in_one_patch))
  temp_patch = 0
  DO i=1,n_node_sidg
      IF(point_patch(i,1) < 6 .AND. point_patch(i,1)/=0)THEN
        IF(re_times(i)==0)THEN 
          DO j=2,tmp_big_patch(i,1)+1
            IF(temp_patch(i,1) ==0)THEN
                k = tmp_big_patch(i,j)
                xx = tmp_big_patch(k,1)
                temp_patch(i,1) = xx
                temp_patch(i,2:xx+1) = tmp_big_patch(k ,2:xx+1)    
            ELSE
                k = tmp_big_patch(i,j)
                xx = tmp_big_patch(k,1)
                length = xx + temp_patch(i,1)
                !WRITE(6,*) 'length',length
                !PAUSE
                ALLOCATE(temp(length))
                temp(1:xx) = tmp_big_patch(k,2:xx+1)
                temp(xx+1:length) = temp_patch(i,2:temp_patch(i,1)+1)

                CALL Bubble_Sort_new(temp, length)
                temp_patch(i,:) = 0                             
                DO x=1, length-1
                    IF(temp(x) /= temp(x+1))THEN
                        temp_patch(i,1) = temp_patch(i,1) + 1
                        yy = temp_patch(i,1)+1
                        temp_patch(i,yy) = temp(x)
                                  
                        IF(x == length-1)THEN
                            temp_patch(i,1) = temp_patch(i,1) + 1
                            yy = temp_patch(i,1)+1
                            temp_patch(i,yy) = temp(x+1)    
                        ENDIF
                    ELSEIF(x == length-1 .AND. temp(x) == temp(x+1))THEN
                        temp_patch(i,1) = temp_patch(i,1) + 1
                        yy = temp_patch(i,1)+1
                        temp_patch(i,yy) = temp(x)
                    ENDIF
                ENDDO
                DEALLOCATE(temp)      
            ENDIF
          ENDDO
        ELSE
          temp_patch(i,:) = temp_patch(re_times(i),:) 
        ENDIF
      ENDIF
    !ENDIF
  ENDDO
  WRITE(*,*) 'extend layers for nodes with less than 6 samples, done' 
!========================finish extending patch for nodes with less than 6 samples=====
  
  DO i=1,n_node_sidg                                    ! 将temp_patch更新到point_patch
    IF(temp_patch(i,1)/=0)THEN 
        DO num=1,temp_patch(i,1)+1
          point_patch(i,num) = temp_patch(i,num)
        ENDDO
    ENDIF   
  ENDDO
  !DEALLOCATE(temp_patch)
  !DO i=1,n_node_sidg
  !  DO j=1,n_node_in_one_patch
  !      tmp_big_patch(i,j) = point_patch(i,j)
  !  ENDDO
  !ENDDO
  
!===========================boundary treatment （include interface）====================================
!=====================================================================================
  WRITE(*,*) 'ready for boundary treatment'
  ALLOCATE(temp_patch(n_node_sidg,n_node_in_one_patch))
  temp_patch(:,:)=0
  DO i=1,n_node_sidg
    IF(point_patch_element(i,6)==2 .OR. point_patch_element(i,6)==3)THEN!to find boundary 
      IF(point_patch(i,1)<12) THEN!can change the results===========================
        IF(re_times(i)==0) THEN
          DO j=2, tmp_big_patch(i,1)+1
            IF(temp_patch(i,1)==0) THEN
              k=tmp_big_patch(i,j)
              xx=tmp_big_patch(k,1)!
              temp_patch(i,1) = xx!
              temp_patch(i,2:xx+1) = tmp_big_patch(k,2:xx+1)
            ELSE
              k=tmp_big_patch(i,j)
              xx=tmp_big_patch(k,1)
              length=xx+temp_patch(i,1)
              IF(length > n_node_in_one_patch) THEN
                !write(*,*) 'n_node_in_one_patch is', n_node_in_one_patch
                !write(*,*) 'It is not big enough. STOP'
                n_node_in_one_patch = n_node_in_one_patch+8
                !DEALLOCATE(temp_patch)
                !goto 10
                !
                !PAUSE
              ENDIF
              
              ALLOCATE(temp(length))
              temp(1:xx)=tmp_big_patch(k,2:xx+1)
              temp(xx+1:length)=temp_patch(i,2:temp_patch(i,1)+1)
          
              CALL Bubble_Sort_new(temp,length)
              temp_patch(i,:)=0
              
              DO x=1,length-1
                IF(temp(x) /= temp(x+1))THEN
                  temp_patch(i,1) = temp_patch(i,1) + 1
                  yy = temp_patch(i,1)+1
                  temp_patch(i,yy) = temp(x)
                  IF(x == length-1)THEN
                    temp_patch(i,1) = temp_patch(i,1) + 1
                    yy = temp_patch(i,1)+1                  
                    temp_patch(i,yy) = temp(x+1)  
                  ENDIF
                ELSEIF(x == length-1 .AND. temp(x) == temp(x+1))THEN
                  temp_patch(i,1) = temp_patch(i,1) + 1
                  yy = temp_patch(i,1)+1
                  temp_patch(i,yy) = temp(x)
                ENDIF
              ENDDO
              DEALLOCATE(temp)
            ENDIF
          ENDDO
        ELSE
          temp_patch(i,:) = temp_patch(re_times(i),:) 
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  WRITE(*,*) 'boundary treatment, done'  
  
  
  
  
!===========================finish boundary treatment ====================================
  DO i=1,n_node_sidg                                    ! 将temp_patch更新到point_patch
    IF(temp_patch(i,1)/=0)THEN 
        DO num=1,temp_patch(i,1)+1
          point_patch(i,num) = temp_patch(i,num)
        ENDDO
    ENDIF   
  ENDDO
  WRITE(*,*) 'ready to recover'
  !=====================================abnormal treatment================================================= 
  
  
  DO i =1,n_node_sidg
    IF(point_patch(i,1)== 0)THEN
      !find nodes without patch
      ALLOCATE(ab_sample_node(16), ab_point_patch(n_node_sidg,17))
      ab_sample_node=0.D0
      ab_point_patch= 0.D0
      exit
    ENDIF
  ENDDO
  
  IF (ALLOCATED(ab_sample_node)) THEN
    DO i =1,n_node_sidg
      IF(point_patch(i,1)== 0)THEN    
        x_loc = HP(1,i)!x_location
        y_loc = HP(2,i)
        CALL ParPosition_ppr_elem(x_loc, y_loc, part_ex_1, part_ex_2, part_ex_3, part_ex_4)
        ab_sample_node(1:4)   = HT(:,part_ex_1)
        ab_sample_node(5:8)   = HT(:,part_ex_2)
        ab_sample_node(9:12)  = HT(:,part_ex_3)
        ab_sample_node(13:16) = HT(:,part_ex_4)
        k=0
        IF(node_index(i) == Out_index) THEN!abnormal node is outside
          DO j=1,size(ab_sample_node)
            IF(point_patch(ab_sample_node(j),1) /= 0)THEN!remove other abnormal nodes
              IF(node_index(ab_sample_node(j)) == Out_index)THEN!adapt outside nodes around 
                ab_point_patch(i,1) = ab_point_patch(i,1) + 1
                ab_point_patch(i,k+2) = ab_sample_node(j)
                k=k+1
              ENDIF
            ENDIF
          ENDDO
        ELSE!abnormal node is inside
          DO j=1,size(ab_sample_node)
            IF(point_patch(ab_sample_node(j),1) /= 0)THEN
              IF(node_index(ab_sample_node(j)) == In_index)THEN  
                ab_point_patch(i,1) = ab_point_patch(i,1) + 1
                ab_point_patch(i,k+2) = ab_sample_node(j)
                k=k+1
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      
        length = ab_point_patch(i,1)
        ALLOCATE(temp(length))
        temp(1:length)=ab_point_patch(i,2:length+1)
        ab_point_patch(i,:)=0
        CALL Bubble_Sort_new(temp,length)
        DO x=1,length-1
          IF (temp(x)/=temp(x+1)) THEN
            ab_point_patch(i,1)=ab_point_patch(i,1)+1
            ab_point_patch(i,ab_point_patch(i,1)+1)=temp(x)
            IF(x==length-1) THEN
              ab_point_patch(i,1)=ab_point_patch(i,1)+1
              ab_point_patch(i,ab_point_patch(i,1)+1)=temp(x+1)
            ENDIF
          ELSEIF(x==length-1 .AND. temp(x)==temp(x+1) ) THEN
            ab_point_patch(i,1)=ab_point_patch(i,1)+1
            ab_point_patch(i,ab_point_patch(i,1)+1)=temp(x)
          ENDIF
        ENDDO
        
        !DO j=2,ab_point_patch(i,1)+1
        !  point_patch(i,j)=ab_point_patch(i,j)
        !ENDDO
        
      ENDIF
    ENDDO
  ENDIF
  !Do i = 1, n_node_sidg
  !  Write(*,*) 'node = ', i
  !  Write(*,*) point_patch(i, 1:18)
  !End Do
  !stop

  ! check the patch ensure the number of sampling node not equal zero
!Do i = 1, n_node_sidg
!    If (point_patch(i,1) == 0) then
!        Write(*,*) i, '  node dont have sampling node pls check'
!        Write(*,*) 'error happen in Get_PointPatch'
!        stop
!    End IF
!End Do
  
ALLOCATE(efx_1(n_node_sidg,1), efy_1(n_node_sidg,1))

Allocate(PPRPointIndex(n_node_sidg))

PPRPointIndex(:) = 0 

Do i = 1, size(element_index,1)
    If (element_index(i) > 0) then
        Do j = 1, 4
            PPRPointIndex(HT(j, i)) = 1
        End Do
    End If
End Do
PPRPointIndex = 1

    
END SUBROUTINE Get_PointPatch
