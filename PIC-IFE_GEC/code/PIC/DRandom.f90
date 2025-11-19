SUBROUTINE DRandom(randum)
    
! this is a version of the random number generator dprandom due to
! c. bingham and the yale computer center, producing numbers
! in the interval (0,1).  written for the sun by viktor k. decyk, ucla

!!! ************************ bjw add for impic 2019-6-3 **********************************************
!REAL(8)       ::  ranum
!real:: lbound,ubound,leng,t,t1,t2
!real:: seed
!!INTEGER :: time
!
!!CALL CPU_TIME(t1)
!lbound = 0.0
!ubound = 1.0
!leng=ubound-lbound  !计算范围大小
!
!!call random_seed()!库函数，生成随机数前调用
!call random_number(t)  !t是0-1之间的随机数
!
!!call random_seed(time)  !根据系统时间生存seed
!!call date_and_time(value=time)
!!seed=time(4)*(360000*time(5)+6000*time(6)+100*time(7)+time(8))
!!call random_seed(put=seed)
!!call random_number(t)  !t是0-1之间的随机数
!
!randum=lbound+leng*t
!!CALL CPU_TIME(t2)
!!WRITE(*,*) t2-t1
!!WRITE(*,*) randum
!return
!!! ************************ bjw add for impic 2019-6-3 **********************************************



integer r1,r2
double precision randum,h1l,h1u,r0,r3,asc,bsc
save r1,r2,h1l,h1u
data r1,r2 /1271199957,1013501921/
data h1l,h1u /65533.0d0,32767.0d0/
isc = 65536
asc = dble(isc)
bsc = asc*asc
i1 = r1 - (r1/isc)*isc
r3 = h1l*dble(r1) + asc*h1u*dble(i1)
i1 = r3/bsc
r3 = r3 - dble(i1)*bsc
bsc = 0.5d0*bsc
i1 = r2/isc
isc = r2 - i1*isc
r0 = h1l*dble(r2) + asc*h1u*dble(isc)
asc = 1.0d0/bsc
isc = r0*asc
r2 = r0 - dble(isc)*bsc
r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
isc = r3*asc
r1 = r3 - dble(isc)*bsc
randum = (dble(r1) + dble(r2)*asc)*asc
return


END   