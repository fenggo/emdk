!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!  ------              WRITE file for tinker       ------                -
!-------------------------------------------------------------------------
!

subroutine tinker_ff
    use mode,only: x,npart,types,ctype
    use modnei,only: atable
    !!use mod_nw,only: atom_type
    implicit none
!
!   Local Variables
!
    integer               i,j,jj
    integer               n,k
    integer               n_c,n_h,n_n,n_o
    integer,allocatable:: atom_type(:)
    character(2)          tc   
    character(6)          cc
    character(80)         line
    real(8)               pi
    logical           ::  dihe_zero

    pi = dacos(-1.0d0)
!
    allocate(atom_type(npart))

    call neighbor
    call findmole

!@@@@@@@@@@@@@@@@       topology analysis      @@@@@@@@@@@@@@@@
    !!! determine atom type  ! MMFF parameters
    do i=1,npart
       if(ctype(types(i))=='C') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          if(n_h==2) atom_type(i) = 1 ! Carbon with two H sp3
       endif
    enddo

    do i=1,npart
       if(ctype(types(i))=='H') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          if(n_c==1) atom_type(i) = 28
       elseif(ctype(types(i))=='N') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          dihe_zero = .false.
          if(n_c==2.and.n_n==1) atom_type(i) = 54  ! N sp3

          if(n_o==2.and.atable(1,i)==3) then  !######## N in nitro-group
           
                  atom_type(i) = 45 ! N sp3 Nitro-group N_R1
          endif                               !########## end analysis N in nitro-group
       elseif(ctype(types(i))=='O') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          if(n_n==1.and.atable(1,i)==1) atom_type(i) = 32 ! sp2 Oxygen
          !!Oxygen in Nitro-group
       endif
    enddo

!@@@@@@@@@@@@@@@@       topology analysis      @@@@@@@@@@@@@@@@


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    open(30,file='card_tinker.xyz',status='unknown')
    !print *,'problem here ...'
    !stop
    write(30,*)npart
    !print *,npart
    !write(30,'(A4,3F12.6,3F12.6,A5)')'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    !print *,'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    do i=1,npart
          jj = atable(1,i)+1
          !print *,jj            problem !!!!! atable
          !if(alive(i)/=0) write(30,*)i,ctype(types(i)),x(1,i),x(2,i),x(3,i), &
            !& atom_type(i), (atable(j,i),j=2,jj)
          write(line,'(I5,A3,3F12.7,I4)')i,ctype(types(i)),x(1,i),x(2,i),x(3,i),atom_type(i)
          do j=2,jj  
             !line = trim(adjustl(line))
             k=len(trim(line)) + 1
             write(cc,'(I6)') atable(j,i)
             cc = trim(adjustl(cc))
             line(k+1:K+3)= cc
          enddo
          write(30,'(A)') line
    enddo


    deallocate(atom_type)
    close(30)

    return
end
