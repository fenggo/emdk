!
!                           WRITED BY GUO F.                             
!                                                                        
!              Institute of Atomic & Molecular Physics of Si-            
!                       chuan University of China.                       
!                                                                        
!                     email: gfeng.alan@foxmail.com                      
!
!
!
!
!
!  develop the slab into two part without split molecules

subroutine develop
    use mode,only: xd,npart,v,vel
    use modfim,only: n_mol, mlist,molelist
    implicit none
 
    integer      j,jj,l,atom
    integer::    larg = 0
    !real(8)::    vel = 8.0
    !logical lwrite

    open(10,file='velocity.txt',status='unknown')
    write(10,*)
    write(10,'(A10)')'Velocities'
    write(10,*)

    !print *,'Please input the velocity:'
    !read(5,*) vel

    call neighbor
    call findmole
    call proj_coord

    do j =1, n_mol
       if(j>1)then
           jj=mlist(j-1) + 1
       else
           jj=1
       endif

       do l=jj, mlist(j) 
          atom = molelist(l)
          if(xd(1,atom)>0.5) larg = 1
       enddo
       do l=jj, mlist(j) 
          atom = molelist(l)
          if(larg == 1)then
                  v(1,atom)= -vel
          elseif(larg == 0)then
                  v(1,atom)= vel
          endif
       enddo
    enddo

    v(2,:) = 0.0
    v(3,:) = 0.0

    do j=1,npart
       write(10,*)v(:,j)
    enddo

    close(10)
    return
end



