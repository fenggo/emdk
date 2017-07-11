!=========================================================================
!                           WRITED BY GUO F.                     =
!                                                                                                                                                 =
!              Institute of Atomic & Molecular Physics of Si-    =
!                       chuan University of China.               =
!                                                                                                                                                 =
!                     email: gfeng.alan@gmail.com                =
!=========================================================================
!
!
! 
!----    shock velocity for molecular (symmetic impact)   ----

subroutine velocity
use mode,   only: xd,v,types,mass,shockv
use modfim, only: molelist,mlist,n_mol
implicit none

    integer     i,iatom,j,jj
    real(8)     mass1,mass2
    real(8) ::  kms=0.0100
    !real(8)     shockv

    call neighbor
    call findmole
    call putmol
    call proj_coord
        
    shockv = shockv*kms
    ! v = 0.0  ! initialize velocity
    mass1 = 0.0
    mass2 = 0.0

    do i=1,n_mol  ! dertermin mass of the two slab

       if(i>1)then
            jj=mlist(i-1)+1
       else
            jj=1
       endif

       do j=jj, mlist(i)
          iatom = molelist(j)
          if (xd(1,molelist(jj))<0.5) then
             mass1 = mass1 + mass(types(iatom))
          else
             mass2 = mass2 + mass(types(iatom))
          endif

       enddo

    enddo

    !print *,mass1,mass2

    do i=1,n_mol

       if(i>1)then
            jj=mlist(i-1)+1
       else
            jj=1
       endif

       do j=jj, mlist(i)
          iatom = molelist(j)
          if (xd(1,molelist(jj))<0.5) then
             v(1,iatom) = v(1,iatom) + shockv
          else
             v(1,iatom) = v(1,iatom)-shockv*mass1/mass2 
          endif

       enddo

    enddo

end subroutine


