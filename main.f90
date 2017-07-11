!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!                     LiaoChen University of China.                      =
!                                                                        =
!                     email: gfeng.alan@foxmail.com                      =
!=========================================================================
!
!
!
!-------------------------------------------------------------------------
!  ------            Effective Material Design Kit  ------               -
!-------------------------------------------------------------------------
!
!

program main
   use mode,    only: movex,movey,movez,alive
   use modctrl, only: iwr,isc,center,irec,istretch,fzoom,lmap,lcut, &
                fremove,ldipole,lmole, istretchfreebond
   use mod_def, only: idef
   !use modnei, only: 
   use modfim, only: mlist,molelist
   use mod_rotate,only: irotate,iswing
   use mod_grd,only: grid_op
   implicit none
!
   integer      i
!   character    iquit
!
   call control
!
   if(isc)call supercell
   
   call write_lp
   call proj_coord
   !print *,'problem here ...'


   !print *,'problem lp ...'
!
   if(iwr==3.or.iwr==4)call neighbor

   if(irec)then
      call neighbor
      call findmole
      call putmol
      if (lmole) then
          alive = 0
          !print *,mlist(1)
          do i=1,mlist(1)
             alive(molelist(i))= 1
             !print *,molelist(i)
          enddo
      endif
   endif
!
   if(center)call zerodata
!

   if(fzoom)call zoom
   if(idef>0) call defects

   if(lcut) call develop

   if(grid_op) call grid

   if(movex/=0.0.or.movey/=0.0.or.movez/=0.0) call move
   
   if (fremove) then
        call neighbor
        call findmole
        call remove
   endif

   if (ldipole) then
        call neighbor
        call findmole
        call dipole        ! compute the dipole moments of every molecules~
   endif

   call writeout
   !print *,'problem here ...'
   !stop
   if(istretch)call stretch_config
   if(istretchfreebond) call stretch_free_bond
   !print *,istretchfreebond,istretch

   if(lmap) call write_map 

   if(irotate) call rotate
   if(iswing)  call swing
!
!  Release Memory
   call release
!
!   print '(A)','Press any key and Enter to exit ...'
!   read(*,*)iquit

!
!
   stop
end
!
!-------------------------------------------------------------------------
