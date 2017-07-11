!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@foxmail.com                     =
!=========================================================================
!
!

subroutine atm_typ
    use modfim,only: n_mol, mlist, molelist
    use modnei,only: mang,mdih
    use mode,  only: npart,types,ctype!,q !, a, b, c ,alive ,x
    use mod_nw,only: atom_type
    implicit none

    integer         mstart, mend
    integer         i,j,jj,ii 
    integer         n_c,n_h,n_n,n_o
    character(7)    tc
    real(8),external  :: torsion_angle
    !logical           :: dihe_vertical = .false.
    logical           :: molnr = .false.

    if(.not. allocated(atom_type)) allocate(atom_type(npart))
    atom_type  = 'null '

!@@@@@@@@@@@@@@@@       topology analysis      @@@@@@@@@@@@@@@@
    !!! determine atom type
    
    !print '(A,I)','Molecule: ',n_mol
    
    do ii=1,n_mol
       if(ii>1)then      ! n molecule ii list start
         mstart=mlist(ii-1) + 1
       else
         mstart=1
       endif
       mend = mlist(ii)
       n_c = 0
       n_h = 0
       n_n = 0
       n_o = 0
       do jj=mstart, mend ! find out molecule types
          i = molelist(jj)
          tc = ctype(types(i))
          if (tc=='H') n_h = n_h + 1
          if (tc=='C') n_c = n_c + 1
          if (tc=='N') n_n = n_n + 1
          if (tc=='O') n_o = n_o + 1
       enddo 
      
       !print '(A,I)','Molecule: ',ii

       if (n_h==6 .and. n_c==6 .and. n_n==12 .and. n_o==12)then 
          call atm_typ_cl20(mstart,mend)
          !print '(I,A)',ii, ' CL20'
       elseif(n_h==8 .and. n_c==4 .and. n_n==8 .and. n_o==8)then
          call atm_typ_hmx(mstart,mend)
          !print '(I,A)',ii, ' HMX'
       elseif(n_h==3 .and. n_c==1 .and. n_n==1 .and. n_o==2)then
          call atm_typ_nm(mstart,mend)
       else
          print '(A)','Warning: Molecule not recognized!'
          !print '(I6,4I3)',n_mol, n_c, n_h, n_o, n_n
          !call writeat
          molnr = .true.
       endif
    enddo

    if(molnr)then
        do j=1,npart
           atom_type(j) = ctype(types(j))
        enddo
    endif
    
    return
end


subroutine write_species
use mode,  only: x,npart
use mod_nw,only: atom_type
implicit none

   integer i


   call atm_typ

   open(36,file='species.txt',status='unknown')

   do i=1,npart
      write(36,'(A5,3F13.8)')atom_type(i),x(1,i),x(2,i),x(3,i)
   enddo

   close(36)

return
end
