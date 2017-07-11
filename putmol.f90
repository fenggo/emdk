!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@foxmail.com                      =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!         把因为周期性边界条件而被分裂的分子再放到一块                      -
!-------------------------------------------------------------------------
!
      subroutine putmol
      use mode  , only: npart
      use modfim, only: cflag
      !use modctrl,only: irec
      implicit none
!
      integer           i,iatom
!  
!
     print *,'Recovering moles ...'
      do i=1,npart
         cflag(i)=1
      enddo

      do i=1,npart
!
         iatom = i

         if(cflag(iatom)==1)then
                cflag(iatom) = 0
                call pull_atom(iatom)
!                print *,'!!Debug!!'
         endif
      enddo
!
!
     return
     end
!
!
!-------------------------------------------------------------------------
!  ------             RECURSIVE pull_atom             ------             -
!-------------------------------------------------------------------------
!
   
    recursive subroutine pull_atom(m)
    use mode,   only: types,x,a,b,c
    use modnei, only: atable,matom,rcutsq
    use modfim, only: cflag
    implicit none
!
    integer      i,ii,m
    integer      iatom,jatom
    integer      itype,jtype
    integer      ina,inb,inc
    real(8)      delx,dely,delz
    real(8)      rsq
    real(8)      r_a,r_b,r_c
!
!
!    if(cflag(m)/=0)return
!
    iatom= m
    ii = atable(1,iatom) + 1

    itype= types(iatom)
!    print *,itype
    cflag(iatom) = 0

    if(ii>1)then
    do i=2,ii
       jatom = atable(i,iatom)
       jtype = types(jatom)
!       print *, jtype
       delx=x(1,jatom)-x(1,iatom)
       dely=x(2,jatom)-x(2,iatom)
       delz=x(3,jatom)-x(3,iatom)
!
       rsq=delx*delx+dely*dely+delz*delz

       if(rsq>rcutsq(itype,jtype))then
             r_c = (a(1)*b(2)-a(2)*b(1))*(b(2)*delz-b(3)*dely) - (a(3)*b(2)-a(2)*b(3))*(b(2)*delx-b(1)*dely)
             r_c = r_c/((a(1)*b(2)-a(2)*b(1))*(c(3)*b(2)-c(2)*b(3)) - (a(3)*b(2)-a(2)*b(3))*(c(1)*b(2)-c(2)*b(1)))

             r_a = b(2)*delx-b(1)*dely - (c(1)*b(2)-c(2)*b(1))*r_c
             r_a = r_a/(a(1)*b(2)-a(2)*b(1))

             r_b = dely- r_a*a(2) - r_c*c(2)
             r_b = r_b/b(2)

             ina = anint(r_a)
             inb = anint(r_b)
             inc = anint(r_c)

             x(1,jatom)= x(1,jatom) - ina*a(1) - inb*b(1) - inc*c(1)
             x(2,jatom)= x(2,jatom) - ina*a(2) - inb*b(2) - inc*c(2)
             x(3,jatom)= x(3,jatom) - ina*a(3) - inb*b(3) - inc*c(3)
       endif

       if(cflag(jatom)==1)call pull_atom(jatom)

    enddo
    endif
!
!
    return
    end
!
!-------------------------------------------------------------------------
!

