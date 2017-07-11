!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@hotmail.com                      =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!  ------             ZERO DATA                 ------                   -
!-------------------------------------------------------------------------
!

subroutine zerodata
   use mode, only: npart,xd,a,b,c,r1,r2,r3,alpha,beta,gamm
   use modctrl,only: ilt
   implicit none
!
   integer i,j
   real(8) center(3)
!
!   real,parameter::a0=0.529117
!
   !real(8)      cubex,cubey,cubez
   real(8)      mol(3),mol1(3),mol2(3),mol3(3)
   real(8)      cy,ay
   real(8)      pi,alph,bet,gam

   real(8),external:: compare_min
   real(8),external:: compare_max

   pi = dacos(-1.0d0)
!
   print *, 'Centering coordinate ...'
   call proj_coord
!
   do i=1,npart
!
       if(i==1)then
          do j=1,3
             mol(j)=xd(j,i)
             mol2(j)=xd(j,i)
          enddo
       elseif(i>1)then
!
         do j=1,3
           mol1(j)=xd(j,i)
           mol3(j)=xd(j,i)
!
           mol(j)  = compare_min(mol(j),mol1(j))
           mol2(j) = compare_max(mol2(j),mol3(j))
!
         enddo
       endif
!
   enddo
!
!   write(*,*)'max boundery(Used by CPMD):'
!
!   cubex=(mol2(1)-mol(1))!/a0
!   cubey=(mol2(2)-mol(2))!/a0
!   cubez=(mol2(3)-mol(3))!/a0
!   print '(2x,A20)','Max Atom positions:'
!   print '(2X,A)',' '
!   write(*,500)'MAX IN X',cubex
!   write(*,500)'MAX IN Y',cubey
!   write(*,500)'MAX IN Z',cubez
!   print !'(A80)','--------------------------------------------------------------------------------'
!
!  used by CPMD
!
!   write(*,*)a,b/a,c/a
!
!
  center(1)= 0.5*(mol2(1)+mol(1))
  center(2)= 0.5*(mol2(2)+mol(2))
  center(3)= 0.5*(mol2(3)+mol(3))


  do i=1,npart
     xd(1,i)=xd(1,i) - center(1) +0.5
     xd(2,i)=xd(2,i) - center(2) +0.5
     xd(3,i)=xd(3,i) - center(3) +0.5
  enddo
  alph = alpha*pi/180.0d0
  bet = beta*pi/180.0d0
  gam= gamm*pi/180.0d0
! a in X axis, b in XY plane

    if(ilt==1)then
!           print '(2x,A)','Defult: a in X axis, b in XY Plane, c arbitary!'
           a = (/r1,0.0d0,0.0d0/)
           b = (/r2*dcos(gam),r2*dsin(gam),0.0d0/)
           cy= r3*dcos(alph)-r3*dcos(bet)*dcos(gam)
           cy= cy/dsin(gam)
           c = (/r3*dcos(bet),cy,dsqrt(r3*r3*dsin(bet)*dsin(bet)-cy*cy)/)
    elseif(ilt==2)then
!  c along Z, b in YZ plane
!           print '(2x,A)','For Reax geo: c along Z, b in YZ plane, a arbitary!'
           c = (/0.0d0,0.0d0,r3/)
           b = (/0.0d0,r2*dsin(alph),r2*dcos(alph)/) 
           ay= r1*dcos(gam)-r1*dcos(alph)*dcos(bet)
           ay= ay/dsin(alph)
           a = (/dsqrt(r1*r1*dsin(bet)*dsin(bet)-ay*ay),ay,r1*dcos(bet)/)
    elseif(ilt==0)then
!  c along Z, b in YZ plane
!           print '(2x,A)','For Reax geo: c along Z, b in YZ plane, a arbitary!'
           c = (/0.0d0,0.0d0,r3/)
           b = (/0.0d0,r2*dsin(alph),r2*dcos(alph)/) 
           !print *,'alpha',alph,sin(alph)
           ay= r1*dcos(gam)-r1*dcos(alph)*dcos(bet)
           ay= ay/dsin(alph)
           a = (/dsqrt(r1*r1*dsin(bet)*dsin(bet)-ay*ay),ay,r1*dcos(bet)/)

    endif

!  print '(12X,A50)','**********    Lattice vector after rotate     **********'

!  print '(3F12.5)'    ,(a(i),i=1,3)
!  print '(3F12.5)'    ,(b(i),i=1,3)
!  print '(3F12.5)'    ,(c(i),i=1,3)

  call cart_coord
!
!
 return
!

end
!
!**********************************
!        subroutine  max,min      *
!**********************************
!

function compare_min(x,y)
    implicit none
    real(8) x,y
    real(8) compare_min
    if(x>y)then
      compare_min = y
    else
      compare_min = x
    endif

    return
end
!

function compare_max(x,y)
   implicit none
   real(8) x,y
   real(8) compare_max
   if(x<y)then
      compare_max = y
   else
      compare_max = x
   endif
   return
end
!
