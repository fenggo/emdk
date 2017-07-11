!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                       email: gfengs@hotmail.com                        =
!=========================================================================
!
!
! analysis neighborring atoms
!
   subroutine  analysis
   use modnei, only: atable,matom
   use mode,   only: x,types,npart
   implicit none

   integer                i,j
   integer                iatom,jatom
   integer                itype
   integer, allocatable:: fst(:)
   integer, allocatable:: snd(:)
   integer                fls,sls
   integer                per
   real(8)                delx,dely,delz
   real(8)                rsq,radi
   real(8)                rt1,rt2
   real(8)                re1,re2
!
   allocate(fst(npart))
   allocate(snd(npart))
!
   re1 =  3.1650
   re2 =  sqrt(3.0)*re1/2
!
!
   do i=1,npart
      iatom = i
      itype=types(i)

      fst(i) = 0
      snd(i) = 0

      do j=2,matom+1

         jatom= atable(j,iatom)
         if(jatom==0)exit  !exit cycles when there no atom.

         delx=x(1,iatom)-x(1,jatom)
         !delx=delx-cube(1)*anint(delx/cube(1))
         dely=x(2,iatom)-x(2,jatom)
         !dely=dely-cube(2)*anint(dely/cube(2))
         delz=x(3,iatom)-x(3,jatom)
         !delz=delz-cube(3)*anint(delz/cube(3))
!
         rsq=delx*delx+dely*dely+delz*delz
         radi=sqrt(rsq)

! First near neighbor or Second.
         if(radi>=re1)then
            snd(i)=snd(i)+1
         elseif(radi<=re2)then
            fst(i)=fst(i)+1
         else
            rt1=radi/re2
            rt2=re1/radi
            if(rt1<=rt2)then
              fst(i)=fst(i)+1
            else
              snd(i)=snd(i)+1
            end if
         end if

      end do
   end do

   open(33,file='NEICAR',status='unknown',access='append')
   write(33,100)'Number of 1st and 2nd Neighbors of every atom'

   fls = 0
   sls = 0
   per = 0

   do i=1,npart
      write(33,200)i,fst(i),snd(i)
      if(fst(i)<8)fls=fls+8-fst(i)
      if(snd(i)<6)sls=sls+6-snd(i)
      if(fst(i)==8.and.snd(i)==6)per=per+1
   end do
      write(33,100)'1st and 2nd neighbor lost'
      write(33,300)fls,sls
      write(33,100)'Atoms with no neighbor lost'
      write(33,310)per

100   format(2x,A)
200   format(3I5)
300   format(2I5)
310   format(I5)

   close(33)
   return
   end
!
