!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                       email: gfengs@foxmail.com                        =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!-----   compute Radia Pair Distribution Function                    -----
!-------------------------------------------------------------------------
!
     subroutine gofr
     use mode, only: types,ctype,npart,x
     use modgr
     implicit none
!
!----  parameters  ----
     real(8),parameter:: bin=0.01
     real(8),parameter:: pi=3.1415926535897932
!     real(8),parameter:: avo=6.022
!
     integer nbin
     integer i,j
     integer itype,jtype
     !integer nt !, id
     integer natm1,natm2
!
!---- icom1,2 atom types to be computed ----
     integer uni
!
     real(8) rdx,rdz,rdy,rd,rsq
!
     real(8) g,tem
!
!
     integer ind
     real(8),allocatable:: gr(:)
!
     real(8) volume,rou
!
     integer sflag
!
     print 500,'Computing g(r) ...'
!
     open(50,file='GOFR.log',status='unknown',access='append')
!
!
!
     ntype(1)=0
     ntype(2)=0
     ntype(3)=0
     ntype(4)=0
!
!
     do i=1,npart
        if(types(i)==1)ntype(1)=ntype(1)+1
        if(types(i)==2)ntype(2)=ntype(2)+1
        if(types(i)==3)ntype(3)=ntype(3)+1
        if(types(i)==4)ntype(4)=ntype(4)+1
     end do
!
     if(ntype(1)/=0)write(*,550)ctype(1),ntype(1)
     if(ntype(2)/=0)write(*,550)ctype(2),ntype(2)
     if(ntype(3)/=0)write(*,550)ctype(3),ntype(3)
     if(ntype(4)/=0)write(*,550)ctype(4),ntype(4)
!
!
550  FORMAT(2X,'Number of:',1X,A2,1X,I7)
!
!
     nbin=anint(gcut/bin)
!
!   write(*,*)'nbin=',nbin
!
     allocate(gr(nbin))
!
     do i=1,nbin
        gr(i)=0.0
     end do
!
!
     if(icom1==icom2)then
         uni=2
     else
         uni=1
     end if
!
!
!
       natm1=ntype(icom1)
       natm2=ntype(icom2)
!
!       write(*,*)natm1,natm2
!
!     write(*,*)'computed number of atom1',natm1
!     write(*,*)'computed number of atom2',natm2
!     write(*,*)'tatol mass of the nano-block',mass
!
!
!---- start of the g(r) function computation ----
!
        volume=0 !!!cube(1)*cube(2)*cube(3)
!
        rou=natm2/volume
!
!
!      n_pair=0
!
      do i=1,npart
!
         itype=types(i)
!
         if(itype/=icom1.and.itype/=icom2)cycle
!
         if(icom1==icom2)then
            sflag=i+1
         else
            sflag=1
         endif
!
         do j=sflag,npart
!
           jtype=types(j)
!
           if(jtype/=icom1.and.jtype/=icom2)cycle
!
           if(icom1/=icom2.and.jtype==itype)cycle
!
!
           rdx=x(1,i)-x(1,j)
       	   !rdx=rdx-cube(1)*anint(rdx/cube(1))
           rdy=x(2,i)-x(2,j)
	   !rdy=rdy-cube(2)*anint(rdy/cube(2))
           rdz=x(3,i)-x(3,j)
	   !rdz=rdz-cube(3)*anint(rdz/cube(3))
           rd=sqrt(rdx*rdx+rdy*rdy+rdz*rdz)
!
!
            if(rd<=(gcut-bin))then
!               n_pair=n_pair+uni
!               write(*,*)n_pair
               ind=anint(rd/bin)+1
               rsq=(ind-1)*bin
               rsq=rsq**2
               gr(ind)=gr(ind)+uni/(rou*4*pi*rsq*bin*natm1)
           end if
        end do
       end do
!
!
! ----the end of compute----
!
   write(50,*)"#g(r)"
!
!
     do i=1,nbin
!
         g=gr(i)
         tem=(i-1)*bin
         write(50,600)tem,g
!
     end do
!
!
!---- close memory ----
!
     deallocate(gr)
!
     close(50)
!
!
500  format(2X,A)
600  format(2F15.5)
!
     return
!
     end
!
