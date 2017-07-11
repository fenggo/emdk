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
!     Read In LAMMPS Trajectory                                          -
!-------------------------------------------------------------------------

     subroutine readltr
     use mode, only: cpart,x,types,q,v,a0,b0,c0,r1,r2,r3, &
                     alpha,beta,gamm,npart
     use modctrl,only:  cfile
     implicit none

     integer        i ,id
     integer        j
     !integer        typ
     real(8)               lo,hi
     real(8)               xlo_bound, xhi_bound, xy, xlo, xhi
     real(8)               ylo_bound, yhi_bound, xz, ylo, yhi
     real(8)               zlo_bound, zhi_bound, yz, zlo, zhi
     real(8)               cos_alpha,cos_beta,cos_gamma
     !character(50)  ltr
     character(80)  line
     character             head

     !print 100,'Input LAMMPS trajectory file name:'
     !read(*,*)ltr
     open(11,file=cfile,status='old')

     read(11,*)head
!
     read(11,*)!step
     read(11,*)
     read(11,*)npart
     read(11,'(A80)')line
!
!
      if(index(line,'ITEM: BOX BOUNDS')/=0)then
           if(index(line,'xy xz yz')/=0)then      
                read(11,*)xlo_bound, xhi_bound, xy
                read(11,*)ylo_bound, yhi_bound, xz
                read(11,*)zlo_bound, zhi_bound, yz
                xlo= xlo_bound - min(0.0,xy,xz,xy+xz)
                xhi= xhi_bound - MAX(0.0,xy,xz,xy+xz)
                ylo= ylo_bound - MIN(0.0,yz)
                yhi= yhi_bound - MAX(0.0,yz)
                zlo= zlo_bound 
                zhi= zhi_bound
                a0 = (/xhi-xlo, 0.0d0, 0.0d0/)
                b0 = (/xy, yhi-ylo, 0.0d0/)
                c0 = (/xz, yz, zhi-zlo/)

                r1 = a0(1)*a0(1)+ a0(2)*a0(2)+ a0(3)*a0(3)
                r1 = dsqrt(r1)
                r2 = b0(1)*b0(1)+ b0(2)*b0(2)+ b0(3)*b0(3)
                r2 = dsqrt(r2)
                r3 = c0(1)*c0(1)+ c0(2)*c0(2)+ c0(3)*c0(3)
                r3 = dsqrt(r3)

                cos_alpha = (b0(1)*c0(1) + b0(2)*c0(2) + b0(3)*c0(3))/(r2*r3)
                cos_beta  = (c0(1)*a0(1) + c0(2)*a0(2) + a0(3)*c0(3))/(r3*r1)
                cos_gamma = (a0(1)*b0(1) + a0(2)*b0(2) + a0(3)*b0(3))/(r1*r2)

                alpha = dacos(cos_alpha)
                beta  = dacos(cos_beta)
                gamm = dacos(cos_gamma)
           else
               read(11,*)lo,hi
               a0(1)= hi-lo
               r1= a0(1)
               read(11,*)lo,hi
               b0(2)= hi-lo
               r2= b0(2)
               read(11,*)lo,hi
               c0(3)= hi-lo
               r3= c0(3)

               alpha = dacos(0.0d0)
               beta  = dacos(0.0d0)
               gamm = dacos(0.0d0)

           endif
     endif


     read(11,*)
     do i=1,cpart
        read(11,*)id,types(id),(x(j,id),j=1,3),q(id),(v(j,id),j=1,3)
     enddo

     close(11)
     return
     end

!-------------------------------------------------------------------------
!     Read In    XYZ  Trajectory                                         -
!-------------------------------------------------------------------------

     subroutine readxyz
     use mode, only: cpart,x,types,ctype
     use modctrl,only:  cfile
     implicit none

     integer        i!,id
     integer        j
     character(4)   cha

!     print 100,'Input XYZ trajectory file name:'
!     read(*,*)cfile
     open(11,file=cfile,status='old')

     read(11,*)
     read(11,*,err=99)!cel(1),cel(2),cel(3)

     goto 999
99   print 100,'Lattice parameter read in input file!'
999  continue

     do i=1,cpart
        read(11,*)cha,(x(j,i),j=1,3)
        !print *,cha
        !print *,ctype
        if(cha==ctype(1))types(i)=1
        if(cha==ctype(2))types(i)=2
        if(cha==ctype(3))types(i)=3
        if(cha==ctype(4))types(i)=4
!        print *,types(i)
     enddo

100  format(2X,A)

     close(11)
     return
     end

!-------------------------------------------------------------------------
!    Read In CIF file  .cif                                              -
!-------------------------------------------------------------------------

subroutine read_cif
     use mode, only: xd,types,ctype,a0,b0,c0,nelem,a,b,c
     implicit none

     integer        i,j,k,l
     integer        eln(20)
     character(80)  line
     real(8)        ac0
     !real(8),allocatable:: xf(:,:)

     open(11,file='POSCAR',status='old')
     read(11,'(A)')line
     print 100,line
     read(11,*)ac0
     read(11,*)a0
     a0(:)= a0(:)*ac0
     read(11,*)b0
     b0(:)= b0(:)*ac0
     read(11,*)c0
     c0(:)= c0(:)*ac0
     read(11,*)(ctype(i),i=1,nelem)
     read(11,*)(eln(i),i=1,nelem)


     read(11,100)line
     k = 0
     do l = 1,nelem
        do i=1,eln(l)
           k = k+1
           read(11,*)(xd(j,k),j=1,3)
           types(k)=l
        enddo
      enddo

    do i=1,3
        a(i)= a0(i)
        b(i)= b0(i)
        c(i)= c0(i)
    enddo

    call cart_coord

100  format(2x,A)
     close(11)
     return
end

!-------------------------------------------------------------------------
!    Read In VASP POSCAR                                                 -
!-------------------------------------------------------------------------

subroutine readposcar
     use mode, only: xd,types,ctype,a0,b0,c0,nelem,a,b,c
     implicit none

     integer        i,j,k,l
     integer        eln(20)
     character(80)  line
     real(8)        ac0
     !real(8),allocatable:: xf(:,:)

     open(11,file='POSCAR',status='old')
     read(11,'(A)')line
     print 100,line
     read(11,*)ac0
     read(11,*)a0
     a0(:)= a0(:)*ac0
     read(11,*)b0
     b0(:)= b0(:)*ac0
     read(11,*)c0
     c0(:)= c0(:)*ac0
     read(11,*)(ctype(i),i=1,nelem)
     read(11,*)(eln(i),i=1,nelem)


     read(11,100)line
     k = 0
     do l = 1,nelem
        do i=1,eln(l)
           k = k+1
           read(11,*)(xd(j,k),j=1,3)
           types(k)=l
        enddo
      enddo

    do i=1,3
        a(i)= a0(i)
        b(i)= b0(i)
        c(i)= c0(i)
    enddo

    call cart_coord

100  format(2x,A)
     close(11)
     return
end
!
!-------------------------------------------------------------------------
!        Read in "CONTCAR"  format coordinate generated by VASP          -
!-------------------------------------------------------------------------
!

subroutine read_contcar
     use mode,    only: x,cpart,types,a0,b0,c0,nelem
     use modctrl, only: cfile
     implicit none

     integer                   i,j,k,l,m
     integer,allocatable::     cpar(:)
     character(80)             line
     real(8)                   ac0
     real(8),allocatable::     xf(:,:)


     allocate(cpar(nelem))

     open(11,file=cfile,status='old')
     read(11,'(A)')line
     print 100,line
     read(11,*)ac0
     read(11,*)a0
     a0(:)= a0(:)*ac0
     read(11,*)b0
     b0(:)= b0(:)*ac0
     read(11,*)c0
     c0(:)= c0(:)*ac0
     read(11,'(80A)')line
     print 100,'element:'
     print 100,line
     read(line,*)(cpar(i),i=1,nelem)

!     print *,cpar(1),cpar(2)

     allocate(xf(3,cpart))

     l= 0
     read(11,100)line
     do j=1,nelem
        do i=1,cpar(j)
           read(11,*)(xf(m,i),m=1,3)

           l= l+1
           types(l)=j
! for Non-orthogonal(triclinic symmetry) simulation boxes

           do k=1,3
              x(k,l)= xf(1,i)*a0(k) + xf(2,i)*b0(k) + xf(3,i)*c0(k) ! x = R*(a,b,c)^T
           enddo
        enddo
     enddo

!     print *,xf(1,21),xf(2,21),xf(3,21)

100  format(2x,A)
     close(11)
     return
end

!-------------------------------------------------------------------------
!     Read in ".car"  format coordinate generated by Material Studio     -
!-------------------------------------------------------------------------

     subroutine readcar
     use mode, only: cpart,x,types,ctype
     use modctrl, only: cfile
     implicit none

     integer        i,id
     integer        j
     !integer        typ
     real           q
!     character(50)  xyz
     character(60)  line
     character(4)   cha,tcha

!     print 100,'Input .car trajectory file name:'
!     read(*,*)xyz
     open(11,file=cfile,status='old')

     read(11,*)
     read(11,*)line
     if((index(line,'PBC').NE.0).and.(index(line,'ON').NE.0))then
     read(11,*)
     endif
     read(11,*)!cel(1),cel(2),cel(3)
     read(11,*)

     do i=1,cpart
        read(11,*)tcha,(x(j,i),j=1,3),tcha,id,tcha,cha,q
        if(cha==ctype(1))types(i)=1
        if(cha==ctype(2))types(i)=2
        if(cha==ctype(3))types(i)=3
        if(cha==ctype(4))types(i)=4
     enddo

     close(11)
     return
     end
!
!-------------------------------------------------------------------------
!                                                                        -
!-------------------------------------------------------------------------
!
    subroutine read_lmpdata
    use mode,    only: cpart,x,v,types,npart,mass
    use modctrl, only: cfile
    implicit none

    integer            ns
    integer            i,j,l
    real(8)            q
    character(80)      line

    open(11,file=cfile,status='unknown')

    do
      read(11,'(A80)',end=990),line
      if    ((index(line,'atoms').NE.0))then
            read(line,*)cpart
            print '(A8,I8,1x,A5)','Read in',cpart,'atoms'
!            allocate(types(npart)) 
!            allocate(x(3,npart))
!            allocate(v(3,npart))
      elseif((index(line,'atom types').NE.0))then 
            read(line,*)ns
            print '(I3,A10)',ns,' atom types'
      elseif((index(line,'Masses').NE.0))then 
            print '(A16)','read in mass ...'
            read(11,*)
            do l=1,ns
               read(11,*)j,mass(j)
            enddo
            print '(A16)','read in mass ...'
      elseif((index(line,'Atoms').NE.0))then 
            print '(A27)','read in atom coordinate ...'
            read(11,*)
            do i=1,npart
               read(11,*)j,types(j),q,(x(l,j),l=1,3)
            enddo
      elseif((index(line,'Velocities').NE.0))then 
            print '(A22)','read in velocities ...'
            read(11,*)
            do i=1,cpart
               read(11,*)j,(v(l,j),l=1,3)
            enddo
      endif
    enddo
990 continue

    close(11)

    return
    end
!
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

     subroutine read_pdb
     use mode, only: cpart,x,types,ctype,nelem
     use modctrl,only:  cfile
     implicit none

     integer        i!,id
     integer        j,ee
     !integer        typ
!     character(50)  xyz
     character(4)   cha
     character(100) line

!     print 100,'Input XYZ trajectory file name:'
!     read(*,*)cfile
     open(11,file=cfile,status='old')

!     read(11,*)
!     read(11,*,err=99)!cel(1),cel(2),cel(3)

     goto 999
     print 100,'Lattice parameter read in input file!'
999  continue

     i= 0
     do 
        read(11,'(A100)',end=990)line
        if(index(line,'ATOM')/=0.or.index(line,'HETATM')/=0)then
              i= i + 1
              read(line,'(30X,3F8.3)')(x(j,i),j=1,3)
              read(line,'(76X,A2)')cha
              cha=adjustl(cha)
              do ee=1, nelem
                if(index(cha,ctype(ee)).NE.0)types(i)=ee
              enddo
        endif

        !print *,ctype(types(i))
        !print *,cha

     enddo
990  continue

     if(i>cpart)print '(A29)','Atom number > defined in .pdb'
     if(i<cpart)print '(A29)','Atom number < defined in .pdb'

100  format(2X,A)

     close(11)
     return
     end
!-------------------------------------------------------------------------
!                   read in gulp input                                   -
!-------------------------------------------------------------------------

subroutine read_gulp
     use mode, only: cpart,xd,types,ctype,r1,r2,r3,alpha,beta,&
                     gamm,a0,b0,c0,a,b,c
     use modctrl,only:  cfile
     implicit none

     integer        i!,id
     integer        j!,ee
     !integer        typ

     real(8)        cy!,ay
     real(8) ::     pi = 3.1415927

     character(4)   cha
     character(100) line

     pi = dacos(-1.0d0)

     open(11,file=cfile,status='old')

     !print *,cfile
     i= 0
     do 
        read(11,'(A100)',end=990)line
        !print *,line
        if(index(line,'cell')/=0)then
            read(11,*,err=990)r1,r2,r3,alpha,beta,gamm
            alpha= alpha*pi/180.0d0
            beta = beta*pi/180.0d0
            gamm = gamm*pi/180.0d0

            a0 = (/r1,0.0d0,0.0d0/)
            b0 = (/r2*dcos(gamm),r2*dsin(gamm),0.0d0/)
            cy= r3*dcos(alpha)-r3*dcos(beta)*dcos(gamm)
            cy= cy/dsin(gamm)
            c0 = (/r3*dcos(beta),cy,dsqrt(r3*r3*dsin(beta)*dsin(beta)-cy*cy)/)
        endif

        if(index(line,'fractional')/=0 .or. index(line,'cartesian')/=0)then
           !print *,line
           do while(index(line,'connect')==0)
              read(11,'(A100)',end=990)line
              !print *,line
              if(index(line,'connect')==0) then
                 i = i + 1
                 read(line,*)cha,(xd(j,i),j=1,3)
                 if (cha==ctype(1)) types(i) = 1
                 if (cha==ctype(2)) types(i) = 2
                 if (cha==ctype(3)) types(i) = 3
                 if (cha==ctype(4)) types(i) = 4
              endif
           enddo
        endif

     enddo
990  continue

     if(i>cpart)print '(A34)','Atom number > defined in gulpinput'
     if(i<cpart)print '(A29)','Atom number < defined in gulpinput'

     !print *,'partical number have found:',i

     do i=1,3
        a(i)= a0(i)
        b(i)= b0(i)
        c(i)= c0(i)
     enddo

     call cart_coord

     close(11)
     return
end
!
!----------------------------------华丽的分割线----------------------------------
!    read in DFTB+ gen format
! 

subroutine read_gen
     use mode, only: cpart,x,types,ctype,nelem,a0,b0,c0
     use modctrl,only:  cfile
     implicit none

     integer        i,id
     integer        j
     !  integer        typ
     !  character(50)  xyz
     character(4)   cha
     !character(100) line


     open(11,file=cfile,status='old')

     read(11,*)cpart,cha
     read(11,*)(ctype(j),j=1,nelem)

     do i=1,cpart
        read(11,*)id,types(i),(x(j,i),j=1,3)
     enddo


     if(cha=='S') then
          read(11,*)
          read(11,*) a0(1),a0(2),a0(3)
          read(11,*) b0(1),b0(2),b0(3)
          read(11,*) c0(1),c0(2),c0(3)
     endif

     close(11)
     return
end
!
!----------------------------------华丽的分割线----------------------------------
!
