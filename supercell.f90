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
!  ------             SUPER CELL        ------                           -
!-------------------------------------------------------------------------
!
!
    subroutine supercell
    use mode, only: npart,x,cpart,sc,types,a0,b0,c0,a,b,c !,cel
    use lattice
    implicit none

    integer i,j,k,l,ii

!   write(6,100)'CREATE SUPE-CELL ...'
!   write(6,*)sc(1), sc(2), sc(3)

    ii=sc(1)*sc(2)*sc(3)*cpart
    if(ii.NE.npart)call err('ERROR: CELL ATOM NUMBER ERROR!')

    do i=1,3
        a(i)= a0(i)*sc(i)
        b(i)= b0(i)*sc(i)
        c(i)= c0(i)*sc(i)
    end do


    do i=1, sc(1)
        do j=1, sc(2)
           do k=1, sc(3)
              do l=1, cpart
                 ii= l+(k-1)*cpart + (j-1)*sc(3)*cpart + (i-1)*sc(2)*sc(3)*cpart
                 x(1,ii)=x(1,l) + (k-1)*c0(1) +(j-1)*b0(1) +(i-1)*a0(1)
                 x(2,ii)=x(2,l) + (k-1)*c0(2) +(j-1)*b0(2) +(i-1)*a0(2)
                 x(3,ii)=x(3,l) + (k-1)*c0(3) +(j-1)*b0(3) +(i-1)*a0(3)
                 types(ii)=types(l)     

                 ! lattice indices
                 ! used by lammps FixPhonon commond
                 Natom(ii) = l-1
                 LatticeIndex(ii,1) = i-1
                 LatticeIndex(ii,2) = j-1
                 LatticeIndex(ii,3) = k-1          
              enddo
           enddo
       enddo
    enddo
    !print *,'in super-cell function ...'
    !print *,x(:,384)
   

    return
    end
!
!-------------------------------------------------------------------------
!
