!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!
! build the neighbor-list
!
!
    subroutine neighbor
    use mode,only: x, types, npart,a,b,c  !, cstep
    use modnei, only: mbond, matom, rcutsq, nbond, atable, bond
    implicit none
!
    integer                 i,j, itype,jtype, cb, ib
    integer                 ina,inb,inc

    real(8)                 delx,dely,delz,rsq
    real(8)                 r_a,r_b,r_c
!
!
!-------------------------------------------------------------------------
!
!    build the density-list(neighbor-list)
!
!    listcou=0
!
     do i=1,npart
        do j=1, matom
          atable(j,i)=0
        end do
     end do
!
     cb = 0
!
     do i=1,npart-1
!
!        ilist=0 ! counter
!
        do j=i+1,npart
!
!    Coordinate vector 
!            R= (R1 - R2) = delx.i + dely.j + delz.k  
!            R.a/abs(a)= Rcos(alpha)

            delx=x(1,i)-x(1,j)
            dely=x(2,i)-x(2,j)
            delz=x(3,i)-x(3,j)

!   apply PBC
        
            r_c = (a(1)*b(2)-a(2)*b(1))*(b(2)*delz-b(3)*dely) - (a(3)*b(2)-a(2)*b(3))*(b(2)*delx-b(1)*dely)
            r_c = r_c/((a(1)*b(2)-a(2)*b(1))*(c(3)*b(2)-c(2)*b(3)) - (a(3)*b(2)-a(2)*b(3))*(c(1)*b(2)-c(2)*b(1)))

            r_a = b(2)*delx-b(1)*dely - (c(1)*b(2)-c(2)*b(1))*r_c
            r_a = r_a/(a(1)*b(2)-a(2)*b(1))
            r_b = dely- r_a*a(2) - r_c*c(2)
            r_b = r_b/b(2)

            ina = anint(r_a)
            inb = anint(r_b)
            inc = anint(r_c)

            delx= delx - ina*a(1) - inb*b(1) - inc*c(1)
            dely= dely - ina*a(2) - inb*b(2) - inc*c(2)
            delz= delz - ina*a(3) - inb*b(3) - inc*c(3)

!             delx=x(1,i)-x(1,j)
!             dely=x(2,i)-x(2,j)
!             delz=x(3,i)-x(3,j)
!  R^2 after PBC
            rsq=delx*delx+dely*dely+delz*delz
!
            itype=types(i)
            jtype=types(j)
!
           if(rsq<=rcutsq(itype,jtype))then
!
               cb=cb+1
!
               if(cb>mbond)call err('ERROR: Bondlist over float')
!
!   --------  construction of local bond list  --------
!
                bond(1,cb)=i
                bond(2,cb)=j
!
!   --------  construction of local bond list  --------
!
               atable(1,i)=atable(1,i)+1
               atable(1,j)=atable(1,j)+1
!
!
               if(atable(1,i)>matom) call err('Too many atoms bonded to one, Modify the maxatom and retry!')

               if(atable(1,j)>matom) call err('Too many atoms bonded to one, Modify the maxatom and retry!')
!
!  -------- count the bond of atom i and j --------
!
               ib=atable(1,i)+1
               atable(ib,i)=j
               ib=atable(1,j)+1
               atable(ib,j)=i
!
!  -------- count the bond of atom i and j --------
!
           end if
!
        end do
!
!  ----atable = find up list for neighbor-list-----
!
     end do
!
     nbond=cb  ! total number of bond
!
    return
    end
!
!-------------------------------------------------------------------------
!
