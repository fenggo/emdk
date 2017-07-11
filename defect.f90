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
!   Subroutine to create circular defects
!

subroutine defects
    use mode,     only:  types,x, npart !,a,b,c  x,
    use mod_def,  only: idef,Defect
    use modfim,   only: n_mol,mlist,molelist
    implicit none

    integer            i,j,l,jj!,ll
    integer            atom
    integer            nd,ndd
    integer            isur   ! if molecules on surface
    integer            def(npart) ,dtyp(npart) ! atom list to delete
    !integer            ina,inb,inc

    real(8)            delx,dely,delz,rd,r2
    !real(8)            r_a,r_b,r_c
    real(8)            xdef(3,npart)

    call neighbor 
    call findmole

    nd = 0
    xdef = x
    dtyp = types

    do j=1,n_mol        ! index from molecues
       if(j>1)then
           jj=mlist(j-1) + 1
       else
           jj=1
       endif

       isur = 0

       do l=jj, mlist(j)      
          do i=1,idef 
             atom = molelist(l)
             delx=x(1,atom)-Defect(1,i)
             dely=x(2,atom)-Defect(2,i)
             delz=x(3,atom)-Defect(3,i)

!   apply PBC

             !r_c = (a(1)*b(2)-a(2)*b(1))*(b(2)*delz-b(3)*dely) - (a(3)*b(2)-a(2)*b(3))*(b(2)*delx-b(1)*dely)
             !r_c = r_c/((a(1)*b(2)-a(2)*b(1))*(c(3)*b(2)-c(2)*b(3)) - (a(3)*b(2)-a(2)*b(3))*(c(1)*b(2)-c(2)*b(1)))

             !r_a = b(2)*delx-b(1)*dely - (c(1)*b(2)-c(2)*b(1))*r_c
             !r_a = r_a/(a(1)*b(2)-a(2)*b(1))
             !r_b = dely- r_a*a(2) - r_c*c(2)
             !r_b = r_b/b(2)

             !ina = anint(r_a)
             !inb = anint(r_b)
             !inc = anint(r_c)

             !delx= delx - ina*a(1) - inb*b(1) - inc*c(1)
             !dely= dely - ina*a(2) - inb*b(2) - inc*c(2)
             !delz= delz - ina*a(3) - inb*b(3) - inc*c(3)

!  R^2 after PBC
             r2 = delx*delx+dely*dely+delz*delz
             rd = sqrt(r2)
             if (rd < defect(4,i)) isur=1
          enddo
       enddo
       if (isur == 1) then
           do l=jj, mlist(j)  
              nd = nd + 1
              def(nd) = molelist(l)
           enddo
       endif
    enddo

! end of check atoms

    ndd  = 0

    do i=1,npart
       isur = 0
       do j=1,nd
          if(i==def(j))then
             isur = 1
             exit
          endif
       enddo
       if (isur==0)then
          ndd = ndd+1
          x(:,ndd)=xdef(:,i)
          types(ndd)=dtyp(i)
       endif
    enddo

    npart = npart - nd
    if (npart/=ndd) then
         print *,'Error found!'
         print *,npart,nd,ndd
    endif

    return
end


!------------- remove water or other molecules -----------------

subroutine remove
    use mode,     only: types,ctype, alive !, npart x, 
    !use mod_def,  only: idef
    use modfim,   only: n_mol,mlist,molelist
    implicit none

    integer     i, j, jj, iatom
    integer     n_c, n_h, n_o, n_n, n_h2o

    !call neighbor 
    !call findmole
    n_h2o = 0
    do i=1,n_mol  ! ! molecular index 

       if(i>1)then
            jj=mlist(i-1)+1
       else
            jj=1
       endif
       n_c = 0
       n_h = 0
       n_n = 0
       n_o =0
       do j=jj, mlist(i)
          iatom = molelist(j)
           if(ctype(types(iatom))=='C')n_c=n_c+1
           if(ctype(types(iatom))=='H')n_h=n_h+1
           if(ctype(types(iatom))=='O')n_o=n_o+1
           if(ctype(types(iatom))=='N')n_n=n_n+1
       enddo

       if(n_c==0.and.n_n==0.and.n_h==2.and.n_o==1) then  ! remove 2/3 H2O molecules
          n_h2o = n_h2o + 1
          !if (mod(n_h2o,3)==0)then
          do j=jj, mlist(i)
             iatom = molelist(j)
             alive(iatom) = 0
          enddo
          !endif
       endif

    enddo

    return
end

!------------- remove water or other molecules -----------------

