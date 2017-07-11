!=========================================================================
!                              WRITED BY GUO F.                              
!                                                                         
!              Institute of Atomic & Molecular Physics of Liao-             
!                             cheng University of China.                        
!                                                                        
!                           email: gfeng.alan@gmail.com                      
!=========================================================================

!----------------- analysis the dipole moments of molecules ------------------

subroutine dipole
    use mode,     only: x, npart,q
    use modfim,   only: n_mol,mlist,molelist
    implicit none

    integer     i, j, jj, iatom
    real(8)      r1,xpc,ypc,zpc,xnc,ync,znc,pc,nc,d
    character(6) c1

    !allocate(q(npart))

    open(222,file='charge.txt',status='old')
    do i=1,npart
          read(222,*) c1,c1,c1,r1,q(i)
    enddo

    open(223,file='dipole.xyz',status='unknown')
    write(223,*) n_mol*2
    write(223,*) 

    do i=1,n_mol  ! ! molecular index 

       if(i>1)then
            jj=mlist(i-1)+1
       else
            jj=1
       endif

       xpc = 0.0
       ypc = 0.0
       zpc = 0.0
       xnc = 0.0
       ync = 0.0
       znc = 0.0
       pc = 0.0
       nc = 0.0 
       do j=jj, mlist(i)
            iatom = molelist(j)

            if (q(iatom)>0) then
                   xpc = xpc + q(iatom)*x(1,iatom)
                   ypc = ypc + q(iatom)*x(2,iatom)
                   zpc = zpc + q(iatom)*x(3,iatom)
                   pc = pc + q(iatom)
            elseif (q(iatom)<0) then
                   xnc = xnc + q(iatom)*x(1,iatom)
                   ync = ync + q(iatom)*x(2,iatom)
                   znc = znc + q(iatom)*x(3,iatom)
                   nc = nc + q(iatom)
            endif
       enddo

       print *,'Mole', i, 'negatively charged:', nc
       print *,'Mole', i, 'positively  charged:', pc
       xpc = xpc /pc
       ypc = ypc/pc
       zpc = zpc/pc
       xnc = xnc /nc
       ync = ync/nc
       znc = znc/nc
       write(223,*)'e+',xpc,ypc,zpc
       write(223,*)'e-',xnc,ync,znc
       d = sqrt((xpc-xnc)**2 + (ypc-ync)**2 + (zpc-znc)**2 )
       print *,'d of the dipole:',d
       print *, 'dipole of the molecule:',d*pc
    enddo

    close(222)
    close(223)

    return
end

!----------------- analysis the dipole moments of molecules ------------------

