!                         WRITED BY GUO F.                 
!                                                                      
!               Institute of Atomic & Molecular Physics of Si-          
!                     chuan University of China.                       
!                                                                       
!                   email: gfeng.alan@hotmail.com                     
!
!
!-------------------- grid operations -------------------
! 
!
!-------------------------------------------------------------------------
!
!

subroutine memory
    use mode, only: npart,x,types,xd,v,q,alive,ctype
    use modnei,only:atable,matom,rcutsq,rcut,ang_value,dih_value,bond_value
    use modfim,only:cflag,molecule,molelist,mlist
    use lattice
    !use mod_def,only: idef,defect
    use mod_grd, only: fix
    implicit none

    integer      ic
    integer      i
    integer      j
    real(8)      rc(4)
!
!
!   data rc /1.3742, 0.6867, 1.3142, 1.2456/
    data rc /1.3742, 0.6867, 1.3142, 1.35/
! set memory

    !ctype = ' '
    !ctype(1) = 'C'
    !ctype(2) = 'H'
    !ctype(3) = 'O'
    !ctype(4) = 'N'

    ic = matom+1

    allocate(atable(ic,npart))
    allocate(x(3,npart))
    !print *,'allocate memory ...'
    allocate(q(npart))
    allocate(v(3,npart))
    allocate(alive(npart))
    allocate(xd(3,npart))
    allocate(types(npart))

    allocate(fix(npart))
!
    allocate(cflag(npart))
    allocate(molecule(npart))
    allocate(mlist(npart))
    allocate(molelist(npart+1))
!
!   lattice index

    allocate(Natom(npart))
    allocate(LatticeIndex(npart,3))

    allocate(bond_value(2*npart))
    allocate(ang_value(3*npart))
    allocate(dih_value(4*npart))

   alive(:) = 1

    do i=1,4
       do j=1,4
          rcut(i,j)=1.35*(0.5*rc(i)+0.5*rc(j))
          !if(ctype(i)=='O'.and.ctype(j)=='O') rcut(i,j)= 1.15 
          ! for O-O in HMX, not a bond, and for O2, abover line should be
          ! Commented 
        enddo
    enddo

    do i=1,4
       do j=1,4
          rcutsq(i,j)= rcut(i,j)*rcut(i,j)
       enddo
    enddo

    return
end
!
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
!
!

subroutine release
   use mode, only: x,types,q,v,xd,alive  !npart,
   use modnei,only:atable,matom,ang_value,dih_value
   use modfim,only:cflag,molecule,molelist,mlist
   use mod_def,only: Defect
   use lattice
   use mod_grd, only: fix
   use mod_nw , only: nw_bond, nw_angle, nw_dihe, nw_bond_val, nw_angle_val, atom_type,&
                       nw_dihe_val
   implicit none

! release memory

   if(allocated(atable))deallocate(atable)
   if(allocated(x))deallocate(x)

   if(allocated(q))deallocate(q)
   if(allocated(v))deallocate(v)
   if(allocated(xd))deallocate(xd)
   if(allocated(types))deallocate(types)
   if(allocated(alive)) deallocate(alive)
   if(allocated(fix)) deallocate(fix)

   if(allocated(cflag))deallocate(cflag)
   if(allocated(molecule))deallocate(molecule)
   if(allocated(mlist))deallocate(mlist)
   if(allocated(molelist)) deallocate(molelist)

   if(allocated(Natom)) deallocate(Natom)
   if(allocated(LatticeIndex)) deallocate(LatticeIndex)
   if(allocated(Defect)) deallocate(Defect)

   if(allocated(ang_value)) deallocate(ang_value)
   if(allocated(dih_value)) deallocate(dih_value)
!!
!! in nwchem module
   if(allocated(atom_type)) deallocate(atom_type)
   if(allocated(nw_bond)) deallocate(nw_bond)
   if(allocated(nw_angle)) deallocate(nw_angle)
   if(allocated(nw_dihe)) deallocate(nw_dihe)
   if(allocated(nw_bond_val)) deallocate(nw_bond_val)
   if(allocated(nw_angle_val)) deallocate(nw_angle_val)
   if(allocated(nw_dihe_val)) deallocate(nw_dihe_val)
   !if(allocated(ang_name)) deallocate(ang_name)
   return
end
!
!-------------------------------------------------------------------------
!
