!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!

subroutine mmfit
   use modnei,only: nbond,bond,n_ang,n_dih,angle,dihedral,bond_value, &
                    ang_value,dih_value
   !use mode,  only: x !, a, b, c ,alive ,ctype ,npart ,q ,types
   use mod_nw,only: atom_type
   use mod_mmfit,only: task
   implicit none

   integer          i
   integer          iatom,jatom,katom,latom


   call  neighbor
   call  findmole

   call  bonds
   call  angles
   call  torsions

   !if(.not. allocated(atom_type)) then
      !allocate(atom_type(npart))
      !atom_type  = 'null '

      !do i=1,npart
         !atom_type(i) = ctype(types(i))
         !print *,atom_type(i)
      !enddo
   !endif
   call atm_typ

   if (task == 1)then   ! bat: bond angle torsion 
      open(35,file='bonds.txt',status='unknown')
      open(36,file='angles.txt',status='unknown')
      open(37,file='tors.txt',status='unknown')
      do i =1,nbond
         iatom = bond(1,i)
         jatom = bond(2,i)
         write(35,'(2A7,F11.6,2(1X,I5))')atom_type(bond(1,i)),atom_type(bond(2,i)),& 
                    & bond_value(i),iatom,jatom
      enddo
      do i =1,n_ang
         iatom = angle(1,i)
         jatom = angle(2,i)
         katom = angle(3,i)
         write(36,'(3A7,F11.6,3(1X,I5))')atom_type(iatom),atom_type(jatom),&
             atom_type(katom),ang_value(i),iatom,jatom,katom
      enddo
      do i =1,n_dih
         iatom = dihedral(1,i)
         jatom = dihedral(2,i)
         katom = dihedral(3,i)
         latom = dihedral(4,i)
         write(37,'(4A7,F11.6,4(1X,I5))')atom_type(iatom),atom_type(jatom),&
             atom_type(katom),atom_type(latom),dih_value(i),iatom,jatom,katom,latom
      enddo
   elseif(task == 2) then ! check bonds angles and torsions varied or not
      open(35,file='bonds.txt',status='old')
      open(36,file='angles.txt',status='old')
      open(37,file='tors.txt',status='old')

      do i=1,nbond
         read(35,*)atom_type(bond(1,i)),atom_type(bond(2,i)),bond_value(i)
      enddo
      do i =1,n_ang
         read(36,*)atom_type(angle(1,i)),atom_type(angle(2,i)),&
             atom_type(angle(3,i)),ang_value(i)
      enddo
      do i =1,n_dih
         read(37,*)atom_type(dihedral(1,i)),atom_type(dihedral(2,i)),&
             atom_type(dihedral(3,i)),atom_type(dihedral(4,i)),dih_value(i)
      enddo
      call w_bond
      call w_angle
      call w_torsion
   endif
   close(35)
   close(36)
   close(37)
   return
end
