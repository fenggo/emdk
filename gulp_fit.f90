!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                     =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!
!-------------------------------------------------------------------------
!        Use gulp fit the Molecular Mechanic Potententials (bondes)      -
!-------------------------------------------------------------------------
!

subroutine fit_bondes
    !use modfim,only: n_mol
    use modnei,only:  mang,mdih
    !use mode,  only: x,types
    use mod_nw,only: nw_bond, nw_bond_val,  &
                     atom_type,nw_bonds
    implicit none

    integer                      i,j,k
    integer                      bond_typs
    character(5),allocatable::   bond_name(:,:)
    character(8)             ::  bond_file = '        ' 
    character(4)             ::  con = '    '
    logical                      new

    allocate(bond_name(2,nw_bonds))

    bond_typs  = 0

    do i=1, nw_bonds
       if(i==1)then
           bond_typs = 1
           bond_name(1,bond_typs) = atom_type(nw_bond(1,i))
           bond_name(2,bond_typs) = atom_type(nw_bond(2,i))
       else 
           !print *,bond_name(1,ang_typs),atom_type(nw_angle(1,i))
           new = .true.
           do k=1,bond_typs
              if(bond_name(1,k)==atom_type(nw_bond(1,i)) .and. &
                 & bond_name(2,k)==atom_type(nw_bond(2,i)) ) then
                 new = .false.
              elseif(bond_name(2,k)==atom_type(nw_bond(1,i)) .and. &
                 & bond_name(1,k)==atom_type(nw_bond(2,i)) ) then
                 new = .false.
                 !print *,bond_name(1,ang_typs),bond_name(2,ang_typs),bond_name(3,ang_typs)
              endif
           enddo
           if(new) then
              bond_typs = bond_typs + 1
              bond_name(1,bond_typs) = atom_type(nw_bond(1,i))
              bond_name(2,bond_typs) = atom_type(nw_bond(2,i))
           endif
       endif
    enddo
    print '(A23,I5)', 'Number of bond types: ', bond_typs



   do i=1,bond_typs
       bond_file(1:5) = 'bond_'
       write(con,'(I4)') i
       con = adjustl(con)
       if (i<10) bond_file(6:6) = con
       if (i>=10 .and. i<100) bond_file(6:7) = con
       bond_file = trim(adjustl(bond_file))
       open(101,file=bond_file,status='unknown')
       write(101,'(A18)') 'Key words in gulp:'
       write(101,'(A16)') 'morse bond intra'
       write(101,'(2(A5,A6),A3,A3,A5)') bond_name(1,i),' core ',bond_name(2,i),' core ', & 
                                    & 'De ',' a ',' r_0 '
       do j=1,nw_bonds
           if(bond_name(1,i)==atom_type(nw_bond(1,j)) .and. &
              & bond_name(2,i)==atom_type(nw_bond(2,j)) ) then
              write(101,'(A5,2(I5,1X,A5),F10.5,I4)') 'Bond ',nw_bond(1,j),atom_type(nw_bond(1,j)), &
                                  nw_bond(2,j), atom_type(nw_bond(2,j)), &
                                  nw_bond_val(j),j
           elseif(bond_name(2,i)==atom_type(nw_bond(1,j)) .and. &
              & bond_name(1,i)==atom_type(nw_bond(2,j))  ) then
              write(101,'(A5,2(I5,1X,A5),F10.5,I4)') 'Bond ',nw_bond(1,j),atom_type(nw_bond(1,j)), &
                                  nw_bond(2,j), atom_type(nw_bond(2,j)), &
                                  nw_bond_val(j),j
           endif
       enddo
       close(101)
    enddo
    !print '(A)','**'
    !print '(A)','** Using keyword geometry  autoz in NWchem'
    !print '(A)','** geometry adjust'
    !print '(A)','**    zcoord'
    !print '(A)','**       bond'
    !print '(A)','**    end'
    !print '(A)','** end'
    !print '(A)','**'

    if(allocated(bond_name)) deallocate(bond_name)
    return
end

!
!-------------------------------------------------------------------------
!        Use gulp fit the Molecular Mechanic Potententials (angles)      -
!-------------------------------------------------------------------------
!

subroutine fit_angles
    !use modfim,only: n_mol
    use modnei,only: mang,mdih
    !use mode,  only: x, a, b, c,npart,types,ctype
    use mod_nw,only:  nw_angle, nw_angle_val, &
                      atom_type, nw_bonds,nw_angles
    implicit none

    integer                      i,j,k
    integer                      ang_typs
    character(5),allocatable::  ang_name(:,:)
    character(8)             ::  ang_file = '        ' 
    character(4)             ::  con = '    '
    logical                      new

    allocate(ang_name(3,nw_angles))

    ang_typs  = 0

    do i=1, nw_angles
       if(i==1)then
           ang_typs = 1
           ang_name(1,ang_typs) = atom_type(nw_angle(1,i))
           ang_name(2,ang_typs) = atom_type(nw_angle(2,i))
           ang_name(3,ang_typs) = atom_type(nw_angle(3,i))
       else 
           !print *,ang_name(1,ang_typs),atom_type(nw_angle(1,i))
           new = .true.
           do k=1,ang_typs
              if(ang_name(1,k)==atom_type(nw_angle(1,i)) .and. &
                 & ang_name(2,k)==atom_type(nw_angle(2,i)) .and. &
                 & ang_name(3,k)==atom_type(nw_angle(3,i)) ) then
                 new = .false.
              elseif(ang_name(3,k)==atom_type(nw_angle(1,i)) .and. &
                 & ang_name(2,k)==atom_type(nw_angle(2,i)) .and. &
                 & ang_name(1,k)==atom_type(nw_angle(3,i)) ) then
                 new = .false.
                 !print *,ang_name(1,ang_typs),ang_name(2,ang_typs),ang_name(3,ang_typs)
              endif
           enddo
           if(new) then
              ang_typs = ang_typs + 1
              ang_name(1,ang_typs) = atom_type(nw_angle(1,i))
              ang_name(2,ang_typs) = atom_type(nw_angle(2,i))
              ang_name(3,ang_typs) = atom_type(nw_angle(3,i))
           endif
       endif
    enddo
    print '(A23,I5)', 'Number of angle types: ', ang_typs


    do i=1,ang_typs
       !print *,ang_name(1,i),ang_name(2,i),ang_name(3,i)
       ang_file(1:5) = 'angl_'
       write(con,'(I4)') i
       con = adjustl(con)
       if (i<10) ang_file(6:6) = con
       if (i>=10 .and. i<100) ang_file(6:7) = con
       ang_file = trim(adjustl(ang_file))
       open(101,file=ang_file,status='unknown')
       write(101,'(A18)') 'Key words in gulp:'
       write(101,'(A16)') 'three bond intra'
       write(101,'(3(A5,A6),A2,A3)') ang_name(2,i),' core ',ang_name(1,i),' core ',ang_name(3,i),& 
                     & ' core ', 'k ','r_0'
       do j=1,nw_angles
           if(ang_name(1,i)==atom_type(nw_angle(1,j)) .and. &
              & ang_name(2,i)==atom_type(nw_angle(2,j)) .and. &
              & ang_name(3,i)==atom_type(nw_angle(3,j)) ) then
              write(101,'(A4,3(I5,1X,A5),F10.5,I4)')'Ang ', nw_angle(1,j),atom_type(nw_angle(1,j)), &
                                  nw_angle(2,j), atom_type(nw_angle(2,j)), &
                                  nw_angle(3,j), atom_type(nw_angle(3,j)),nw_angle_val(j),j+nw_bonds
           elseif(ang_name(3,i)==atom_type(nw_angle(1,j)) .and. &
              & ang_name(2,i)==atom_type(nw_angle(2,j)) .and. &
              & ang_name(1,i)==atom_type(nw_angle(3,j)) ) then
              write(101,'(A4,3(I5,1X,A5),F10.5,I4)')'Ang ', nw_angle(1,j),atom_type(nw_angle(1,j)), &
                                  nw_angle(2,j), atom_type(nw_angle(2,j)), &
                                  nw_angle(3,j), atom_type(nw_angle(3,j)),nw_angle_val(j),j+nw_bonds
           endif
       enddo
       close(101)
    enddo
    print '(A)','**'
    print '(A)','** Using keyword geometry  autoz in NWchem'
    print '(A)','** geometry adjust'
    print '(A)','**    zcoord'
    print '(A)','**       bond angle torsion'
    print '(A)','**    end'
    print '(A)','** end'
    print '(A)','**'

    if(allocated(ang_name)) deallocate(ang_name)
    return
end
!
!-------------------------------------------------------------------------
!       Use gulp fit the Molecular Mechanic Potententials (torsiones)    -
!-------------------------------------------------------------------------
!

subroutine fit_torsions
    !use modfim,only: n_mol
    use modnei,only: mang,mdih
    !use mode,  only: x
    use mod_nw,only:  nw_dihe, &
                      nw_dihe_val,atom_type, nw_bonds,nw_angles,nw_dihes
    implicit none

    integer                       i,j,k
    integer                       tors_typs
    character(5),allocatable::   tors_name(:,:)
    character(8)             ::   tors_file = '        ' 
    character(4)             ::   con = '    '
    logical                       new

    allocate(tors_name(4,nw_dihes))

    tors_typs  = 0

    do i=1, nw_dihes
       if(i==1)then
           tors_typs = 1
           tors_name(1,tors_typs) = trim(adjustl( atom_type(nw_dihe(1,i)) ))
           tors_name(2,tors_typs) = trim(adjustl( atom_type(nw_dihe(2,i)) ))
           tors_name(3,tors_typs) = trim(adjustl( atom_type(nw_dihe(3,i)) ))
           tors_name(4,tors_typs) = trim(adjustl( atom_type(nw_dihe(4,i)) ))
       else 
           new = .true.
           do k=1,tors_typs
              if(  tors_name(1,k)==atom_type(nw_dihe(1,i)) .and. &
                 & tors_name(2,k)==atom_type(nw_dihe(2,i)) .and. &
                 & tors_name(3,k)==atom_type(nw_dihe(3,i)) .and. &
                 & tors_name(4,k)==atom_type(nw_dihe(4,i)) ) then
                 new = .false.
                 !print *,tors_name(1,k), atom_type(nw_dihe(1,i))
                 !print *,tors_name(2,k), atom_type(nw_dihe(2,i))
                 !print *,tors_name(3,k), atom_type(nw_dihe(3,i))
                 !print *,tors_name(4,k), atom_type(nw_dihe(4,i))
              elseif(tors_name(4,k)==atom_type(nw_dihe(1,i)) .and. &
                    & tors_name(3,k)==atom_type(nw_dihe(2,i)) .and. &
                    & tors_name(2,k)==atom_type(nw_dihe(3,i)) .and. &
                    & tors_name(1,k)==atom_type(nw_dihe(4,i)) ) then
                 new = .false.
              endif
           enddo
           if( new ) then
              tors_typs = tors_typs + 1
              tors_name(1,tors_typs) = trim(adjustl( atom_type(nw_dihe(1,i)) ))
              tors_name(2,tors_typs) = trim(adjustl( atom_type(nw_dihe(2,i)) ))
              tors_name(3,tors_typs) = trim(adjustl( atom_type(nw_dihe(3,i)) ))
              tors_name(4,tors_typs) = trim(adjustl( atom_type(nw_dihe(4,i)) ))
              !print *,tors_name(1,tors_typs) ,tors_name(2,tors_typs) ,tors_name(3,tors_typs) ,tors_name(4,tors_typs) 
              !print *,new
           endif
       endif
    enddo
    print '(A31,I5)', 'Number of torsion angle types: ', tors_typs


    do i=1,tors_typs
       !print *,ang_name(1,i),ang_name(2,i),ang_name(3,i)
       tors_file(1:5) = 'tors_'
       write(con,'(I4)') i
       con = adjustl(con)
       if (i<10) tors_file(6:6) = con
       if (i>=10 .and. i<100) tors_file(6:7) = con
       tors_file = trim(adjustl(tors_file))
       open(101,file=tors_file,status='unknown')
       write(101,'(A18)') 'Key words in gulp:'
       write(101,'(A18)') 'torsion bond intra'
       write(101,'(4(A5,A6),A2,A2,A5)') tors_name(1,i),' core ',tors_name(2,i),' core ',tors_name(3,i),& 
                     & ' core ',tors_name(4,i),' core ', 'k ','n ','phi0 '
       do j=1,nw_dihes
           if(tors_name(1,i)==atom_type(nw_dihe(1,j)) .and. &
              & tors_name(2,i)==atom_type(nw_dihe(2,j)) .and. &
              & tors_name(3,i)==atom_type(nw_dihe(3,j)) .and. &
              & tors_name(4,i)==atom_type(nw_dihe(4,j))  ) then
              write(101,'(A4,4(I5,1X,A5),F10.5,I4)')'Tor ',nw_dihe(1,j),atom_type(nw_dihe(1,j)), &
                                  nw_dihe(2,j), atom_type(nw_dihe(2,j)), &
                                  nw_dihe(3,j), atom_type(nw_dihe(3,j)), &
                                  nw_dihe(4,j), atom_type(nw_dihe(4,j)), &
                                  nw_dihe_val(j),j+nw_bonds+nw_angles
           elseif(tors_name(4,i)==atom_type(nw_dihe(1,j)) .and. &
              & tors_name(3,i)==atom_type(nw_dihe(2,j)) .and. &
              & tors_name(2,i)==atom_type(nw_dihe(3,j)) .and. & 
              & tors_name(1,i)==atom_type(nw_dihe(4,j))    ) then
              write(101,'(A4,4(I5,1X,A5),F10.5,I4)')'Tor ',nw_dihe(1,j),atom_type(nw_dihe(1,j)), &
                                  nw_dihe(2,j), atom_type(nw_dihe(2,j)), &
                                  nw_dihe(3,j), atom_type(nw_dihe(3,j)), &
                                  nw_dihe(4,j), atom_type(nw_dihe(4,j)), &
                                  nw_dihe_val(j),j+nw_bonds+nw_angles
           endif
       enddo
       close(101)
    enddo
    !print '(A)','**'
    !print '(A)','** Using keyword geometry  autoz in NWchem'
    !print '(A)','** geometry adjust'
    !print '(A)','**    zcoord'
    !print '(A)','**       torsion'
    !print '(A)','**    end'
    !print '(A)','** end'
    !print '(A)','**'

    deallocate(tors_name)

    return
end
