!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!

subroutine dihedrals
    !use modfim,only: n_mol
    use modnei,only: n_ang,n_dih,angle,dihedral,mang,mdih
    !use mode,  only: x!,npart !, a, b, c
    implicit none


    integer         i,j
    !integer         n_atom
    !real(8)      :: pi=3.1415926535897932
    real(8)         dihe
    !logical,allocatable:: bond_no_dih(:,:)
    real(8),external:: torsion_angle

!  Determining dihedrals
    !allocate(bond_no_dih(nbond,nbond))

    !bond_no_dih(:,:) = .true.
    n_dih = 0
    do i=1,n_ang-1
       do j=i+1,n_ang  
          if(n_dih>=mdih)call err('Error: Too many dihedrals, Modify mdih and retry!')
          if(angle(2,i)==angle(1,j).and.angle(2,j)==angle(1,i))then
            !if(bond_no_dih(angle(1,i),angle(2,i)))then
                !bond_no_dih(angle(1,i),angle(2,i)) = .false.
                !bond_no_dih(angle(2,i),angle(1,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(3,j)
                dihedral(2,n_dih)=angle(1,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(3,i)
            !endif
          elseif(angle(2,i)==angle(3,j).and.angle(2,j)==angle(1,i))then
            !if(bond_no_dih(angle(1,i),angle(2,i)))then
                !bond_no_dih(angle(1,i),angle(2,i)) = .false.
                !bond_no_dih(angle(2,i),angle(1,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(1,j)
                dihedral(2,n_dih)=angle(1,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(3,i)
            !endif
          elseif(angle(2,i)==angle(1,j).and.angle(2,j)==angle(3,i))then
            !if(bond_no_dih(angle(2,i),angle(3,i)))then
                !bond_no_dih(angle(2,i),angle(3,i)) = .false.
                !bond_no_dih(angle(3,i),angle(2,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(3,j)
                dihedral(2,n_dih)=angle(3,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(1,i)
            !endif
          elseif(angle(2,i)==angle(3,j).and.angle(2,j)==angle(3,i))then
            !if(bond_no_dih(angle(2,i),angle(3,i)))then
                !bond_no_dih(angle(2,i),angle(3,i)) = .false.
                !bond_no_dih(angle(3,i),angle(2,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(1,j)
                dihedral(2,n_dih)=angle(3,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(1,i)
            !endif
          endif
       enddo
    enddo
    ! still have problems
    !print '(A56)','** The determing of dihedral angle still have problems, '
    !print '(A29)','** you need set it manuly ...'
    do i=1,n_dih
       dihe  = torsion_angle(dihedral(1,i),dihedral(2,i),dihedral(3,i),dihedral(4,i))
       !print *,dihe
    enddo

    return
end
