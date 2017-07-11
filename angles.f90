!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!

!@@@@@@@@@@@@@@@@@@@@@@       analysing  angle    @@@@@@@@@@@@@@@@@@@@@@

subroutine angles
    use modnei,only: nbond,bond,n_ang,angle,mang,ang_value
    use mode,  only: x!, a, b, c
    implicit none

    integer         i,j
    !integer         n_ang
    real(8)         r1(3),r2(3),radi,rsq1,rsq2
    real(8)         ang,cos_ang,ave_ang
    real(8)      :: pi=3.1415926535897932

    n_ang = 0
    do i=1,nbond-1
       do j=i+1,nbond
          if(n_ang>=mang)call err('Error: Too many angles, Modify mang and retry!')
          if(bond(1,i)==bond(1,j))then
             n_ang=n_ang+1
             angle(2,n_ang)=bond(1,i)
             angle(1,n_ang)=bond(2,i)
             angle(3,n_ang)=bond(2,j)
          elseif(bond(1,i)==bond(2,j))then
             n_ang=n_ang+1
             angle(2,n_ang)=bond(1,i)
             angle(1,n_ang)=bond(2,i)
             angle(3,n_ang)=bond(1,j)
          elseif(bond(2,i)==bond(1,j))then
             n_ang=n_ang+1
             angle(2,n_ang)=bond(2,i)
             angle(1,n_ang)=bond(1,i)
             angle(3,n_ang)=bond(2,j)
          elseif(bond(2,i)==bond(2,j))then
             n_ang=n_ang+1
             angle(2,n_ang)=bond(2,i)
             angle(1,n_ang)=bond(1,i)
             angle(3,n_ang)=bond(1,j)
          endif
       enddo
    enddo

    ave_ang = 0.0
    ang_value = 0
    do i=1,n_ang
       r1(1)=x(1,angle(1,i))-x(1,angle(2,i))
       r1(2)=x(2,angle(1,i))-x(2,angle(2,i))
       r1(3)=x(3,angle(1,i))-x(3,angle(2,i))

       r2(1)=x(1,angle(3,i))-x(1,angle(2,i))
       r2(2)=x(2,angle(3,i))-x(2,angle(2,i))
       r2(3)=x(3,angle(3,i))-x(3,angle(2,i))

       rsq1= r1(1)**2+r1(2)**2+r1(3)**2
       rsq2= r2(1)**2+r2(2)**2+r2(3)**2
       radi= sqrt(rsq1*rsq2)

       cos_ang=r1(1)*r2(1)+r1(2)*r2(2)+r1(3)*r2(3)
       cos_ang=cos_ang/radi
       ang    = acos(cos_ang)
       ang_value(i) = ang*180.0/pi
       ave_ang = ave_ang + ang
       !write(45,'(4I5,3F12.8,3I5)') angle(1,i),angle(2,i),angle(3,i),0,bond_length(angle(1,i),angle(2,i)),ang,0.0,0,0,0
       
    enddo
    !ave_ang = ave_ang/n_ang
    !ave_ang = 180.0*ave_ang/pi
    !print '(A2)', '**'
    !print '(A24,F10.5)', '** average angle value: ',ave_ang
    !print '(A2)', '**'

    return
end

!@@@@@@@@@@@@@@@@@@@@@@       analysing  bond    @@@@@@@@@@@@@@@@@@@@@@

subroutine bonds
    use modnei,only: nbond,bond,bond_value
    use mode,  only: x!, a, b, c
    implicit none

    integer         i
    real(8)         r1(3),radi,rsq1


    do i=1,nbond
       r1(1)=x(1,bond(1,i))-x(1,bond(2,i))
       r1(2)=x(2,bond(1,i))-x(2,bond(2,i))
       r1(3)=x(3,bond(1,i))-x(3,bond(2,i))

       rsq1= r1(1)**2+r1(2)**2+r1(3)**2
       radi= sqrt(rsq1)
       bond_value(i) = radi
    enddo

    return
end

subroutine torsions
    !use modfim,only: n_mol
    use modnei,only: n_dih,dihedral,dih_value,bond,nbond,atable
    !use mode,  only: x!, a, b, c
    implicit none

    integer            i,j,k
    integer            iatom,jatom,katom,latom
    real(8),external:: torsion_angle


    n_dih = 0
    ! determin dihedral angles
    !print *,nbond,' bonds'
    do i=1,nbond 
       jatom = bond(1,i)
       katom = bond(2,i)
       if (atable(1,jatom)>1.and.atable(1,katom)>1) then
          do j=2,atable(1,jatom)+1
             iatom = atable(j,jatom)
             if(iatom/=katom)then
                do k=2,atable(1,katom)+1
                   latom = atable(k,katom)
                   if(latom/=jatom) then
                      n_dih = n_dih + 1
                      dihedral(1,n_dih) = iatom
                      dihedral(2,n_dih) = jatom
                      dihedral(3,n_dih) = katom
                      dihedral(4,n_dih) = latom
                      dih_value(n_dih)  = torsion_angle(iatom,jatom,katom,latom)
                   endif
                enddo
             endif
          enddo
       endif
    enddo

    !print *,n_dih,' dihedral angles.'
    return
end
