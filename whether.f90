!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                     =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!

!@@@@@@@@@@@@@@@@@@@@@@       analysing  angle    @@@@@@@@@@@@@@@@@@@@@@

subroutine w_angle
    !use modfim,only: n_mol
    use modnei,only: n_ang,angle,mang,ang_value
    use mode,  only: x!, a, b, c
    use mod_nw,only: atom_type
    implicit none

    integer         i,iatom,jatom,katom
    !integer         n_ang
    real(8)         r1(3),r2(3),radi,rsq1,rsq2
    real(8)         ang,cos_ang,ave_ang
    real(8)      :: pi=3.1415926535897932

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
       ang = ang*180.0/pi

       if (ang_value(i)-ang>0.01 .or. ang_value(i)-ang<-0.01) then
          iatom = angle(1,i)
          jatom = angle(2,i)
          katom = angle(3,i)
          print '(A7,3(I4,1X,A7,1X),F8.4,1X,F8.4)','Angle: ',iatom,atom_type(iatom),jatom,atom_type(jatom),& 
                 katom,atom_type(katom),ang_value(i),ang
       endif 
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

subroutine w_bond
    use modnei,only: nbond,bond,bond_value
    use mode,  only: x!, a, b, c
    use mod_nw,only: atom_type
    implicit none

    integer         i,iatom,jatom
    real(8)         r1(3),radi,rsq1


    do i=1,nbond
       r1(1)=x(1,bond(1,i))-x(1,bond(2,i))
       r1(2)=x(2,bond(1,i))-x(2,bond(2,i))
       r1(3)=x(3,bond(1,i))-x(3,bond(2,i))

       rsq1= r1(1)**2+r1(2)**2+r1(3)**2
       radi= sqrt(rsq1)
       if (bond_value(i)-radi>0.001 .or. bond_value(i)-radi<-0.001) then
          iatom = bond(1,i)
          jatom = bond(2,i)
          print '(A6,2(I4,1X,A7,1X),F10.6,1X,F10.6)','Bond: ',iatom,atom_type(iatom),jatom,atom_type(jatom), &
                                         bond_value(i),radi 
       endif 
    enddo

    return
end

subroutine w_torsion
    !use modfim,only: n_mol
    use modnei,only: n_dih,dihedral,dih_value
    !use mode,  only: x!, a, b, c
    use mod_nw,only: atom_type
    implicit none

    integer            i
    integer            iatom,jatom,katom,latom
    real(8),external:: torsion_angle
    real(8)            dihv

    !n_dih = 0
    ! determin dihedral angles
    !print *,nbond,' bonds'
    do i=1,n_dih 
      iatom = dihedral(1,i)
      jatom = dihedral(2,i)
      katom = dihedral(3,i) 
      latom = dihedral(4,i)
      dihv = torsion_angle(iatom,jatom,katom,latom)

      !print *,dihv,dih_value(i)
      if (dih_value(i)-dihv>0.01 .or. dih_value(i)-dihv<-0.01) then
         print '(A9,4(I4,1X,A7,1X),F8.4,1X,F8.4)','Torsion: ', iatom,atom_type(iatom),&
                                jatom,atom_type(jatom), &
                   katom,atom_type(katom),latom,atom_type(latom), dih_value(i),dihv
      endif 
    enddo

    !print *,n_dih,' dihedral angles.'
    return
end
