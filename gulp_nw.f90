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
!          Write FIELD file for Molecular Mechanic Simulations           -
!-------------------------------------------------------------------------
!

subroutine nw_info
    !use modfim,only: n_mol
    use modnei,only: nbond,bond,n_ang,mang,mdih,atable
    use mode,  only: x,npart,types,ctype
    use mod_nw,only: nw_bond, nw_angle, nw_dihe, nw_bond_val, nw_angle_val, nw_dihe_val,&
                      atom_type, nw_angles, nw_bonds, nw_dihes,latmtyp
    implicit none

    integer         i,j,k,ti,l,m,n,o,p
    real(8),external:: torsion_angle

    character(7)         tc,tcc
    character(100)       line

    !logical           :: dihe_zero = .false.

    call  neighbor
    call  findmole

!  write gulp force filed, using NWchem output

!  Determining Angles

    n_ang = 0
    !print *,nbond
    if (latmtyp)then
       call atm_typ
    else
       if(.not. allocated(atom_type)) allocate(atom_type(npart))
       atom_type  = 'null '

       do i = 1,npart
          atom_type(i) = ctype(types(i))
       enddo
    endif
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   Read in nwchem output of optimized z-matrix
    open(33,file='nw.log',status='old')
    !print *,'** you need autoz options in nwchem input file!'
    !do i=1,7
       !read(33,*)
    !enddo

    nw_bonds = 0
    nw_angles= 0
    nw_dihes = 0
    do
        read(33,'(A100)',err=99,end=99) line
        if (index(line,'Stretch').NE.0 .or. index(line,'Bend').NE.0 &
            .or. index(line,'Torsion').NE.0 ) then
           read(line,*)ti,tc
           if (tc=='Stretch') nw_bonds= nw_bonds + 1 
           if (tc=='Bend')    nw_angles= nw_angles + 1
           if (tc=='Torsion') nw_dihes= nw_dihes + 1
        endif
    enddo

99  continue
    print '(A2)','**'
    print '(A2,I5,A6)','**',nw_bonds,' bonds'
    print '(A2,I5,A7)','**',nw_angles,' angles'
    print '(A2,I5,A16)','**',nw_dihes,' dihedral angles'
    print '(A25)','** in NWchem output file.'
    print '(A2)','**'
!   
!   prepare memory
    allocate(nw_bond(2,nw_bonds))
    allocate(nw_angle(3,nw_angles))
    allocate(nw_dihe(4,nw_dihes))
    allocate(nw_bond_val(nw_bonds))
    allocate(nw_angle_val(nw_angles))
    allocate(nw_dihe_val(nw_dihes))

!
    rewind(33)
    !do i=1,7
       !read(33,*)
    !enddo
    i=0
    j=0
    k=0
    do
        read(33,'(A100)',err=9,end=9) line
        if (index(line,'Stretch').NE.0 .or. index(line,'Bend').NE.0 &
            .or. index(line,'Torsion').NE.0 ) then
           read(line,*) ti,tc
           if (tc=='Stretch') then
               i = i + 1
               read(line,*)ti,tcc,nw_bond(1,i),nw_bond(2,i),nw_bond_val(i)
           elseif (tc=='Bend') then
               j = j + 1
               read(line,*)ti,tcc,nw_angle(1,j),nw_angle(2,j),nw_angle(3,j),nw_angle_val(j)
           elseif (tc=='Torsion') then
               k = k + 1
               read(line,*)ti,tcc,nw_dihe(1,k),nw_dihe(2,k),nw_dihe(3,k),nw_dihe(4,k),nw_dihe_val(k)
           endif
        endif
    enddo
9   continue
    if(i/=nw_bonds) print '(A)','** bond read error!'
    if(j/=nw_angles) print '(A)','** angle read error!'
    if(k/=nw_dihes) print '(A)','** dihedral read error!'
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@       analysing  bond    @@@@@@@@@@@@@@@@@@@@@@


    open(34,file='nw_cn.txt',status='unknown')
    open(40,file='nw_ch.txt',status='unknown')
    open(41,file='nw_no.txt',status='unknown')
    k = 0
    l = 0
    m = 0
    do i=1,nw_bonds
       if( (types(nw_bond(1,i))==1.and.types(nw_bond(2,i))==4) .or. &
         &  (types(nw_bond(1,i))==4.and.types(nw_bond(2,i))==1) )then ! CN
          k = k + 1
          !@@@       topology analysis      @@@@
          if ( types(nw_bond(1,i))==1) then
              o=atable(1,nw_bond(1,i))
              line = ' '
              line(1:10) = ctype(1)//' bond to '
              n=12
              do p=2,o+1
                 !print *,ctype(types((atable(p,nw_bond(1,i)))))
                 ti = len(ctype(types((atable(p,nw_bond(1,i))))))
                 line(n:n+ti) = ctype(types((atable(p,nw_bond(1,i)))))
                 n = n + ti
              enddo
          endif
          write(34,*) k, nw_bond_val(i),line ! for plot
       elseif( (types(nw_bond(1,i))==1.and.types(nw_bond(2,i))==2) .or. &
         &  (types(nw_bond(1,i))==2.and.types(nw_bond(2,i))==1) )then ! CH
          l = l + 1
          write(40,*) l, nw_bond_val(i)! for plot
       elseif( (types(nw_bond(1,i))==3.and.types(nw_bond(2,i))==4) .or. &
         &  (types(nw_bond(1,i))==4.and.types(nw_bond(2,i))==3) )then ! NO
          m = m + 1
          write(41,*) m, nw_bond_val(i)! for plot
       endif
    enddo
    close(34)
    close(40)
    close(41)
    call plot('nw_cn.txt','bondNumb','bondlength','CN','cn.eps')
    call plot('nw_ch.txt','bondNumb','bondlength','CH','ch.eps')
    call plot('nw_no.txt','bondNumb','bondlength','NO','no.eps')


!@@@@@@@@@@@@@@@@@@@@@@       analysing  angle    @@@@@@@@@@@@@@@@@@@@@@

    call angles
    !print '(A11)','All Angles:'
    !do i=1,nw_angles 
       !print '(3A5)',atom_type(nw_angle(1,i)),atom_type(nw_angle(2,i)),atom_type(nw_angle(3,i))
       !print '(3A5)',atom_type(angle(1,i)),atom_type(angle(2,i)),atom_type(angle(3,i))   
    !enddo

    !call write_dlfield
    ! angle C_3 N_3 C_3
    k = 0
    open(40,file='nw_CNC.txt',status='unknown')    ! C_3 N_3 C_3
    do i=1,nw_angles 
       !print '(A5)',atom_type(nw_angle(2,i))
       if(atom_type(nw_angle(1,i))=='C3'.and.atom_type(nw_angle(2,i))=='N3' &
          & .and.atom_type(nw_angle(3,i))=='C3') then
           k = k + 1
           write(40,'(I5,F10.5)') k, nw_angle_val(i)
           !write(40,'(I5,F10.5)') k, ang_value(i)
       endif
    enddo
    close(40)
    call plot('nw_CNC.txt','AngNumb','AngValue','CNC','cnc.eps')

    ! angle H_ N_3 H_
    k = 0
    open(40,file='nw_HCH.txt',status='unknown')    ! H_  C_3 H_
    do i=1,nw_angles 
       !print '(A5)',atom_type(nw_angle(2,i))
       if(atom_type(nw_angle(1,i))=='H'.and.atom_type(nw_angle(2,i))=='C3' &
          & .and.atom_type(nw_angle(3,i))=='H') then
           k = k + 1
           write(40,'(I5,F10.5)') k, nw_angle_val(i)
       endif
    enddo
    close(40)
    call plot('nw_HCH.txt','AngNumb','AngValue','HCH','hch.eps')

    ! angle N N O
    k = 0
    open(40,file='nw_NNO.txt',status='unknown')    ! NNO
    do i=1,nw_angles 
       !print '(A5)',atom_type(nw_angle(2,i))
       if(atom_type(nw_angle(1,i))=='N3'.and.atom_type(nw_angle(2,i))=='N32' &
          & .and.atom_type(nw_angle(3,i))=='O2') then
           k = k + 1
           write(40,'(I5,F10.5)') k, nw_angle_val(i)
       endif
    enddo
    close(40)
    call plot('nw_NNO.txt','AngNumb','AngValue','NNO','nno.eps')


    ! angle O N O
    k = 0
    open(40,file='nw_ONO.txt',status='unknown')    ! NNO
    do i=1,nw_angles 
       !print '(A5)',atom_type(nw_angle(2,i))
       if(atom_type(nw_angle(1,i))=='O2'.and.atom_type(nw_angle(2,i))=='N32' &
          & .and.atom_type(nw_angle(3,i))=='O2') then
           k = k + 1
           write(40,'(I5,F10.5)') k, nw_angle_val(i)
       endif
    enddo
    close(40)
    call plot('nw_ONO.txt','AngNumb','AngValue','NNO','ono.eps')

!@@@@@@@@@@@@@@@@@@@@@@       analysing  angle    @@@@@@@@@@@@@@@@@@@@@@
! dreiding harmonic/torsion bond kcal dreiding single/three-body bond/
! torsion bond kcal dreiding double/torsion bond kcal dreiding resonant
! epsilon/inversion
! lennard
    open(42,file='atoms.xyz',status='unknown')
    write(42,*)npart
    write(42,*)
    do i=1,npart
       write(42,*)atom_type(i),x(1,i),x(2,i),x(3,i)
    enddo
    close(42)

    !!! input for gulp
    open(42,file='inp',status='unknown')
    write(42,'(A41)')'opti conse qok nomod pres pot efg noauto '
    write(42,'(A5)') 'title'
    write(42,'(A16)') 'GULP calculation'
    write(42,'(A3)') 'end'
    write(42,'(A16)') 'library dreiding'
    write(42,'(A13)') 'dump gulp.grs'
    write(42,'(A21)') 'output movie xyz gulp' 
    write(42,*)
    write(42,'(A7)')  'Species'
    write(42,'(A12)') 'C3 core C3'
    write(42,'(A12)') 'H  core H '
    write(42,'(A12)') 'N3 core N3'
    write(42,'(A12)') 'O2 core O2'
    write(42,*)
    write(42,'(A9)') 'cartesian'
    do i=1,npart
       write(42,'(A5,A6,3F12.8,I3,2F10.5,3I3)')atom_type(i),' core ',x(1,i),x(2,i),x(3,i),0,1.0,0.0,1,1,1
    enddo
    do i=1,nbond
       write(42,'(A8,2I8)')'connect ',bond(1,i),bond(2,i)
    enddo
    close(42)

    call fit_bondes
    call fit_angles
    call fit_torsions

    return
end


