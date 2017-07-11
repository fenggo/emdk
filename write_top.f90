!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!  ------             WRITE topology file for Gromacs      ------        -
!-------------------------------------------------------------------------
!

subroutine write_top(filename)
    use mode,only: npart,types,ctype, mass,q
    use mod_nw,only: atom_type,latmtyp  
    use modnei,only: nbond,bond,n_ang,angle,n_dih,dihedral
    implicit none
!
!   Local Variables
!
    integer          i,j,k
    real(8)          pi
    character(*)     filename
    character(6),allocatable::     atomtype(:)
    logical          hav

    pi = dacos(-1.0d0)
!
    call  neighbor
    call  findmole

    allocate(atomtype(npart))
    
    if (latmtyp) then
       call atm_typ
    else
       if(.not. allocated(atom_type)) allocate(atom_type(npart))
       do j=1,npart
          atom_type(j) = ctype(types(j))
       enddo
    endif

    open(36,file=filename,status='unknown') !,access='append'
    write(36,'(A21)')';  GROMACS FF for HMX'
    write(36,*)
    write(36,*)
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,'(A21)')'; Intermolecular part'
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,*)
    write(36,*)
    write(36,'(A12)')'[ defaults ]'
    write(36,'(A50)')'; nbfunc   comb-rule   gen-pairs  fudgeLJ  fudgeQQ'
    write(36,'(A50)')'      1       3           no          1.0     0.0'
    write(36,*)
    write(36,*)
    write(36,'(A13)')'[ atomtypes ]'
    write(36,'(A57)')'; name atomnb mass charge ptype sigma(nm) epsilon(kJ/mol)'

    j = 0
    do i =1,npart
       if (j>=1) then
           hav = .false.
           do k = 1,j
              if (atom_type(i) == atomtype(k)) hav = .true.
           enddo
           
           if (.not.hav) then
              j = j + 1
              atomtype(j) = (atom_type(i)) 
              if (ctype(types(i))=='H') write(36,'(1X,A5,I3,A42)')atomtype(j),1,'   1.0079   0.2433   A    0.2420     0.1255 '
              if (ctype(types(i))=='C') write(36,'(1X,A5,I3,A42)')atomtype(j),6,'  12.0110  -0.480    A    0.3750     0.43932 '
              if (ctype(types(i))=='O') write(36,'(1X,A5,I3,A42)')atomtype(j),8,'  15.9994  -0.270    A    0.2960     0.87864 '
              if (ctype(types(i))=='N') write(36,'(1X,A5,I3,A42)')atomtype(j),7,'  14.0067   0.280    A    0.3250     0.71128 '

              !print *,j,' : ',atomtype(j)
           endif 
       endif

       if(j==0)then
           j = 1
           atomtype(j) = atom_type(i) 
              if (ctype(types(i))=='H') write(36,'(1X,A5,I3,A42)')atomtype(j),1,'   1.0079   0.2433   A    0.2420     0.1255 '
              if (ctype(types(i))=='C') write(36,'(1X,A5,I3,A42)')atomtype(j),6,'  12.0110  -0.480    A    0.3750     0.43932 '
              if (ctype(types(i))=='O') write(36,'(1X,A5,I3,A42)')atomtype(j),8,'  15.9994  -0.270    A    0.2960     0.87864 '
              if (ctype(types(i))=='N') write(36,'(1X,A5,I3,A42)')atomtype(j),7,'  14.0067   0.280    A    0.3250     0.71128 '
       endif
       
       if (ctype(types(i))=='H') q(i) =  0.2433    
       if (ctype(types(i))=='C') q(i) = -0.480    
       if (ctype(types(i))=='O') q(i) = -0.270    
       if (ctype(types(i))=='N') q(i) =  0.280    
    enddo

    write(36,*)
    write(36,*)
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,'(A49)')'; Intramolecular part (to be fitted from QM data)'
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,*)
    write(36,*)
    write(36,'(A16)')'[ moleculetype ]'
    write(36,'(A15)')'; Name      HMX'
    write(36,'(A13)')'  hmx       3'
    write(36,*)
    write(36,*)
    write(36,'(A9)')'[ atoms ]'
    write(36,'(A53)')';   nr   type  resnr residue  atom  cgnr charge  mass'
    do i=1,npart
        write(36,'(I7,A6,I5,A6,A6,I5,F10.6,F12.4)') i,atom_type(i),1,'hmx',atom_type(i),i,0.0,mass(types(i))
    enddo

    write(36,*)
    write(36,*)

    call bonds
    write(36,'(A9)') '[ bonds ]'
    write(36,'(A36)')';  ai    aj   type    r0          kr'
    do i=1,nbond
        write(36,'(2I7,I7,F10.6,F12.4)') bond(1,i),bond(2,i), 1, 0.0,0.0 ! 3=morse,1=harm
    enddo

    write(36,*)
    write(36,*)

    call angles
    write(36,'(A10)') '[ angles ]'
    write(36,'(A46)')'; ai    aj   ak   type     theta0           kt'
    do i=1,n_ang
        write(36,'(3I7,I7,F10.6,F12.4)') angle(1,i),angle(2,i),angle(3,i), 1,0.0, 0.0 ! 3=morse,1=harm
    enddo
    write(36,*)
    write(36,*)
    
    call torsions
    write(36,'(A13)') '[ dihedrals ]'
    write(36,'(A50)')';  ai    aj    ak    al   type     phi    k   mult'
    do i=1,n_dih
        write(36,'(4I7,I7,F10.6,F12.4,I4)') dihedral(1,i),dihedral(2,i),dihedral(3,i), dihedral(4,i),1,0.0, 0.0, 2! 3=morse,1=harm
    enddo

    write(36,*)
    write(36,*)
 
    write(36,'(A10)') '[ system ]'
    write(36,'(A6)') '; Name'
    write(36,'(A5)') '  hmx'
    write(36,*)
    write(36,*)
    write(36,'(A13)') '[ molecules ]'
    write(36,'(A22)') '; Compound       #mols'
    write(36,'(A19)') '   hmx            1'
    write(36,*)
    write(36,*)


    close(36)

    return
end

!-------------------------------------------------------------------------c
!  ------           prepare coordinate in gromacs format      ------      -
!-------------------------------------------------------------------------c
!

subroutine write_gro
    use mode,   only: x,v,npart,types,a,b,c,r1,r2,r3,alpha,gamm,beta,nelem, alive,mass
    use modctrl,only: iwr
    use mod_nw,only: atom_type !,latmtyp 
    use modfim,only: n_mol,molecule
    implicit none
!
!   Local Variables
!
    integer         i,np,j
    real(8)         cy
    real(8)         pi, alph, bet, gam
!
    pi = acos(-1.0d0)

    call  neighbor
    call  findmole

    call proj_coord
    call atm_typ
!  For gromacs: 

    gam = gamm*pi/180
    bet = beta*pi/180
    alph= alpha*pi/180

    a = (/r1,0.0d0,0.0d0/)
    b = (/r2*dcos(gam),r2*dsin(gam),0.0d0/)
    cy= r3*dcos(alph)-r3*dcos(bet)*dcos(gam)
    cy= cy/dsin(gam)
    c = (/r3*dcos(bet),cy,dsqrt(r3*r3*dsin(bet)*dsin(bet)-cy*cy)/)

    call cart_coord

    open(31,file='coord.gro',status='unknown')
    write(31,'(A16)')'Coordinate file for gromacs chemical simulaitons'
    write(31,'(I5)')npart 

    np = 0
    do i=1,npart
         if(alive(i)/=0) then
              np = np + 1                  !molecule(i)
              write(31,'(i5,2a5,i5,3f8.3)')molecule(i),'OCNM ',atom_type(i),np,x(1,i)/10.0,x(2,i)/10.0,x(3,i)/10.0
         endif
    end do
    write(31,'(3f10.5)')a(1)/10.0,b(2)/10.0,c(3)/10.0
    if(b(1)>0.1.or.b(1)<-0.01.or.c(1)>0.01.or.c(1)<-0.01.or.c(2)>0.01.or.c(2)<-0.01) then
         write(31,'(3F13.6)')b(1)/10.0,c(1)/10.0,c(2)/10.0
    endif

!
!   shock velocity
!      if(iwr==19)then  
!         write(31,*)
!         write(31,'(A10)')'Velocities'
!         write(31,*)
!         do i=1,npart
!             write(31,*)i,(v(j,i),j=1,3)
!         enddo
!      endif
!

   close(31)
!
   return
end

!-------------------------------------------------------------------------c
!  ------          prepare topological in gromacs format      ------      -
!-------------------------------------------------------------------------c
!

subroutine write_gromacs_top(filename)
    use mode,only: npart,types,ctype, mass, q
    use mod_nw,only: atom_type!,latmtyp  
    use modnei,only: nbond,bond,n_ang,angle,n_dih,dihedral
    use modfim,only: n_mol,molecule
    implicit none
!
!   Local Variables
!
    integer          i,j,k
    integer          ia,ja,ka,la
    real(8)          pi
    real(8)          r0,kr,theta0,kt,dihe

    character(*)     filename
    character(6),allocatable::     atomtype(:)
    logical          hav

    real(8),external:: torsion_angle


    pi = dacos(-1.0d0)
!
    !call  neighbor
    !call  findmole

    allocate(atomtype(npart))
    

    open(36,file=filename,status='unknown') !,access='append'
    write(36,'(A30)')';  GROMACS FF for Nitromethane'
    write(36,*)
    write(36,*)
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,'(A21)')'; Intermolecular part'
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,*)
    write(36,*)
    write(36,'(A12)')'[ defaults ]'
    write(36,'(A50)')'; nbfunc   comb-rule   gen-pairs  fudgeLJ  fudgeQQ'
    write(36,'(A50)')'      1       3          yes          1.0     0.0'
    write(36,*)
    write(36,*)
    write(36,'(A13)')'[ atomtypes ]'
    write(36,'(A57)')'; name atomnb mass charge ptype sigma(nm) epsilon(kJ/mol)'

    j = 0
    do i =1,npart
       if (j>=1) then
           hav = .false.
           do k = 1,j
              if (atom_type(i) == atomtype(k)) hav = .true.
           enddo
           
           if (.not.hav) then
              j = j + 1
              atomtype(j) = (atom_type(i)) 
              if (ctype(types(i))=='H') write(36,'(1X,A5,I3,A42)')atomtype(j),1,'   1.0079   0.2433   A    0.2420     0.1255 '
              if (ctype(types(i))=='C') write(36,'(1X,A5,I3,A42)')atomtype(j),6,'  12.0110  -0.480    A    0.3750     0.43932 '
              if (ctype(types(i))=='O') write(36,'(1X,A5,I3,A42)')atomtype(j),8,'  15.9994  -0.270    A    0.2960     0.87864 '
              if (ctype(types(i))=='N') write(36,'(1X,A5,I3,A42)')atomtype(j),7,'  14.0067   0.280    A    0.3250     0.71128 '

              !print *,j,' : ',atomtype(j)
           endif 
       endif

       if(j==0)then
           j = 1
           atomtype(j) = atom_type(i) 
           if (ctype(types(i))=='H') write(36,'(1X,A5,I3,A42)')atomtype(j),1,'   1.0079   0.2433   A    0.2420     0.1255 '
           if (ctype(types(i))=='C') write(36,'(1X,A5,I3,A42)')atomtype(j),6,'  12.0110  -0.480    A    0.3750     0.43932 '
           if (ctype(types(i))=='O') write(36,'(1X,A5,I3,A42)')atomtype(j),8,'  15.9994  -0.270    A    0.2960     0.87864 '
           if (ctype(types(i))=='N') write(36,'(1X,A5,I3,A42)')atomtype(j),7,'  14.0067   0.280    A    0.3250     0.71128 '
       endif
       
       if (ctype(types(i))=='H') q(i) =  0.2433    
       if (ctype(types(i))=='C') q(i) = -0.480    
       if (ctype(types(i))=='O') q(i) = -0.270    
       if (ctype(types(i))=='N') q(i) =  0.280    
    enddo


    write(36,*)
    write(36,*)
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,'(A49)')'; Intramolecular part (to be fitted from QM data)'
    write(36,'(A78)')';-----------------------------------------------------------------------------'
    write(36,*)
    write(36,*)
    write(36,'(A16)')'[ moleculetype ]'
    write(36,'(A15)')'; Name     OCNM'
    write(36,'(A12)')'  OCNM     3'
    write(36,*)
    write(36,*)
    write(36,'(A9)')'[ atoms ]'
    write(36,'(A53)')';   nr   type  resnr residue  atom  cgnr charge  mass'
    do i=1,npart
        write(36,'(I7,A6,I5,A6,A6,I5,F10.6,F12.4)') i,atom_type(i),molecule(i),'OCNM',atom_type(i),i,q(i),mass(types(i))
    enddo

    write(36,*)
    write(36,*)

    call bonds
    write(36,'(A9)') '[ bonds ]'
    write(36,'(A36)')';  ai    aj   type    r0          kr'
    do i=1,nbond
       ia =  bond(1,i)
       ja =  bond(2,i)
       if ((atom_type(ia) == 'N_R' .and. atom_type(ja) == 'O_2') .or. &
         & (atom_type(ia) == 'O_2' .and. atom_type(ja) == 'N_R'))then
          r0 = 0.122 
          kr = 518141.439
       elseif ((atom_type(ia) == 'N_R' .and. atom_type(ja) == 'C_3') .or. & 
             & (atom_type(ia) == 'C_3' .and. atom_type(ja) == 'N_R'))then
          r0 = 0.1503 
          kr = 135777.68
       elseif ((atom_type(ia) == 'H_' .and. atom_type(ja) == 'C_3') .or. &
             & (atom_type(ia) == 'C_3' .and. atom_type(ja) == 'H_'))then
          r0 = 0.108733333333 
          kr = 328131.568
       else
          print *,'Warning: bond type not recognized!'
          print *,atom_type(ia),atom_type(ja) 
          r0 = 0.0
          kr = 0.0
       endif
       write(36,'(2I5,I5,F10.6,F14.5)') ia,ja, 1, r0, kr ! 3=morse,1=harm
    enddo

    write(36,*)
    write(36,*)

    call angles
    write(36,'(A10)') '[ angles ]'
    write(36,'(A46)')'; ai    aj   ak   type     theta0           kt'
    do i=1,n_ang
        ia = angle(1,i)
        ja = angle(2,i)
        ka = angle(3,i)
        if ((atom_type(ia) == 'C_3' .and. atom_type(ja) == 'N_R' .and. atom_type(ka)=='O_2') .or. &
          & (atom_type(ia) == 'O_2' .and. atom_type(ja) == 'N_R' .and. atom_type(ka) == 'C_3'))then
            theta0 =  116.99
            kt = 370.93135
        elseif (atom_type(ia) == 'H_' .and. atom_type(ja) == 'C_3' .and. atom_type(ka)=='H_')then
            theta0 =  111.333333333
            kt =  297.2443
        elseif ((atom_type(ia) == 'H_' .and. atom_type(ja) == 'C_3' .and. atom_type(ka)=='N_R') .or. &
          & (atom_type(ia) == 'N_R' .and. atom_type(ja) == 'C_3' .and. atom_type(ka) == 'H_'))then
            theta0 =   107.53
            kt = 480.135
        elseif (atom_type(ia) == 'O_2' .and. atom_type(ja) == 'N_R' .and. atom_type(ka)=='O_2')then
            theta0 =  126.0
            kt =   1421.237
       else
          print *,'Warning: bond type not recognized!'
          print *,atom_type(ia),atom_type(ja),atom_type(ka)
          theta0 = 0.0
          kt = 0.0
        endif

        write(36,'(3I6,I4,F12.6,F12.6)') ia,ja,ka, 1,theta0, kt ! 3=morse,1=harm
    enddo
    write(36,*)
    write(36,*)
    
    call torsions
    write(36,'(A13)') '[ dihedrals ]'
    write(36,'(A50)')';  ai    aj    ak    al   type     phi    k   mult'
    do i=1,n_dih
        ia = dihedral(1,i)
        ja = dihedral(2,i)
        ka = dihedral(3,i)
        la = dihedral(4,i)
        dihe  = torsion_angle(ia,ja,ka,la)
        !print *,dihe
        if ((dihe <360.0 .and. dihe >350.0) .or. (dihe <10.0 .and. dihe >0.0)) then
           write(36,'(4I6,I4,F11.6,F10.4,I4)') ia,ja,ka,la,1,358.10, 9.375, 2! 3=morse,1=harm
        endif
    enddo

    write(36,*)
    write(36,*)
 
    !write(36,'(A13)') '[ pairs ]'
    !write(36,'(A71)')'; name  bond_type    mass    charge   ptype          sigma      epsilon'
    !do i=1,n_mol-1
    !   do j = i+1,n_mol
    !      if ((dihe <360.0 .and. dihe >350.0) .or. (dihe <10.0 .and. dihe >0.0)) then
    !         write(36,'(4I6,I4,F11.6,F10.4,I4)') ia,ja,ka,la,1,358.10, 9.375, 2! 3=morse,1=harm
    !      endif
    !   enddo
    !enddo

    write(36,*)
    write(36,*)

    write(36,'(A10)') '[ system ]'
    write(36,'(A6)') '; Name'
    write(36,'(A5)') ' OCNM'
    write(36,*)
    write(36,*)
    write(36,'(A13)') '[ molecules ]'
    write(36,'(A22)') '; Compound       #mols'
    write(36,'(A14,I5)') '  OCNM        ',n_mol
    write(36,*)
    write(36,*)

    close(36)

    return
end

!-------------------------------------------------------------------------c
!  ------           prepare index file for gromacs           ------      -
!-------------------------------------------------------------------------c
!

subroutine write_gromacs_ndx
    use mode,only: npart,cpart,sc,a,b,c !,types,ctype, mass, q
!   use mod_nw,only: atom_type!,latmtyp  
!   use modnei,only: nbond,bond,n_ang,angle,n_dih,dihedral
    implicit none
!
!   Local Variables
!
    integer          i,j,k
    integer          nqm,nu,qme,qms   ! QM atoms in index

    character(75)    line
    character(5)     ci 


    open(37,file='index.ndx',status='unknown') !,access='append'
    write(37,'(A10)') '[ System ]'

    j = 0
    do i = 1, npart        ! writing  system groups to file
       j = j+1
       write(ci,'(I5)') i
       line((j-1)*5+1:j*5) = adjustl(ci)
       if (j == 15 .or. i==npart) then
          line = adjustl(line)
          write(37,'(A75)') line
          line = '                                                                           '
          j = 0
       endif
    enddo

    write(37,*)
    write(37,*)

    write(37,'(A8)') '[ OCNM ]'

    j = 0
    do i = 1, npart        ! writing  groups to file
       j = j+1
       write(ci,'(I5)') i
!      print *,(j-1)*5+1,' to ',j*5,' is ',ci

       line((j-1)*5+1:j*5) = adjustl(ci)
       if (j == 15 .or. i==npart) then
          line = adjustl(line)
          write(37,'(A75)') line
          line = '                                                                           '
          j = 0
       endif
    enddo

    write(37,*)
    write(37,*)
    write(37,'(A6)') '[ QM ]'

    i = int(sc(1)/2.0) 
    j = int(sc(2)/2.0) 
    k = int(sc(3)/2.0)

    nu = i*sc(2)*sc(3) + j*sc(3) + k

    j = 0
    nqm = cpart                     ! one unit cell
    qms = (nu)*nqm+1
    qme = (nu+1)*nqm
    do i = qms, qme                 ! writing  qm region atoms to file
       j = j+1
       write(ci,'(I5)') i
       line((j-1)*5+1:j*5) = adjustl(ci)
       if (j == 15 .or. i==qme) then
          line = adjustl(line)
          write(37,'(A75)') line
          line = '                                                                           '
          j = 0
       endif
    enddo
    write(37,*)
    write(37,*)
    write(*,'(A24)') '* QM region boundaries *'
    write(*,'(A24)') '  (in unit Angstrom)'
    write(*,'(3f10.5)')a(1)/sc(1),b(2)/sc(2),c(3)/sc(3)
    if(b(1)>0.1.or.b(1)<-0.01.or.c(1)>0.01.or.c(1)<-0.01.or.c(2)>0.01.or.c(2)<-0.01) then
         write(*,'(3F13.6)')b(1)/sc(1),c(1)/sc(1),c(2)/sc(2)
    endif

    close(37)

    return
end

