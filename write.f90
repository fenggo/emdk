!=========================================================================
!                             WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!

subroutine writeout
    use mode,only: x,npart !,ctype ,mass ,types
    use modctrl,only: iwr
    implicit none

    integer    i

    do i=1,npart     ! zero 
       if(abs(x(1,i))<1.0e-8) x(1,i)=0.0d0
       if(abs(x(2,i))<1.0e-8) x(2,i)=0.0d0
       if(abs(x(3,i))<1.0e-8) x(3,i)=0.0d0
    enddo

    if(iwr==0)then
          call write_df
    elseif(iwr==1)then
          call writexyz('card.xyz')
    elseif(iwr==2)then
          call writegeo
    elseif(iwr==3)then
          call writenei
    elseif(iwr==4)then
          call writenei
          call analysis
    elseif(iwr==5)then
          !print *,'in write gulp'
          call writegulp
    elseif(iwr==6)then
          call writepdb
    elseif(iwr==7)then
          call writedlpoly
    elseif(iwr==8)then
          call nw_info
    elseif(iwr==9)then
          call gofr
    elseif(iwr==10)then
          call write_poscar
    elseif(iwr==11)then
          call write_siesta
    elseif(iwr==12)then
          call write_gen
    elseif(iwr==13)then
          call write_cpmd
    elseif(iwr==18)then
          call write_cp2k
!     elseif(iwr==14)then
!           call stretch
    elseif(iwr==15)then
            !do i =1,npart
               !print *,i,mass(types(i)),ctype(types(i)),x(1,i),x(2,i),x(3,i),0,0,0
            !enddo
            call write_cfg
    elseif(iwr==16)then
            call neighbor
            call findmole
            call writeat  ! atom table
    elseif(iwr==17)then
            call write_moltop
    elseif(iwr==19)then
          call velocity
          call write_df
    elseif(iwr==21)then
          call write_nwchem('inp-nw')
    elseif(iwr==22)then
          call write_zmatrix
    elseif(iwr==23)then
          call gulp_ff
    elseif(iwr==24)then
          call tinker_ff
    elseif(iwr==25)then
          call write_top('card.top')
    elseif(iwr==26)then
          call write_uspex('uspex.xyz')
    elseif(iwr==27)then
          call mmfit
          call writeat  ! atom table
    elseif(iwr==28)then
          call neighbor
          call findmole
          call write_species
    elseif(iwr==29)then
          call write_cif
    elseif(iwr==30)then
          call write_gro
          call write_gromacs_top('molecule.top')
          call write_gromacs_ndx
    elseif(iwr==99)then
          call write_db
    elseif(iwr==100)then
          call auto_rotate
    elseif(iwr==101)then
          call auto_swing
    else
          call err('ERROR: Out put type not supported!')
    end if

   return
end

!-------------------------------------------------------------------------c
!  ------             WRITE Lammps data TO FILE              ------      -
!-------------------------------------------------------------------------c
!

    subroutine write_df
    use mode,   only: x,v,npart,types,a,b,c,r1,r2,r3,alpha,gamm,beta,nelem, alive,mass
    use modctrl,only: iwr
    implicit none
!
!   Local Variables
!
    integer         i,np,j
    real(8)         cy
    real(8)         pi, alph, bet, gam
!
    pi = acos(-1.0d0)

    call proj_coord

!  For lammps: 
    print '(2x,A)','For lammps: a in X axis, b in XY Plane, c arbitary!'

    gam = gamm*pi/180
    bet = beta*pi/180
    alph= alpha*pi/180

    a = (/r1,0.0d0,0.0d0/)
    b = (/r2*dcos(gam),r2*dsin(gam),0.0d0/)
    cy= r3*dcos(alph)-r3*dcos(bet)*dcos(gam)
    cy= cy/dsin(gam)
    c = (/r3*dcos(bet),cy,dsqrt(r3*r3*dsin(bet)*dsin(bet)-cy*cy)/)

    call cart_coord

    open(31,file='coord.df',status='unknown')
    write(31,'(A16)')'# Writed by EMDK'
    write(31,*)
    write(31,'(I7,A7)')npart, '  atoms'
    write(31,'(I7,A12)')nelem, '  atom types'
    write(31,*)
    write(31,'(2f10.5,A8)')0.0,a(1),'xlo xhi'
    write(31,'(2f10.5,A8)')0.0,b(2),'ylo yhi'
    write(31,'(2f10.5,A8)')0.0,c(3),'zlo zhi'
    write(31,*)
    if(b(1)>0.1.or.b(1)<-0.01.or.c(1)>0.01.or.c(1)<-0.01.or.c(2)>0.01.or.c(2)<-0.01) then
    	  write(31,'(3F13.6,A9)')b(1),c(1),c(2),' xy xz yz'
    endif
    !write(31,*)a(1),a(2),a(3)
    !write(31,*)b(1),b(2),b(3)
    !write(31,*)c(1),c(2),c(3)

    write(31,*)
    write(31,'(A6)')'Masses'
    write(31,*)
    do i = 1, nelem
       write(31,*)i,mass(i)
    enddo
    write(31,*)
    write(31,'(A5)')'Atoms'
    write(31,*)

    np = 0
    do i=1,npart
         if(alive(i)/=0) then
              np = np + 1
              write(31,500)np,types(i),0.0,x(1,i),x(2,i),x(3,i)
         endif
    end do
!
!   shock velocity
      if(iwr==19)then  
         write(31,*)
         write(31,'(A10)')'Velocities'
         write(31,*)

         do i=1,npart
             write(31,*)i,(v(j,i),j=1,3)
         enddo
      endif

500 FORMAT(I6,I3,F10.5,3F11.5)

   close(31)
!
   return
   end
!
!-------------------------------------------------------------------------
!                Write DL_PLOY format coordinate(CONFIG)                 -
!-------------------------------------------------------------------------
!
    subroutine writedlpoly
    use mode,only: x,npart,types,ctype,a,b,c
    implicit none
!
    integer i
!
    print 100,'Write our coordinate for DL_POLY ...'

    open(31,file='CONFIG.x',status='unknown')
    write(31,100)'Gnerated by symmetry Program! Writed by G.F.(email:gfeng.alan@gmail.com)'
    write(31,300)2,2,npart
    write(31,500)a(1),0.0,0.0
    write(31,500)0.0,b(2),0.0
    write(31,500)0.0,0.0,c(3)

    do i=1,npart
       write(31,400)ctype(types(i)),i
       write(31,500)x(1,i)-0.5*a(1),x(2,i)-0.5*b(2),x(3,i)-0.5*c(3)
       write(31,500)0.0,0.0,0.0
       write(31,500)0.0,0.0,0.0
    end do
!
100 format(2X,A)
300 format(2I4,I6)
400 format(A4,I6)
500 FORMAT(3F11.5)

    close(31)
!
    return
    end
!
!-------------------------------------------------------------------------
!  ------             WRITE DATA TO FILE   XYZ FORMAT      ------        -
!-------------------------------------------------------------------------
!
subroutine writexyz(filename)
    use mode,only: x,npart,types,ctype,r1,r2,r3,alpha,beta,gamm,alive
    implicit none
!
!   Local Variables
!
    integer          i
    !integer          icoor
    !real(8)          icub(3)
    !real(8)          bx,ay
    real(8)          pi
    character(*)    filename

    pi = dacos(-1.0d0)
!
    open(30,file=filename,status='unknown',access='append')

!     print 100,'Choosing Coordinate type for Output:'
!     print 100,'1. Fractional coordinate.'
!     print 100,'2. Real Coordinate.'
!     print 100,'In put 1 or 2:'
!     read(*,*)icoor

!  c along Z, a in XZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      a = (/r1*dsin(beta),0.0d0,r1*dcos(beta)/)
!      bx= r2*dcos(gamm)-r2*dcos(alpha)*dcos(beta)
!      bx= bx/dsin(beta)
!      b = (/bx,dsqrt(r2*r2*dsin(alpha)*dsin(alpha)-bx*bx),r2*dcos(gamm)/)

!  c along Z, b in YZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      b = (/0.0d0,r2*dsin(alpha),r2*dcos(alpha)/) 
!      ay= r1*dcos(gamm)-r1*dcos(alpha)*dcos(beta)
!      ay= ay/dsin(alpha)
!      a = (/dsqrt(r1*r1*dsin(beta)*dsin(beta)-ay*ay),ay,r1*dcos(beta)/)

    write(30,*)npart
    !print *,npart
    write(30,'(A4,3F12.6,3F12.6,A5)')'PBC: ',r1,r2,r3, &
              alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    !print *,'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    do i=1,npart
          if(alive(i)/=0) write(30,*)ctype(types(i)),x(1,i),x(2,i),x(3,i)
    end do
!    print 100,'Coordinates are saved in:',filename


    close(30)

    return
end


!-------------------------------------------------------------------------
!  ------             WRITE DATA TO FILE   GEO  FORMAT      -----        -
!-------------------------------------------------------------------------
!
subroutine writegeo
     use mode, only: x,npart,types,ctype,a,b,c,alpha,beta,gamm,r1,r2,r3
     implicit none
!
     integer               ::  ibgfversion =200
     integer                   i,j
     character(4)              cty
     real(8) ::                pi = 3.1415944
     real(8) ::                ay !bx,

     pi = dacos(-1.0d0)
 
     open(32,file='card.geo',status='unknown')

     write(32,'(A6,I4)')'XTLGRF',ibgfversion
     write(32,'(A6,A60)')'DESCRP','Information about coord'
     write(32,500)'CRYSTX',r1,r2,r3,alpha*180.0d0/pi,beta*180.0d0/pi,gamm*180.0d0/pi
     write(32,510)

     call proj_coord

!  c along Z, a in XZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      a = (/r1*dsin(beta),0.0d0,r1*dcos(beta)/)
!      bx= r2*dcos(gamm)-r2*dcos(alpha)*dcos(beta)
!      bx= bx/dsin(beta)
!      b = (/bx,dsqrt(r2*r2*dsin(alpha)*dsin(alpha)-bx*bx),r2*dcos(gamm)/)

!  For Reax geo
!     c along Z, b in YZ plane
     c = (/0.0d0,0.0d0,r3/)
     b = (/0.0d0,r2*dsin(alpha),r2*dcos(alpha)/) 
     ay= r1*dcos(gamm)-r1*dcos(alpha)*dcos(beta)
     ay= ay/dsin(alpha)
     a = (/dsqrt(r1*r1*dsin(beta)*dsin(beta)-ay*ay),ay,r1*dcos(beta)/)

     call cart_coord

     do i=1,npart
        cty=ctype(types(i))
        write(32,520)i,cty,(x(j,i),j=1,3),cty,types(i),1,0.00
     end do
     write(32,530)
     write(32,540)
!
 500 format(A6,6f11.5)
 510 format ('FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,',&
              '3f10.5,1x,a5,i3,i2,1x,f8.5)')
 520 format ('HETATM',1x,i5,1x,a2,3x,1x,3x,1x,1x,1x,5x,3f10.5,1x,&
             a5,i3,i2,1x,f8.5)
 530 format('FORMAT CONECT (a6,12i6)')
 540 format('END')

     close(32)
!
     return
end

!-------------------------------------------------------------------------
!             Write out neighbors                                        !
!-------------------------------------------------------------------------

   subroutine writenei
   use modnei,only: atable,matom
   use mode,  only: npart
   implicit none

   integer       i
   integer       j
   integer       ic

   open(33,file='NEICAR',status='unknown')
   write(33,100)'Neighbor of every atom'

   ic = matom +1

   do i=1,npart
      write(33,200)i,(atable(j,i),j=1,ic)
   end do

100   format(2x,A)
200   format(2I5,20I5)

   close(33)
   return
   end

!
!-------------------------------------------------------------------------
!        Write out Atom-table                                            -
!-------------------------------------------------------------------------
!
     subroutine writeat
     use mode, only: npart,ctype,types
     use modnei,only: matom, atable
     implicit none
     integer i,j,ns
     character(80) line
     character(5)  nm 
!
     open(14,file='atable.log',status='unknown')
!     write(14,*)'Atom Bond Table'
     write(14,*)'ATOM TABLE'
!
     do i=1,npart
        line = ' '
        do j =2,atable(1,i)+1
           if (atable(j,i) /= 0) then
              write(nm,'(I5)') atable(j,i)
              !print *,nm
              line(1:3) = ctype(types(i))
              ns = len_trim(line)
              line(ns+1:ns+5)= nm
              line(ns+6:ns+6)= ' '
              ns=ns+6
              line(ns+1:ns+3) = ctype(types(atable(j,i)))
              !print *,line
           endif
        enddo
        !write(14,500)i,(atable(j,i),j=1,8)!matom)
        write(14,'(I5,A1,A)')i,' ',line
     end do
!
     close(14)
!
     return
     end

!-------------------------------------------------------------------------
!     Write out atom coordinates in PDB format                           !
!-------------------------------------------------------------------------

    subroutine writepdb
    use mode,only: r1,r2,r3,alpha,beta,gamm,x,npart,types,ctype,q,alive
    use modfim,only: molecule
    implicit none
!
!   Local Variables
!
    integer        i,l
    integer        n_e(4)
    !character(3)   res_c(4)
    character(3):: atm_tp(4) = (/'C_3','H_ ','O_2','N_R'/)
    character(12)  str(4)
    real(8) ::     pi = 3.1415944
!
!    do i=1,4
!       print *,'Input residue name for atom type',i,':'
!       read(*,*)res_c(i)
!    end do
     open(31,file='card.pdb',status='unknown')

! 列  	数据类型  	字段名称  	定义描述
! 1 - 6  	Record name  	"CRYST1" 	
! 7 - 15  	Real(9.3)  	a 	a(埃)
! 16 - 24  	Real(9.3)  	b 	b(埃)
! 25 - 33  	Real(9.3)  	c 	c(埃)
! 34 - 40 	Real(7.2) 	alpha  	α(度)
! 41 - 44  	Real(7.2) 	beta  	β(度)
! 48 - 54  	Real(7.2) 	gamm  	γ(度)
! 56 - 66  	LString  	sGroup  	空间基团
! 67 - 70  	Integer 	z 	Z 值
     write(31,400)'CRYST1',r1,r2,r3,alpha*180.0d0/pi,beta*180.0d0/pi,gamm*180.0d0/pi,'P1'
! PDB ATOM 记录格式
! 列  	数据类型  	字段名称  	定义描述
! 1 - 6  	Record name  	"ATOM " 	
! 7 - 11   	Integer  	serial 	原子序列号
! 13 - 16  	Atom  	name  	原子名
! 17 	Character 	altLoc  	交替位点标识
! 18 - 20  	Residue name  	resName  	残基名
! 22 	Character 	chainID  	链标识符
! 23 - 26  	Integer  	resSeq  	残基序列号
! 44 	AChar 	iCode  	残基的插入代码
! 31 - 38  	Real(8.3) 	x 	X坐标
! 39 - 46  	Real(8.3) 	y 	Y坐标
! 44 - 54  	Real(8.3) 	z 	Z坐标
! 55 - 60  	Real(6.2)  	occupancy 	空间大小
! 61 - 66  	Real(6.2)  	tempFactor 	温度因数
! 73 - 76  	LString(4)  	segID  	段标识符，左对齐
! 77 - 78  	LString(2)  	element 	元素符号，右对齐
! 79 - 80  	LString(2)  	charge  	原子的价位
     n_e(:) = 0
     do i=1,npart
        if(alive(i)/=0) then
        n_e(types(i))=n_e(types(i))+1
        write(str(types(i)),*)n_e(types(i))
        str(types(i)) = adjustl(str(types(i)))
!        write(31,500)'ATOM  ',i,ctype(types(i)),str(types(i)),ctype(types(i)),1,x(1,i),x(2,i),x(3,i),&
!                     1.00,0.00,ctype(types(i))
        !write(31,'(A6,I5,1X,A4,I4,2X,I4,4X,3F8.3,2F6.2,10X,A2)')'ATOM  ',molecule(i),ctype(types(i)),&
                     !i,1,(x(l,i),l=1,3),1.00,q(i),ctype(types(i)) ! Oringle!!!
        write(31,'(A6,I5,1X,A4,A4,2X,I4,4X,3F8.3,2F6.2,10X,A2)')'ATOM  ',i,atm_tp(types(i)),&
                     'CNM',molecule(i),(x(l,i),l=1,3),1.00,q(i),ctype(types(i)) ! HMX-CL20
!        write(31,500)'ATOM  ',i,'CJ',res_c(types(i)),1,x(1,i),x(2,i),x(3,i),&
!                     1.00,0.00,ctype(types(i))
        endif
     end do
!
400 FORMAT(A6,3F9.3,3F7.2,A11,4X)


    close(31)
!
    return
    end
!
!-------------------------------------------------------------------------
!  ------             Write coordinate of gulp input format              -
!-------------------------------------------------------------------------
!

    subroutine writegulp
    use mode,only: xd,npart,types,ctype,q, a, b, c,alive !, nelem,x,
    !use mod_grd, only: fix
    use modctrl,only: run_type
    implicit none
!
!   Local Variables
!
    integer i
!
    call proj_coord

    open(31,file='inp-gulp',status='unknown')
     !print *,'cannot open here ?'
    if(run_type==0) then
         write(31,'(A9)')'opti conp'
         write(31,*) !opti conv nosymmetry molecule
         write(31,'(A11)')'maxcyc 2000'
         !write(31,'(A9)')'# for MD:'
    elseif(run_type==1) then
         write(31, '(A7)') 'md conv'
         write(31, '(A12)') 'ensemble nvt'

         write(31, '(A23)') 'tau_thermostat  0.05 ps'
         write(31, '(A24)') 'temperature 298.000000 K'

         write(31, '(A20)') 'timestep    0.001 ps'
         write(31, '(A20)') 'production     20 ps'
         write(31, '(A20)') 'equilibration   0 ps'
         write(31, '(A17)') 'write          50'
         write(31, '(A17)') 'sample         50'
    elseif(run_type==2) then

         write(31, '(A7)') 'md conp'
         write(31, '(A26)') 'integrator leapfrog verlet'
         write(31, '(A24)') 'ensemble npt 0.005 0.005'
         write(31, '(A17)') 'temperature 300 K'
         write(31, '(A17)') 'pressure 0.00 GPa'
         write(31, '(A22)') 'tau_barostat    0.1 ps'
         write(31, '(A22)') 'tau_thermostat  0.1 ps'

         write(31, '(A1)')'#'
         write(31, '(A1)')'#'


         write(31, '(A20)') 'timestep        1 fs'
         write(31, '(A20)') 'production     20 ps'
         write(31, '(A20)') 'equilibration 0.5 ps'
         write(31, '(A17)') 'write          50'
         write(31, '(A17)') 'sample         50'
    elseif(run_type==0) then
         write(31,'(A42)')'gradient nosymmetry conp molmec #conjugate'
    endif

    write(31,*)
    write(31,'(A5)') 'title'
    write(31,'(A16)') 'GULP calculation'
    write(31,'(A3)') 'end'
    write(31,*)

    !write(31,'(A10)') 'spacegroup'
    !write(31,'(A4)') ' P 1'
!library
    write(31,'(A11)') 'library bop'
    write(31,*)
    write(31,'(A24)') 'output movie xyz his.xyz'
    write(31,'(A19)') 'dump 50 restart.grs'
    write(31,*)
    write(31,'(A7)') 'vectors'
    write(31,'(3F16.8)') (a(i),i=1,3)
    write(31,'(3F16.8)') (b(i),i=1,3)
    write(31,'(3F16.8)') (c(i),i=1,3)
    write(31,*)
    write(31,'(A10)') 'fractional'
     !print *,'problem here ...'
     !stop
    do i=1,npart
         if(alive(i)/=0) then
            write(31,600)ctype(types(i)),'core',xd(1,i),xd(2,i),xd(3,i),q(i),1.0,0.0,1,1,1 ! 1=vary 0=fix
         endif
    end do
!
     !write(31,'(A7)')'Species'
     !do i=1,nelem
       ! write(31,'(A2,A6,A4)')ctype(i),' core ',ctype(i)
     !enddo 
    write(31,*)
     !write(31,'(A10)') 'spacegroup'
     !write(31,'(A3)') 'P 1'
!

 600 FORMAT(A3,A5,3F15.6,3F10.4,3I3)

    close(31)
!
    return
    end
!-------------------------------------------------------------------------
!  ------        Write coordinate of nwchem input format                 -
!-------------------------------------------------------------------------

subroutine write_nwchem(filename)
    use mode,only: x,npart,types,ctype,alive !,a,b,c ,r1,r2,r3
    use mod_nw
    implicit none
!
!   Local Variables
!
    integer          i
    integer          mult,mm !,icoor
    !real(8)          icub(3)
    !real(8)          bx,ay
    real(8)          pi
    character(*)     filename
    !character(*)     basis

    pi = dacos(-1.0d0)
!
    !print *,filename
    open(30,file=filename,status='unknown')

    !write(30,*)npart
    !write(30,'(A4,3F12.6,3F12.6,A5)')'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    !print *,'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    write(30,'(A14)') 'geometry autoz'
    mult = 0
    do i=1,npart
          if(alive(i)/=0) then
             write(30,*)ctype(types(i)),x(1,i),x(2,i),x(3,i)
             !print *,ctype(types(i))
             !print *,types(i)
             !print *,ctype
             if (ctype(types(i))=='C') mult = mult + 6
             if (ctype(types(i))=='H') mult = mult + 1
             if (ctype(types(i))=='N') mult = mult + 7
             if (ctype(types(i))=='O') mult = mult + 8
          endif
    end do
    mm = mod(mult,2)
    if (mm == 0) then
        mult = 1
    else 
        mult = 2
    endif
    write(30,'(A3)') 'end'
    write(30,*) 
    write(30,*) 
    write(30,'(A5)') 'basis'
    write(30,'(A12,A10)') '  * library ', basis
    write(30,'(A3)') 'end'
    write(30,*) 
    write(30,*) 
    write(30,'(A3)')  'dft'
    write(30,'(A6,I3)')  '  mult',mult
    write(30,'(A10)') '  xc b3lyp'
    write(30,'(A26)') '  convergence energy 1e-07'
    write(30,'(A3)')  'end'
    write(30,*) 
    write(30,*) 
    !write(30,'(A15)') 'task dft energy'

    close(30)

    return
end

!-------------------------------------------------------------------------
!  ------       Write coordinate of gaussian input format                -
!-------------------------------------------------------------------------
subroutine write_gaussian(filename)
    use mode,only: x,npart,types,ctype,alive !,alpha,beta,gamm ,r1,r2,r3
    implicit none
!
!   Local Variables
!
    integer          i,k
    !integer          icoor
    !real(8)          icub(3)
    real(8)          pi
    character(*)     filename
    character(20)    inf
    character(40)    line

    pi = dacos(-1.0d0)
!
    open(33,file=filename,status='unknown')

    !write(30,*)npart
    !write(30,'(A4,3F12.6,3F12.6,A5)')'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    !print *,'PBC: ',r1,r2,r3, &
              !alpha*180.0/pi,gamm*180.0/pi,beta*180.0/pi,'(P1)'
    k = len(trim(filename))
    inf(1:k-4) = filename(1:k-4) 
    !inf()
    !print *,inf
    line = '                                        '
    line(1:5) = '%chk='
    line(6:k+1) = inf
    line(k+2:k+5) = '.chk'
    line = trim(line)
    write(33,'(1x,A40)') line
    !write(*,*) '%chk=',inf,'.chk'
    write(33,'(A16)') '# b3lyp/6-311++g'
    write(33,*) 
    write(33,'(A10)') 'hmx-energy'
    write(33,*) 
    write(33,'(A3)') '0 1'
    do i=1,npart
          if(alive(i)/=0) write(33,*)ctype(types(i)),x(1,i),x(2,i),x(3,i)
    end do

    write(33,*) 

    close(33)

    return
end
!
!-------------------------------------------------------------------------
!         Write zmatrix for CALYPSO Solid structurePrediction            -
!-------------------------------------------------------------------------
!
subroutine write_zmatrix
    !use modfim,only: n_mol
    use modnei,only: nbond,bond,n_ang,n_dih,angle,dihedral,mang,mdih
    use mode,  only: x,npart !, a, b, c
    implicit none

    integer         i,j!,k
    !integer         n_atom
    !integer,allocatable:: ang_numb(:)
    !real(8)      :: pi=3.1415926535897932
    !real(8)      :: mor
    !real(8)      :: hrm
    !real(8)      :: r0=1.4000
    !real(8)      :: kj2ev=1.03643E-02
    !real(8)         kc
    !real(8)         kcth
    !real(8)         thet
    !real(8)         kdih
    !real(8)         dih
    real(8)         r1(3),r2(3),r3(3),ra(3),rb(3)
    !real(8)         r_a,r_b,r_c !,ina,inb,inc
    real(8)         rsq,radi,rsq1,rsq2
    real(8)         delx,dely,delz
    real(8)         ang,cos_ang,dihe !,ave_ang
    logical,allocatable:: bond_no_dih(:,:)
    real(8),allocatable:: bond_length(:,:)
    real(8),allocatable:: ang_value(:)

    call  neighbor
    call  findmole
!      call writeat
!  Determining Angles

    n_ang = 0
    !print *,nbond

    allocate(bond_no_dih(nbond,nbond))
    allocate(bond_length(nbond,nbond))
    allocate(ang_value(3*npart))

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

    !print *,n_ang
!  Determining dihedrals
    bond_no_dih(:,:) = .true.

    n_dih = 0
    do i=1,n_ang-1
       do j=i+1,n_ang  
          if(n_dih>=mdih)call err('Error: Too many dihedrals, Modify mdih and retry!')
          if(angle(2,i)==angle(1,j).and.angle(2,j)==angle(1,i))then
            if(bond_no_dih(angle(1,i),angle(2,i)))then
                bond_no_dih(angle(1,i),angle(2,i)) = .false.
                bond_no_dih(angle(2,i),angle(1,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(3,j)
                dihedral(2,n_dih)=angle(1,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(3,i)
            endif
          elseif(angle(2,i)==angle(3,j).and.angle(2,j)==angle(1,i))then
            if(bond_no_dih(angle(1,i),angle(2,i)))then
                bond_no_dih(angle(1,i),angle(2,i)) = .false.
                bond_no_dih(angle(2,i),angle(1,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(1,j)
                dihedral(2,n_dih)=angle(1,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(3,i)
            endif
          elseif(angle(2,i)==angle(1,j).and.angle(2,j)==angle(3,i))then
            if(bond_no_dih(angle(2,i),angle(3,i)))then
                bond_no_dih(angle(2,i),angle(3,i)) = .false.
                bond_no_dih(angle(3,i),angle(2,i)) = .false.

                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(3,j)
                dihedral(2,n_dih)=angle(3,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(1,i)
            endif
          elseif(angle(2,i)==angle(3,j).and.angle(2,j)==angle(3,i))then
            if(bond_no_dih(angle(2,i),angle(3,i)))then
                bond_no_dih(angle(2,i),angle(3,i)) = .false.
                bond_no_dih(angle(3,i),angle(2,i)) = .false.
                n_dih=n_dih+1
                dihedral(1,n_dih)=angle(1,j)
                dihedral(2,n_dih)=angle(3,i)
                dihedral(3,n_dih)=angle(2,i)
                dihedral(4,n_dih)=angle(1,i)
            endif
          endif
       enddo
    enddo

    open(45,file='zmatrix.txt',status='unknown')

    write(45,'(I3)')nbond+n_ang+n_dih+1   ! the number of varibles
    write(45,'(10I3)') bond(1,1),0,0,0,0,0,0,0,0,0
    do i=1,nbond
       delx=x(1,bond(1,i))-x(1,bond(2,i))
       dely=x(2,bond(1,i))-x(2,bond(2,i))
       delz=x(3,bond(1,i))-x(3,bond(2,i))
       rsq = delx*delx + dely*dely + delz*delz
       radi = sqrt(rsq)
       bond_length(bond(1,i),bond(2,i)) = radi 
       bond_length(bond(2,i),bond(1,i)) = radi 
       write(45,'(2I5,2I5,F11.8,3I5)') bond(1,i),bond(2,i),0,0,radi, 0,0,0
    enddo
   
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
       ang_value(angle(1,i)+angle(2,i)+angle(3,i)) = ang
       write(45,'(4I5,3F12.8,3I5)') angle(1,i),angle(2,i),angle(3,i),0,bond_length(angle(1,i),angle(2,i)),ang,0.0,0,0,0
    enddo

    ! still have problems
    print '(A56)','** The determing of dihedral angle still have problems, '
    print '(A29)','** you need set it manuly ...'
    do i=1,n_dih
       r1(1)=x(1,dihedral(2,i))-x(1,dihedral(1,i))
       r1(2)=x(2,dihedral(2,i))-x(2,dihedral(1,i))
       r1(3)=x(3,dihedral(2,i))-x(3,dihedral(1,i))

       r2(1)=x(1,dihedral(3,i))-x(1,dihedral(2,i))
       r2(2)=x(2,dihedral(3,i))-x(2,dihedral(2,i))
       r2(3)=x(3,dihedral(3,i))-x(3,dihedral(2,i))

       r3(1)=x(1,dihedral(4,i))-x(1,dihedral(3,i))
       r3(2)=x(2,dihedral(4,i))-x(2,dihedral(3,i))
       r3(3)=x(3,dihedral(4,i))-x(3,dihedral(3,i))

       ! cross product of vectors
       ra(1) = r1(2)*r2(3) - r1(3)*r2(2)
       ra(2) = r1(3)*r2(1) - r1(1)*r2(3)
       ra(3) = r1(1)*r2(2) - r1(2)*r2(1)

       rb(1) = r2(2)*r3(3) - r2(3)*r3(2)
       rb(2) = r2(3)*r3(1) - r2(1)*r3(3)
       rb(3) = r2(1)*r3(2) - r2(2)*r3(1)

       !dot  product of vectors
       rsq1= ra(1)*ra(1)+ra(2)*ra(2)+ra(3)*ra(3)
       rsq2= rb(1)*rb(1)+rb(2)*rb(2)+rb(3)*rb(3)
       radi= sqrt(rsq1*rsq2)

       cos_ang=ra(1)*rb(1) + ra(2)*rb(2) + ra(3)*rb(3)
       !print *,cos_ang
       cos_ang=-cos_ang/radi
       !print *,cos_ang
       dihe    = acos(cos_ang)
       !print *,dihe

       write(45,'(4I5,3F12.8,3I5)')dihedral(1,i),dihedral(2,i),dihedral(3,i),dihedral(4,i), &
                               bond_length(dihedral(1,i),dihedral(2,i)), &
                              ang_value(dihedral(1,i)+dihedral(2,i)+dihedral(3,i)),&
                              dihe,0,0,0
    enddo

    close(45)

    return
end
!
!
!-------------------------------------------------------------------------
!  ------            WRITE DATA TO FILE POSCAR FORMAT      ------        -
!-------------------------------------------------------------------------
!
subroutine write_poscar
    use mode,only: x,npart,types,ctype,nelem
    implicit none
!
!   Local Variables
!
    integer i
    integer j
    !integer icoor
    integer elem_n(20)
    real(8) icub(3)
!
    open(30,file='coord.POSCAR',status='unknown')

    !icub = cube

    write(30,100)'You_name_it'
    write(30,*)1.0
    !write(30,*)cube(1),0.0,0.0
    !write(30,*)0.0,cube(2),0.0
    !write(30,*)0.0,0.0,cube(3)

    elem_n= 0

    do i=1,npart
       do j=1, nelem
          if(ctype(types(i))==ctype(j))elem_n(j)=elem_n(j)+1
       enddo
    enddo

    write(30,*)(ctype(i),i=1,nelem),(elem_n(j),j=1,nelem)
    write(30,100)'Cartesian'

    do j=1,nelem
       do i=1,npart
          if(types(i)==j)write(30,*)x(1,i)/icub(1),x(2,i)/icub(2),x(3,i)/icub(3)
       enddo
    enddo

    print 100,'Coordinates are saved in coord.POSCAR!'

100 format(2X,A)
    close(30)

    return
end 
!

!-------------------------------------------------------------------------
!  ------             WRITE Coordinate for Siesta Input      ------      -
!-------------------------------------------------------------------------
!
    subroutine write_siesta
    use mode,only: x,npart,types,ctype,a,b,c,r1,r2,r3,nelem
    implicit none
!
!   Local Variables
!
    integer i
    !real(8)        bx,ay
    real(8)        comp(3)
    real(8)        minr
!
    open(30,file='coord.fdf',status='unknown')

    print *,'Warning: element keyword should be present in in.txt to write siesta inpput!'

!  c along Z, a in XZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      a = (/r1*dsin(beta),0.0d0,r1*dcos(beta)/)
!      bx= r2*dcos(gamm)-r2*dcos(alpha)*dcos(beta)
!      bx= bx/dsin(beta)
!      b = (/bx,dsqrt(r2*r2*dsin(alpha)*dsin(alpha)-bx*bx),r2*dcos(gamm)/)

!  c along Z, b in YZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      b = (/0.0d0,r2*dsin(alpha),r2*dcos(alpha)/) 
!      ay= r1*dcos(gamm)-r1*dcos(alpha)*dcos(beta)
!      ay= ay/dsin(alpha)
!      a = (/dsqrt(r1*r1*dsin(beta)*dsin(beta)-ay*ay),ay,r1*dcos(beta)/)

    comp(1)=r1
    comp(2)=r2
    comp(3)=r3
    
    minr= r1
    do i=1,3
       if(minr>=(comp(i))) minr= comp(i)  
    enddo

    write(30,*)'# Species and atoms'
    write(30,*)'NumberOfSpecies',nelem
    write(30,*)'NumberOfAtoms',npart
    write(30,*)
    write(30,*)'%block Chemical_Species_Label'
    do i=1,nelem
       write(30,*)i,ctype(i)
    enddo
    write(30,*)'%endblock Chemical_Species_Label'
    write(30,*)
    write(30,*)
    write(30,*)'LatticeConstant       1.00 Ang'  
    write(30,*)
    write(30,*)'%block LatticeVectors'
    write(30,'(3F16.8)') (a(i),i=1,3)
    write(30,'(3F16.8)') (b(i),i=1,3)
    write(30,'(3F16.8)') (c(i),i=1,3)
    write(30,*)'%endblock LatticeVectors'
    write(30,*)
    write(30,*)
    write(30,'(A12,F16.8,A4)')'KgridCutoff',minr,'Ang'
    write(30,*)
    write(30,*)
    write(30,*)'AtomicCoordinatesFormat Ang'
    write(30,*)
    write(30,*)'%block AtomicCoordinatesAndAtomicSpecies'

    do i=1,npart
          write(30,*)x(1,i),x(2,i),x(3,i),types(i)
    enddo

    write(30,*)'%endblock AtomicCoordinatesAndAtomicSpecies'
    write(30,*)
    print 100,'Coordinates are saved in coord.fdf!'

 100 format(2X,A)

    close(30)

    return
    end
!
!--------------------------------dftb gen format---------------------------------
!

    subroutine write_gen
    use mode,only: a,b,c,x,npart,types,nelem,ctype
    implicit none

!   Local Variables

    integer        i,j
    real(8)        minx,miny,minz

     open(31,file='card.gen',status='unknown')

     write(31,*)npart,'S'
     write(31,*)(ctype(j),j=1,nelem)
     do i=1,npart
        write(31,500)i,types(i),x(1,i),x(2,i),x(3,i)
        if(i>1)then
            if(minx>x(1,i))minx=x(1,i)
            if(miny>x(2,i))miny=x(2,i)
            if(minz>x(3,i))minz=x(3,i)
        else
            minx=x(1,i)
            miny=x(2,i)
            minz=x(3,i)
        endif
     enddo

     write(31,*)minx,miny,minz
     write(31,*)(a(j),j=1,3)
     write(31,*)(b(j),j=1,3)
     write(31,*)(c(j),j=1,3)

 500 FORMAT(I6,I3,3F11.5)

    close(31)

    return
    end
!
!
!-------------------------------------------------------------------------
!  ------             WRITE DATA TO FILE   CPMD input FORMAT      ------ -
!-------------------------------------------------------------------------
!

subroutine write_cpmd
    use mode,only: x,npart,types,ctype,a,b,c,r1,r2,r3,alpha,beta,gamm,&
                   nelem
    use modctrl,only: run_type
    implicit none
!
!   Local Variables
!
    integer         i
    integer         j
    integer         nt(6)
    real(8)         cy !bx,ay,
    real(8)         pi,d2r,gamm_r,alpha_r,beta_r
    
    pi = dacos(-1.0d0)
!
    open(30,file='coord.cpmd',status='unknown')


!  c along Z, a in XZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      a = (/r1*dsin(beta),0.0d0,r1*dcos(beta)/)
!      bx= r2*dcos(gamm)-r2*dcos(alpha)*dcos(beta)
!      bx= bx/dsin(beta)
!      b = (/bx,dsqrt(r2*r2*dsin(alpha)*dsin(alpha)-bx*bx),r2*dcos(gamm)/)

!  c along Z, b in YZ plane
!      c = (/0.0d0,0.0d0,r3/)
!      b = (/0.0d0,r2*dsin(alpha),r2*dcos(alpha)/) 
!      ay= r1*dcos(gamm)-r1*dcos(alpha)*dcos(beta)
!      ay= ay/dsin(alpha)
!      a = (/dsqrt(r1*r1*dsin(beta)*dsin(beta)-ay*ay),ay,r1*dcos(beta)/)

   d2r = pi/180.0
   gamm_r = gamm*d2r     ! Degree to Rad
   alpha_r = alpha*d2r
   beta_r = beta*d2r

   call proj_coord
   a = (/r1,0.0d0,0.0d0/)
   b = (/r2*dcos(gamm_r),r2*dsin(gamm_r),0.0d0/)
   !print *,gamm,beta
   cy= r3*dcos(alpha_r)-r3*dcos(beta_r)*dcos(gamm_r)
   cy= cy/dsin(gamm_r)
   c = (/r3*dcos(beta_r),cy,dsqrt(r3*r3*dsin(beta_r)*dsin(beta_r)-cy*cy)/)
   call cart_coord

   write(30,'(2X,A43)')'!###RESTART WAVEFUNCTION COORDINATES LATEST'
   write(30,'(2X,A84)')'!###RESTART WAVEFUNCTION COORDINATES VELOCITIES CELL ACCUMULATORS NOSEE NOSEP LATEST'
   write(30,'(A5)')'&CPMD'

   if(run_type==0) then
      write(30,'(2X,A61)')'OPTIMIZE WAVEFUNCTION ! GEOMETRY XYZ mult task can be asigned' 
      write(30,'(2X,A21)')'STEEPEST DESCENT CELL'
      write(30,'(2X,A23)')'PCG MINIMIZE ! or ODIIS'
      write(30,'(2X,A5)')'EMASS'
      write(30,'(4X,A3)')'250'
      write(30,'(2X,A7)')'MAXSTEP'
      write(30,'(4X,A3)')'200'
      write(30,'(2X,A10)')'MAXCPUTIME'
      write(30,'(4X,A5)')'20000'
      write(30,'(2X,A20)')'CONVERGENCE ORBITALS'
      write(30,'(4X,A6)')'1.0d-7'
   elseif(run_type>=1) then       !! NVE 
      write(30,'(2X,A21)')'MOLECULAR DYNAMICS CP'
      write(30,'(2X,A39)')'RESTART WAVEFUNCTION COORDINATES LATEST'
      write(30,'(2X,A5)')'EMASS'
      write(30,'(4X,A3)')'270'
      if (run_type==2) then       !! NVT
         write(30,'(2X,A11)')'TEMPERATURE'
         write(30,'(4X,A3)')'300'
         write(30,'(2X,A9)')'NOSE IONS'
         write(30,'(4X,A8)')'300  10'
         write(30,'(2X,A14)')'NOSE ELECTRONS'
         write(30,'(4X,A11)')'0.002 10000'
      endif
      write(30,'(2X,A8)')'TIMESTEP'
      write(30,'(4X,A21)')'2.5      #### 0.06 fs'
      write(30,'(2X,A7)')'MAXSTEP'
      write(30,'(4X,A4)')'1000'
      write(30,'(2X,A17)')'VDW CORRECTION ON'
      write(30,'(2X,A10)')'MAXCPUTIME'
      write(30,'(4X,A5)')'20000'
      write(30,'(2X,A20)')'CONVERGENCE ORBITALS'
      write(30,'(4X,A6)')'1.0d-8'
      write(30,'(2X,A13)')'STRESS TENSOR'
      write(30,'(4X,A2)')'10'
      write(30,'(2X,A5)')'STORE'
      write(30,'(4X,A3)')'100'
   endif

   write(30,'(A4)')'&END'
   write(30,*)

   write(30,'(A4)')'&DFT'
   write(30,'(2X,A22)')'FUNCTIONAL PBE  # BLYP'
   write(30,'(2X,A9)')'GC-CUTOFF'
   write(30,'(2X,A6)')'1.0D-8'
   write(30,'(A4)')'&END'
   write(30,*)

   write(30,'(A4)')'&VDW'
   write(30,'(2X,A14)')'VDW CORRECTION'
   write(30,'(2X,A10)')'ALL DFT-D2'
   write(30,'(A4)')'&END'
   write(30,*)

   write(30,'(A7)')'&SYSTEM'
   write(30,'(2X,A8)')'ANGSTROM'
   write(30,'(2X,A4,1X,A7)')'CELL','VECTORS'
   write(30,'(3F15.6)')(a(i),i=1,3)
   write(30,'(3F15.6)')(b(i),i=1,3)
   write(30,'(3F15.6)')(c(i),i=1,3)
   write(30,'(2X,A6)')'CUTOFF'
   write(30,'(3X,A4)') '75.'
   write(30,'(A4)')'&END'
   write(30,*)

   nt = 0


   do i=1,npart
      if(types(i)==1)nt(1)=nt(1)+1
      if(types(i)==2)nt(2)=nt(2)+1
      if(types(i)==3)nt(3)=nt(3)+1
      if(types(i)==4)nt(4)=nt(4)+1
      if(types(i)==5)nt(5)=nt(5)+1
      if(types(i)==6)nt(6)=nt(6)+1
   enddo


   write(30,'(A6)')'&ATOMS'
   do j=1,nelem
      write(30,'(A1,A2,A17)')'*',adjustl(ctype(j)),' KLEINMAN-BYLANDE'
      write(30,'(A12)')'LMAX=P LOC=P'
      write(30,'(I5)')nt(j)
      do i=1,npart
         if(types(i)==j)write(30,'(3F20.8)')x(1,i),x(2,i),x(3,i)
      enddo
   enddo
   write(30,'(A4)')'&END'

   print 100,'Coordinates are saved in coord.cpmd!'

100 format(2X,A)

   close(30)

   return
end

!
!-------------------------------------------------------------------------
!  ------             WRITE DATA TO FILE   cp2k input FORMAT      ------ -
!-------------------------------------------------------------------------
!

subroutine write_cp2k
    use mode,only: x,npart,types,ctype,a,b,c,r1,r2,r3,alpha,beta,gamm,&
                   nelem
    implicit none
!
!   Local Variables
!
    integer         i
    integer         j
    integer         nt(4)
    real(8)         cy
    real(8)         pi,d2r,gamm_r,alpha_r,beta_r

    pi = dacos(-1.0d0)
!
    open(30,file='inp-cp2k',status='unknown')
    

    d2r = pi/180.0
    gamm_r = gamm*d2r     ! Degree to Rad
    alpha_r = alpha*d2r
    beta_r = beta*d2r


    call proj_coord
    !cy= r3*dcos(alpha_r)
    a = (/r1,0.0d0,0.0d0/)
    b = (/r2*dcos(gamm_r),r2*dsin(gamm_r),0.0d0/)
    !print *,gamm,beta
    cy= r3*dcos(alpha_r)-r3*dcos(beta_r)*dcos(gamm_r)
    cy= cy/dsin(gamm_r)
    c = (/r3*dcos(beta_r),cy,dsqrt(r3*r3*dsin(beta_r)*dsin(beta_r)-cy*cy)/)
    call cart_coord

    write(30,'(A7)')'&GLOBAL'
    write(30,'(2X,A11)')'RUN_TYPE MD'
    write(30,'(2X,A15)')'PRINT_LEVEL LOW'
    write(30,'(2X,A8)')'&TIMINGS'
    write(30,'(2X,A4)')'&END'
    write(30,'(A11)')'&END GLOBAL'
    write(30,*)

    write(30,*)
    write(30,'(A7)')'&MOTION'
    write(30,'(2X,A3)')'&MD'
    write(30,'(3X,A12)')'ENSEMBLE NVE'
    write(30,'(3X,A8)')'STEPS 10'
    write(30,'(3X,A12)')'TIMESTEP 0.5'
    write(30,'(3X,A17)')'TEMPERATURE 300.0'
    write(30,'(2x,A7)')'&END MD'
    write(30,'(A11)')'&END MOTION'
    write(30,*)

    write(30,*)
    write(30,'(A11)')'&FORCE_EVAL'
    write(30,'(2X,A9)')'METHOD QS'
    write(30,'(2X,A4)')'&DFT'
    write(30,'(2X,A36)')'BASIS_SET_FILE_NAME GTH_BASIS_SETS'
    write(30,'(2X,A32)')'POTENTIAL_FILE_NAME POTENTIAL'
    write(30,'(2X,A8)')'&END DFT'



    write(30,'(2X,A7)')'&SUBSYS'

    nt(1)= 0
    nt(2)= 0
    nt(3)= 0
    nt(4)= 0

    do i=1,npart
       if(types(i)==1)nt(1)=nt(1)+1
       if(types(i)==2)nt(2)=nt(2)+1
       if(types(i)==3)nt(3)=nt(3)+1
       if(types(i)==4)nt(4)=nt(4)+1
    enddo


    do j=1,nelem
       write(30,'(3X,A5,1X,A2)')'&KIND',adjustl(ctype(j))
       write(30,'(4X,A20)')'BASIS_SET TZV2P-GTH'
       write(30,'(4X,A22)')'POTENTIAL GTH-PADE-q1'
       write(30,'(3X,A9)')'&END KIND'
    enddo

    write(30,'(3X,A5)')'&CELL'
    write(30,'(4X,A4,3F12.5)')'ABC',r1,r2,r3
    write(30,'(3X,A9)')'&END CELL'

    write(30,'(3X,A6)')'&COORD'
    do j=1,nelem
       !write(30,'(A1,A2)')'*',adjustl(ctype(j))
       !write(30,'(I5)')nt(j)
       do i=1,npart
          if(types(i)==j)write(30,'(A2,3F20.8)')adjustl(ctype(j)),x(1,i),x(2,i),x(3,i)
       enddo
    enddo
    write(30,'(3X,A10)')'&END COORD'
    write(30,'(2X,A11)')'&END SUBSYS'
    write(30,'(A15)')'&END FORCE_EVAL'
    write(30,*)

    print 100,'Coordinates are saved in inp-cp2k!'

100 format(2X,A)

    close(30)

    return
end

!
!
!
!-------------------------------------------------------------------------
!

    subroutine write_map
    use lattice
    use mode,only: cpart,sc,npart
    implicit none

    integer      i,j

    open(112,file='map.in',status= 'unknown')

    write(112,'(4I4)')(sc(j),j=1,3),cpart
    write(112,'(A22)')'generated by symmetry'
    
    do i=1,npart
       write(112,'(4I3,I7)')(LatticeIndex(i,j),j=1,3),Natom(i),i
    enddo

    close(112)

    return
    end


!
!-------------------------------------------------------------------------
!            Write out cfg for atomeye viewer                        -
!-------------------------------------------------------------------------
!

subroutine write_cfg
    use mode,only:npart,ctype,types,xd,mass,a,b,c,ialive,alive !,x
    !use mod_pbc,only: a,b,c
    implicit none
!
!   Local Variables
!
    !integer lent
    integer i
    !integer ialive
!
!   
    !ialive = 0
    open(44,file= 'card.cfg', status='unknown')
!
    do i=1, npart
       if(alive(i)==1) ialive =ialive + 1
    enddo
    write(44,'(A22,I6)')'Number of particles = ',ialive
    write(44,*)
    write(44,'(A37)') 'A = 1.0 Angstrom (basic length-scale)'
    write(44,*)
    !print *,'in cfg ...'
    write(44,'(A10,F11.5,A2)')'H0(1,1) = ',a(1),' A'
    write(44,'(A10,F11.5,A2)')'H0(1,2) = ',a(2),' A'
    write(44,'(A10,F11.5,A2)')'H0(1,3) = ',a(3),' A'
    write(44,*)
    write(44,'(A10,F11.5,A2)')'H0(2,1) = ',b(1),' A'
    write(44,'(A10,F11.5,A2)')'H0(2,2) = ',b(2),' A'
    write(44,'(A10,F11.5,A2)')'H0(2,3) = ',b(3),' A'
    write(44,*)
    write(44,'(A10,F11.5,A2)')'H0(3,1) = ',c(1),' A'
    write(44,'(A10,F11.5,A2)')'H0(3,2) = ',c(2),' A'
    write(44,'(A10,F11.5,A2)')'H0(3,3) = ',c(3),' A'
    write(44,*)
    !print *,'projecting to fraction coordinate...'
    call proj_coord
    !print *,'fraction coordinate...'
    do i=1,npart
           if(alive(i)/=0) write(44,500)mass(types(i)),ctype(types(i)),xd(1,i),xd(2,i),xd(3,i),0,0,0
    enddo
!
500 FORMAT(F9.4,A3,3F20.6,3I2)
    
    close(44)
!
    return
end
!


!
!-------------------------------------------------------------------------
!            Write out cif for crystal information                       -
!-------------------------------------------------------------------------
!

subroutine write_cif
    use mode,only:npart,ctype,types,xd,ialive,alive, &
                & r1,r2,r3,alpha,beta,gamm !,mass ,x
    !use mod_pbc,only: a,b,c
    implicit none
!
!   Local Variables
!
    !integer lent
    integer i
    integer ne(8)
    !integer ialive
    character(8)  id
    character(2)  ctyp
!
!   
    !ialive = 0
    open(44,file= 'card.cif', status='unknown')
!
    do i=1, npart
       if(alive(i)==1) ialive =ialive + 1
    enddo
    write(44,'(A11)') 'data_image0'
    write(44,'(A14,1x,F12.8)') '_cell_length_a', r1 
    write(44,'(A14,1x,F12.8)') '_cell_length_b', r2
    write(44,'(A14,1x,F12.8)') '_cell_length_c', r3 
    write(44,'(A17,1x,F12.8)') '_cell_angle_alpha', alpha 
    write(44,'(A16,1x,F12.8)') '_cell_angle_beta', beta 
    write(44,'(A17,1x,F12.8)') '_cell_angle_gamma', gamm 
    write(44,*)
    !print *,'in cfg ...'
    write(44,'(A5)')'loop_'
    write(44,'(A17)')' _atom_site_label'
    write(44,'(A21)')' _atom_site_occupancy'
    write(44,'(A19)')' _atom_site_fract_x'
    write(44,'(A19)')' _atom_site_fract_y'
    write(44,'(A19)')' _atom_site_fract_z'
    write(44,'(A33)')' _atom_site_thermal_displace_type'
    write(44,'(A26)')' _atom_site_B_iso_or_equiv'
    write(44,'(A23)')' _atom_site_type_symbol'

    !print *,'projecting to fraction coordinate...'
    call proj_coord
    !print *,'fraction coordinate...'
    ne = 0
    do i=1,npart
       ne(types(i)) = ne(types(i)) + 1
       write(id,'(I8)') ne(types(i))
       id=adjustl(id)
       ctyp=adjustr(ctype(types(i)))
       if(alive(i)/=0) write(44,500)ctyp,id,1.00,xd(1,i),xd(2,i),xd(3,i),' Biso ',1.00,ctype(types(i))
    enddo
!
500 format(A2,A8,F5.2,3F12.8,A6,F5.2,1X,A3)
    close(44)
!
    return
end

