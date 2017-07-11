!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!------                    Read in configure                        ------
!-------------------------------------------------------------------------
!

subroutine control
     use mode,   only: cpart, sc,npart,a0,b0,c0,&
                       movex,movey,movez,ctype,nelem,r1,r2,r3,alpha,beta,gamm,vel,mass,&
                       shockv,struc_name !,a,b,c ,x
     use modctrl,only: center,iwr,irec,cfile,ilt,istretch,fzoom,lmap,lcut,fremove,&
                       ldipole,lmole, run_type,istretchfreebond
     use modgr,  only: icom1,icom2
     use mod_stretch
     use mod_rotate
     use mod_nw,only: latmtyp,basis
     use mod_def,only: idef,defect
     use mod_mmfit,only: task
     use mod_grd,only:grid_op
     implicit none
!
     integer        i,j !,id
     integer        nat
     !integer        typ
     integer        iline
!     integer        lt
     !character      cty
     character(80)  line
     character(40)  c1,c2,c3,c4
     character(40)  coorf
     !character(40)  struc_name
     real(8)        cy,ay
     real(8) ::     pi = 3.1415927

     logical ::     lstruc=.FALSE.

     pi = dacos(-1.0d0)
!     character(40)  key
!
!
!     open(10,file='in.txt',status='old')
!
     iline = 0
     do
        read(5,200,end=99,err=999)line
        iline= iline + 1
        line = trim(adjustl(line))

        if(len(line)==0)cycle
        if(index(line,'#').EQ.1)cycle

        if(index(line,'%coord').NE.0)then
            read(line,*,err=999)c1,coorf
!            print *,c1,coorf
        elseif(index(line,'%structure').NE.0)then
            read(line,*,err=999)c1,struc_name
            lstruc =.true.
        elseif(index(line,'%molecularstructure').NE.0)then
            read(line,*,err=999)c1,struc_name
            lstruc =.true.
            lmole = .true.
            irec  = .true.
        elseif(index(line,'%file').NE.0)then
            read(line,*,err=999)c1,c2
            cfile=trim(adjustl(c2))
        elseif(index(line,'%center').NE.0)then
            read(line,*,err=999)c1,center
        elseif(index(line,'%supercell').NE.0)then
            read(line,*,err=999)c1,(sc(i),i=1,3)
        elseif(index(line,'%totgrid').NE.0)then
            grid_op = .true.
        elseif(index(line,'%cellpart').NE.0)then
            read(line,*,err=999)c1,cpart
        elseif(index(line,'%totpart').NE.0)then
            read(line,*,err=999)c1,npart
        elseif(index(line,'%recover').NE.0)then !iptx,ipty,iptz: 是否将分裂的原子再放在一起
            read(line,*,err=999)c1,irec
        elseif(index(line,'%output').NE.0)then
            read(line,*,err=999)c1,c2
            c2= trim(adjustl(c2))
!            print '(A4)',c2
            if(c2=='lammpsdata')then
               iwr=0
            elseif(c2=='xyz')then
               iwr=1
            elseif(c2=='geo')then
               iwr=2
            elseif(c2=='NEICAR')then
               iwr=3
            elseif(c2=='analysis')then
               iwr=4
            elseif(c2=='gulp')then
               iwr=5
               read(line,*,err=999)c1,c2,c3
               if(c3=='opt')    run_type = 0
               if(c3=='md-nvt') run_type = 1
               if(c3=='md-npt') run_type = 2
               if(c3=='gradient') run_type = 3
            elseif(c2=='pdb')then
               iwr=6
            elseif(c2=='dlpolyconfig')then
               iwr=7
            elseif(c2=='ffield'.or.c2=='FFIELD'.or.c2=='gulp_field')then
               read(line,*,err=999)c1,c2,c3
               iwr=8
               if (c3=='atom_type') latmtyp = .true.
            elseif(c2=='gofr')then
               iwr=9
               read(line,*,err=999)c1,c3,icom1,icom2
            elseif(c2=='POSCAR'.or.c2=='poscar')then
               iwr= 10
            elseif(c2=='SIESTA'.or.c2=='siesta')then
               iwr= 11
            elseif(c2=='dftb'.or.c2=='DFTB')then
               iwr= 12
            elseif(c2=='cpmd'.or.c2=='CPMD')then
               iwr= 13
               read(line,*,err=999)c1,c2,c3
               if(c3=='opt') run_type = 0
               if(c3=='nve') run_type = 1
               if(c3=='nvt') run_type = 2
               if(c3=='npt') run_type = 3
            elseif(c2=='cfg'.or.c2=='CFG')then
               iwr= 15
            elseif(c2=='table'.or.c2=='TABLE')then
               iwr= 16
            elseif(c2=='moltop'.or.c2=='MOLTOP')then
               iwr= 17
            elseif(c2=='cp2k'.or.c2=='CP2K')then
               iwr= 18
            elseif(c2=='lammpsshock')then
               iwr=19
               read(line,*,err=999)c1,c2,shockv
            elseif(c2=='nwchem'.or.c2=='NWCHEM')then
               iwr= 21
               read(line,*,err=999)c1,c2,basis
            elseif(c2=='ZMATRIX'.or.c2=='zmatrix')then
               iwr= 22
            elseif(c2=='gulp_ff'.or.c2=='gulp_force_field'.or.c2=='gulp_ffield')then
               iwr= 23
               read(line,*,err=999)c1,c2,c3
               if(c3=='opt')    run_type = 0
               if(c3=='md-nvt') run_type = 1
               if(c3=='md-npt') run_type = 2
               if(c3=='gradient') run_type = 3
            elseif(c2=='tinker'.or.c2=='TINKER')then
               iwr= 24
            elseif(c2=='gromacs'.or.c2=='GROMACS')then
               iwr= 25
               read(line,*,err=999)c1,c2,c3
               if (c3=='atom_type') then
                  latmtyp = .true.
               else
                  latmtyp = .false.
               endif   
            elseif(c2=='uspex'.or.c2=='USPEX')then
               iwr= 26
            elseif(c2=='mmfit'.or.c2=='MMFIT')then
               iwr= 27
               read(line,*,err=999)c1,c2,c3
               if(c3 == 'bat')   task = 1
               if(c3 == 'check') task = 2
            elseif(c2=='species')then
               iwr=28
               latmtyp = .true.
            elseif(c2=='cif'.or.c2=='CIF')then
               iwr= 29
            elseif(c2=='gro'.or.c2=='GRO')then
               iwr= 30
            elseif(c2=='database'.or.c2=='DATABASE')then
               iwr= 99
            elseif(c2=='autorotate'.or.c2=='auto_rotate')then
               iwr= 100
            elseif(c2=='autoswing'.or.c2=='auto_swing')then
               iwr= 101
            else 
               call err('ERROR: Output type not suported!')
               print '(2x,A)','In line:'
               print '(I6)',iline
            end if
        elseif(index(line,'%cellvector').NE.0)then
            read(5,*,err=999)(a0(i),i=1,3)
            read(5,*,err=999)(b0(i),i=1,3)
            read(5,*,err=999)(c0(i),i=1,3)
        elseif(index(line,'%MappingLattice').NE.0)then
            read(line,*,err=999)c1,lmap
        elseif(index(line,'%cellpara').NE.0)then
            read(line,*)c1,ilt
            read(5,*,err=999)r1,r2,r3,alpha,beta,gamm
            alpha= alpha*pi/180.0d0
            beta = beta*pi/180.0d0
            gamm= gamm*pi/180.0d0
            print '(A80)',&
                  '--------------------------------------------------------------------------------'
            if(ilt==1)then
                     print '(2x,A)','Defult: a in X axis, b in XY Plane, c arbitary!'
                     a0 = (/r1,0.0d0,0.0d0/)
                     b0 = (/r2*dcos(gamm),r2*dsin(gamm),0.0d0/)
                     cy= r3*dcos(alpha)-r3*dcos(beta)*dcos(gamm)
                     cy= cy/dsin(gamm)
                     c0 = (/r3*dcos(beta),cy,dsqrt(r3*r3*dsin(beta)*dsin(beta)-cy*cy)/)
            elseif(ilt==2)then
!                    c along Z, b in YZ plane
                     print '(2x,A)','For Reax geo: c along Z, b in YZ plane, a arbitary!'
                     c0 = (/0.0d0,0.0d0,r3/)
                     b0 = (/0.0d0,r2*dsin(alpha),r2*dcos(alpha)/) 
                     ay= r1*dcos(gamm)-r1*dcos(alpha)*dcos(beta)
                     ay= ay/dsin(alpha)
                     a0 = (/dsqrt(r1*r1*dsin(beta)*dsin(beta)-ay*ay),ay,r1*dcos(beta)/)
            endif
        elseif(index(line,'%move').NE.0)then
            read(line,*,err=999)c1,movex,movey,movez
        elseif(index(line,'%element').NE.0)then       ! must have this commond~
            ctype = '  '
            read(line,*,err=99)c1,nelem,(ctype(i),i=1,nelem)
        elseif(index(line,'%masses').NE.0)then
            read(line,*,err=99)c1, (mass(i),i=1,nelem)
        elseif(index(line,'%zoom').NE.0)then
            fzoom = .true.
        elseif(index(line,'%remov').NE.0)then
            fremove = .true.
        elseif(index(line,'%dipole').NE.0)then
            ldipole = .true.
!  remove h2o molecules or other 
        elseif(index(line,'%stretch').NE.0)then    
!  read in stretch bonds
!  usage:

!  stretch 3 bonds
!  nwchem or gaussian
!  dr 0.02 elongation 8 shorten 8
!  fix-end 1  move-end 3 to-be-moved 2 atoms 5 4 
            read(line,*)c1,stretch_name,c2,c3 !stretch_name: head of file name 
            if (c2=='comp') comp_var = .true.
            if (c3=='free') then
               istretchfreebond = .true.
               istretch = .false.
            elseif(c3/='free') then
               istretch = .true.
               istretchfreebond = .false.
            endif
            read(5,*)inptype,species1,species2,basis
            read(5,*)c1,mv_dr,c2,nenlongation,c3,nshorten
            read(5,*,err=999)c1,mv%a,c2,mv%b,c3,mv%num,c4,(mv%ind(j),j=1,mv%num)
        elseif(index(line,'%rotate').NE.0)then    
            irotate = .true.
            read(line,*)c1,inptype,basis!,nrotate,atom1,atom2,atom3,num_rot,drot  !!! nrotate: n atoms to rotate
            read(5,*,err=999) nrotate,atom1,atom2,atom3,num_rot,drot
            do i=1,nrotate
               read(5,*,err=999) rota(i)!,dih(i)  ! rota: atom to rotate
            enddo
        elseif(index(line,'%swing').NE.0)then    
            iswing = .true.
            read(line,*)c1,inptype,basis
            read(5,*,err=999) nrotate,atom1,atom2,atom3,num_rot,drot  !!! nrotate: n atoms to rotate
            do i=1,nrotate
               read(5,*,err=999) rota(i)!,dih(i)  ! rota: atom to rotate
            enddo
        elseif(index(line,'%defect').NE.0)then   ! creat circular defects
            read(line,*)c1, idef ! idef: # of defect
            allocate(defect(4,idef))
            do i=1,idef
               !centre of the circle and radius of the circle
               read(5,*)defect(1,i),defect(2,i),defect(3,i),defect(4,i) 
            enddo
        elseif(index(line,'%develop').NE.0)then
            read(line,*)c1, lcut,vel
        end if
     end do

999  call err('Read Input File Error, In Line:')
     print '(I6)',iline
99   continue
     
     if(coorf=='pdb'.or.coorf=='PDB')then
          open(11,file=cfile,status='old')
          nat = 0
          do
             read(11,'(A80)',end=9)line
             if(index(line,'ATOM')/=0)nat=nat+1
          enddo
9         continue
          cpart= nat
          close(11)
     endif

! -------------- initialize memory --------------

     if (lstruc) then
! @@@@@@@@@@@@@@@@@@@@@ calling database @@@@@@@@@@@@@@@@@@@@

        if(struc_name=='tatb'.or.struc_name=='TATB') call tatb
        if(struc_name=='n2co'.or.struc_name=='N2CO'.or.struc_name=='CON2'.or.struc_name=='con2')   call n2co
        if(struc_name=='rdx'.or.struc_name=='RDX')   call rdx
        if(struc_name=='nitromethane'.or.struc_name=='NM')   call nitromethane
        if(struc_name=='hmx-cl20'.or.struc_name=='cl20-hmx')   call hmx_cl20
        if(struc_name=='hmx'.or.struc_name=='HMX')  call hmx
        if(struc_name=='cl20'.or.struc_name=='CL20')  call cl20
        if(struc_name=='cl20mol'.or.struc_name=='molcl20')  call cl20mol
        if(struc_name=='hmxmol'.or.struc_name=='molhmx')  call hmxmol
        if(struc_name=='ch4'.or.struc_name=='methane'.or.struc_name=='CH4')  call methane
        if(struc_name=='nh3'.or.struc_name=='ammonia'.or.struc_name=='NH3')  call nh3
        if(struc_name=='nmmol'.or.struc_name=='molnm'.or.struc_name=='nitromethanemol')  call nmmol
        if(struc_name=='nitroethane')  call nitroethane
        if(struc_name=='ethane')  call ethane
        if(struc_name=='ethylene' .or. struc_name=='ch2ch2')  call ethylene
        if(struc_name=='no2' .or. struc_name=='NO2')  call no2
        if(struc_name=='co2' .or. struc_name=='CO2')  call co2
        if(struc_name=='c3h6n2o4')  call c3h6n2o4
        if(struc_name=='fox7'.or.struc_name=='fox-7'.or.struc_name=='FOX-7')  call fox7
        if(struc_name=='Si'.or.struc_name=='si')  call si

! @@@@@@@@@@@@@@@@@@@@ calling database @@@@@@@@@@@@@@@@@@@@@
     else
         npart=cpart*sc(1)*sc(2)*sc(3)
         call memory
     endif
! -------------- initialize memory --------------

     coorf= trim(adjustl(coorf))
     if(coorf=='POSCAR'.or.coorf=='poscar')then
        call readposcar
     elseif(coorf=='lammpstrj'.or.coorf=='LAMMPSTRJ')then
        call readltr
     elseif(coorf=='xyz'.or.coorf=='XYZ')then
        call readxyz
     elseif(coorf=='car'.or.coorf=='CAR')then
        call readcar
     elseif(coorf=='CONTCAR'.or.coorf=='contcar')then
        call read_contcar
     elseif(coorf=='lammpsdata')then
        call read_lmpdata
     elseif(coorf=='pdb'.or.coorf=='PDB')then
        call read_pdb
     elseif(coorf=='gen')then
        call read_gen  ! DFTB+ .gen file format
     elseif(coorf=='gulp'.or.coorf=='GULP')then
        call read_gulp
     !else
        !call err('ERROR: coord file type not supportted!')
     end if
!
!
!100  format(2x,A)
200  format(A80)
!
!     close(5)
     return
end





