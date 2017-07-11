!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!

subroutine auto_swing
  use mode,only: x,types,ctype,npart
  use mod_stretch,only: mv,stretch_name, &
                        species1,species2
  use modnei,only: bond,nbond 
  use modfim, only: cflag
  use mod_nw,only: atom_type  
  use modnei, only: atable
  implicit none

  integer       i,j,is,ie,c,ii,jj
  integer       lent,lent1
  integer       iatom,jatom,iring,jring
  integer       ax(2)
  integer,allocatable::ring(:)
  !real(8)       dr
  real(8),allocatable::xold(:,:)
  character(80) tbm
  character(5)  i2c

  character(10),external::int_to_str
  logical    :: enl = .false.


  stretch_name = '                              '

  !that atoms at a bond end or in a group that not bonded to others

  allocate(ring(npart))
  allocate(xold(3,npart))

  xold = x

  ring = 0

  if(.not. allocated(atom_type)) then
     allocate(atom_type(npart))
     atom_type  = 'null '

     do i=1,npart
        atom_type(i) = ctype(types(i))
     enddo
  endif
  !dr = mv_dr
  

  call neighbor
  call findmole
  call bonds
  call atm_typ


  c = 0
  species1 = adjustl(species1)
  species2 = adjustl(species2)
  lent = len_trim(species1)
  lent1= len_trim(species2)
  stretch_name(1:3) = 'ar_'
  stretch_name(4:3+lent) = species1
  stretch_name(4+lent:4+lent) = '_'
  stretch_name(5+lent:4+lent+lent1) = species2
  lent= 5+lent+lent1
  stretch_name(lent:lent) = '_'


  do i=1,nbond
     do ii =1,2
        if (ii ==1) then
           iatom = bond(1,i)
           jatom = bond(2,i)
        else
           iatom = bond(2,i)
           jatom = bond(1,i) 
        endif
        !print *,'bond',i,': ',ctype(types(iatom)),ctype(types(jatom))
        i2c=int_to_str(c)
        stretch_name(lent+1:lent+1)=i2c
        stretch_name(lent+2:lent+2)='_'
        mv%a=iatom
        mv%b=jatom


        cflag(:) = 1
        cflag(iatom) = 0 !

        mv%num = 1              ! jaotm is to be rotated or swinged
        mv%ind(mv%num) = jatom  ! 1st is jatom

        iring = 0
        jring = 0

        call track(iatom,jatom,iring,jring,jatom)  !!!! not state iring,jring

        ring(iatom) = iring
        ring(jatom) = jring 

        enl = .False.
        if (iring == 1 .and. jring ==1) enl = .true. ! if ring structure
        !print *,'bond: ',i,iatom,jatom,enl
        tbm = '                                                                            '
        if (.not. enl) then ! if not ring structure
           !print *,'bond: ',i,iatom,jatom,mv%num
           if (mv%num>=int(0.5*npart) .or. mv%num<2)  then
              !print *,'bond: ',i,iatom,jatom
              cycle
           endif

           do j = 1,mv%num
              i2c=int_to_str(mv%ind(j))
              is = len_trim(tbm)+2
              ie = len_trim(tbm)+1 + len_trim(i2c)
              if (ie>80) exit
              tbm(is-1:is-1) = ' '
              tbm(is:ie) = i2c
           enddo

           ! find axis
           if (atable(1,iatom)/=3) cycle
           jj = 0
           do j = 1,3
              if (atable(j+1,iatom)/=jatom) then
                 jj = jj +1
                 if (jj>=3) print *,'Error! atoms bond to iatom large than 3'
                 ax(jj) = atable(j+1,iatom)
              endif
           enddo 

           print '(2I4,A1,A2,A1,A2,A6,I4,A2,A6,2I4,A16,A)',iatom,jatom,' ',adjustr(ctype(types(iatom))),'-',&
               adjustl(ctype(types(jatom))), ' swing',i,': ','axis: ',ax(1),ax(2),' to be swinged  ',tbm

           !k = k + 1 + nenlongation + nshorten
           !call stretch
           c = c + 1
        endif
     enddo
  enddo

return
end  subroutine auto_swing

!
!==============================================================================
!

subroutine swing
   use mode,only: x
   use mod_rotate,only: nrotate,atom1,atom2,atom3,num_rot,rota,drot
   implicit none

   integer          i,j
   character(16) :: inp = '                '
   character(16) :: cc = ''
   real(8)          xold(3,20)


   call  neighbor
   call  findmole

   call  bonds
   call  angles
   call  torsions
   
   call  gulp_ff

   inp(1:5) = 'tors_'
   inp(6:16) = '          '

   do i=1, nrotate
      xold(1,i) = x(1,rota(i))
      xold(2,i) = x(2,rota(i))
      xold(3,i) = x(3,rota(i))
   enddo
   
   !print *,drot
   do i=1, nrotate
      !print *,rota(i)
      !print *,drot
      call swing_atom(atom1,atom2,atom3,rota(i), drot)
   enddo

   print '(A21)','** Variables changed: '
   call w_bond
   call w_angle
   call w_torsion

   do i=1, nrotate               ! go back to oringinal
      x(1,rota(i)) = xold(1,i)
      x(2,rota(i)) = xold(2,i)
      x(3,rota(i)) = xold(3,i)
   enddo

   do j=0, num_rot
      do i=1, nrotate
         call swing_atom(atom1,atom2,atom3,rota(i),j*drot)
      enddo

      call writexyz('mol_config')

      inp(6:16) = '          '
      !inp(6:6) = '+'

      write(cc,'(I8)') j
      cc = trim(adjustl(cc))
      if(j<10) then
         inp(6:6) = cc
         inp(7:10)= '.inp'
      elseif(j<100.and.j>9)then
         inp(6:7) = cc
         inp(8:11)= '.inp'
      else
         inp(6:8) = cc
         inp(9:12)= '.inp'
      endif
      call write_nwchem(inp)
      
   enddo

   do i=1, nrotate               ! go back to oringinal
      x(1,rota(i)) = xold(1,i)
      x(2,rota(i)) = xold(2,i)
      x(3,rota(i)) = xold(3,i)
   enddo


   do j=1, num_rot
      do i=1, nrotate
         !print *,-j*drot
         call swing_atom(atom1,atom2,atom3,rota(i),-j*drot)
      enddo
      call writexyz('mol_config')

      inp(6:16) = '          '
      inp(6:6) = '-'
      
      write(cc,'(I8)') j
      cc = trim(adjustl(cc))
      if(j<10) then
         inp(7:7) = cc
         inp(8:11)= '.inp'
      elseif(j<100.and.j>9)then
         inp(7:8) = cc
         inp(9:12)= '.inp'
      else
         inp(7:9) = cc
         inp(10:13)= '.inp'
      endif
      call write_nwchem(inp)

   enddo

   return
end
!                                  
!=========================================================================
!

subroutine swing_atom(atom1,atom2,atom3,atom4,dihe)
!
!
!          ^
!          |
!          | 
!          |              /  axis 
!          |             /  
!          |    atom3   /
!          |           *              atom4  to be swing
!          |          /              /
!          |         /              /
!          |        /     ------- atom1
!          |       * atom2          \
!          |      /                  \
!          |     /                    \   swing of a group atoms
!          |                          atom4' to be swing 
!          |-------------------------->
!
   use mode,only: x
   implicit none

   !dih is target torsion angle
   !--atom2--atom3--  the axis
   !--atom1: point move the axis to

   integer        atom1,atom2,atom3,atom4
   real(8)        dihe,dih
   real(8)        v0(3),vt(3) ! vt target vector
   real(8)        v2(3),v3(3)
   real(8)        vi(3),vj(3)! ! **** vk = uv2 ****
   real(8)        uv2(3) 
   real(8)        b23, rdot,v4(3),bv4
   real(8)     :: pi=3.1415926535897932
   !real(8),external:: torsion_angle

   dih = dihe*pi/180.0000      ! degree to rad

   ! compute torsion angle

   v2(1)=x(1,atom3)-x(1,atom2)
   v2(2)=x(2,atom3)-x(2,atom2)
   v2(3)=x(3,atom3)-x(3,atom2)

   v3(1)=x(1,atom4)-x(1,atom1)
   v3(2)=x(2,atom4)-x(2,atom1)
   v3(3)=x(3,atom4)-x(3,atom1)

   !compute bond length
   b23 = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3)
   b23 = sqrt(b23)
   !print *,b23

   uv2(1) = v2(1)/b23
   uv2(2) = v2(2)/b23
   uv2(3) = v2(3)/b23

   !print *,uv2(1),uv2(2),uv2(3)

   ! point of atom4 project on axis
   rdot = v3(1)*uv2(1) + v3(2)*uv2(2) + v3(3)*uv2(3)
   !print *,rdot
   v0(1) = x(1,atom1) + rdot*uv2(1)
   v0(2) = x(2,atom1) + rdot*uv2(2)
   v0(3) = x(3,atom1) + rdot*uv2(3)

   v4(1) = x(1,atom4) - v0(1)
   v4(2) = x(2,atom4) - v0(2)
   v4(3) = x(3,atom4) - v0(3)

   bv4 = v4(1) * v4(1) + v4(2)*v4(2) + v4(3)*v4(3)
   bv4 = sqrt(bv4)
   !if(atom4==23) print *,bv4,rdot

   vi(1) = v4(1)/bv4
   vi(2) = v4(2)/bv4
   vi(3) = v4(3)/bv4 

   vj(1) = uv2(2)*vi(3) - uv2(3)*vi(2)
   vj(2) = uv2(3)*vi(1) - uv2(1)*vi(3)
   vj(3) = uv2(1)*vi(2) - uv2(2)*vi(1)

   ! target vector
   vt(1) = vi(1)*cos(dih) + vj(1)*sin(dih)
   vt(2) = vi(2)*cos(dih) + vj(2)*sin(dih)
   vt(3) = vi(3)*cos(dih) + vj(3)*sin(dih)
   !print *,cos(dih)

   ! target coordinate
   x(1,atom4) = v0(1) + vt(1)*bv4 
   x(2,atom4) = v0(2) + vt(2)*bv4 
   x(3,atom4) = v0(3) + vt(3)*bv4 
   !print *,x(1,atom4),x(2,atom4),x(3,atom4)

   return
end