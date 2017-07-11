!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                     email: gfeng.alan@foxmail.com                      =
!=========================================================================
!

subroutine stretch_free_bond
  use mode,only: x,types,ctype,npart
  use mod_stretch,only: nenlongation,nshorten,mv_dr,mv,stretch_name, &
                        species1,species2 !,comp_var
  !use mod_rotate,only: inptype
  use modnei,only: bond,nbond!,bond_value
  use modfim, only: cflag
  use mod_nw,only: atom_type  
  implicit none

  integer       i,j,k,is,ie,c!,jj ,ws
  integer       lent,lent1
  integer       iatom,jatom,iring,jring
  integer,allocatable::ring(:)
  real(8)       dr
  real(8),allocatable::xold(:,:)
  character(80) tbm
  character(5)  i2c
  !character(1)  ityp,jtyp
  character(10),external::int_to_str
  logical    :: enl = .false.

  !logical,allocatable::(bond_as(:)) ! already have stretched
  !logical,allocatable::(tobes(:))   ! atoms to be stretch(move): 

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
  dr = mv_dr
  

  call neighbor
  call findmole
  call bonds
  call atm_typ

  k = 0
  c = 0
  species1 = adjustl(species1)
  species2 = adjustl(species2)
  lent = len_trim(species1)
  lent1= len_trim(species2)
  stretch_name(1:3) = 'fb_'
  stretch_name(4:3+lent) = species1
  stretch_name(4+lent:4+lent) = '_'
  stretch_name(5+lent:4+lent+lent1) = species2
  lent= 5+lent+lent1
  stretch_name(lent:lent) = '_'

  !print *, stretch_name
  !jj = 0
  do i=1,nbond

     iatom = bond(1,i)
     jatom = bond(2,i)
     
     !print *,iatom,atom_type(iatom),jatom,atom_type(iatom)

     if( (atom_type(iatom)==species1 .and. atom_type(jatom)==species2) .or. &
        (atom_type(iatom)==species2 .and. atom_type(jatom)==species1)  ) then 
        
        !print *,atom_type(iatom),atom_type(jatom)
     
        i2c=int_to_str(c)
        stretch_name(lent+1:lent+1)=i2c
        stretch_name(lent+2:lent+2)='_'
        mv%a=iatom
        mv%b=jatom

        !if (ctype(types(iatom)) == 'H' .or. ctype(types(jatom)) =='H')then
           !mv_dr = 0.5*dr
        !else
        !   mv_dr = dr
        !endif

        cflag(:) = 1
        cflag(iatom) = 0 !

        mv%num = 1
        mv%ind(mv%num) = jatom
        iring = 0
        jring = 0
        !print *,iatom,jatom,ctype(types(iatom)),ctype(types(jatom))

        call track(iatom,jatom,iring,jring,jatom)  !!!! not state iring,jring

        ring(iatom) = iring
        ring(jatom) = jring 
        if (iring == 1 .and. jring ==1) enl = .true. ! if ring structure
        
        !print *,iring,jring,'there'

        tbm = '                                                                            '
        if (.not. enl) then ! if not ring structure
           do j = 1,mv%num
              i2c=int_to_str(mv%ind(j))
              is = len_trim(tbm)+2
              ie = len_trim(tbm)+1 + len_trim(i2c)
              if (ie>80) exit
              tbm(is-1:is-1) = ' '
              tbm(is:ie) = i2c
           enddo

           print '(2I4,A1,A2,A1,A2,A5,I4,A2,A13,A80)',iatom,jatom,' ',ctype(types(iatom)),'-',&
               ctype(types(jatom)), ' bond',i,': ','to be moved  ',tbm
           !print '(2I4,A1,A2,A1,A2,2I6)',iatom,jatom,' ',ctype(types(iatom)),'-',&
           !    ctype(types(jatom)), k+1,k + 1 + nenlongation + nshorten

           k = k + 1 + nenlongation + nshorten
           call stretch
           c = c + 1
        endif
     endif
     !print *,'there'
    !jj = jj+1
    !print *,jj
  enddo
  !print *,'tmd!',jj
  !if (enl) ! if ring structure

  return
end

subroutine stretch_config
  use mode,only: x,types,ctype,npart
  use mod_stretch,only: nenlongation,nshorten,mv_dr,mv
  use mod_rotate,only: inptype
  use modnei,only: bond,nbond !,bond_value
  use modfim, only: cflag
  use mod_nw,only: atom_type  
  implicit none

  integer       i,j,k,is,ie
  integer       iatom,jatom,iring,jring
  integer,allocatable::ring(:)
  real(8)       dr
  real(8),allocatable::xold(:,:)
  character(80) tbm
  character(5)  i2c
  character(1)  ityp,jtyp
  character(10),external::int_to_str
  logical    :: enl = .false.

  !logical,allocatable::(bond_as(:)) ! already have stretched
  !logical,allocatable::(tobes(:))   ! atoms to be stretch(move): 


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
        !print *,atom_type(i)
     enddo
  endif
  dr = mv_dr

  if(inptype=='nwchem' .or. inptype=='gaussian')then
    call stretch 
  elseif (index(inptype,'xyz').NE.0)then
    call neighbor
    call findmole
    call bonds

    k = 0
    do i=1,nbond
       iatom = bond(1,i)
       jatom = bond(2,i)
       ityp = inptype(5:5)
       jtyp = inptype(6:6)
       if((ityp==ctype(types(iatom)).and.jtyp==ctype(types(jatom)) ).or. &
         (jtyp==ctype(types(iatom)).and.ityp==ctype(types(jatom)) ) .or. &
         index(inptype,'all').NE.0 .or. index(inptype,'ALL').NE.0) then 
            continue
       else
            cycle
       endif
       mv%a=iatom
       mv%b=jatom

       if (ctype(types(iatom)) == 'H' .or. ctype(types(jatom)) =='H')then
          mv_dr = 0.5*dr
       else
          mv_dr = dr
       endif

       cflag(:) = 1
       cflag(iatom) = 0 !

       mv%num = 1
       mv%ind(mv%num) = jatom
       iring = 0
       jring = 0
       !print *,iatom,jatom,ctype(types(iatom)),ctype(types(jatom))

       call track(iatom,jatom,iring,jring,jatom)  !!!! not state iring,jring

       ring(iatom) = iring
       ring(jatom) = jring 
       if (iring == 1 .and. jring ==1) enl = .true. ! ring structure
       !print *,iring,jring
       tbm = '                                                                            '
       if (.not. enl) then  ! if not ring structure
          do j = 1,mv%num
             i2c=int_to_str(mv%ind(j))
             is = len_trim(tbm)+2
             ie = len_trim(tbm)+1 + len_trim(i2c)
             tbm(is-1:is-1) = ' '
             tbm(is:ie) = i2c
          enddo

          print '(2I4,A1,A2,A1,A2,A5,I4,A2,A13,A80)',iatom,jatom,' ',ctype(types(iatom)),'-',&
              ctype(types(jatom)), ' bond',i,': ','to be moved  ',tbm
          print '(2I4,A1,A2,A1,A2,2I6)',iatom,jatom,' ',ctype(types(iatom)),'-',&
              ctype(types(jatom)), k+1,k + 1 + nenlongation + nshorten

          k = k + 1 + nenlongation + nshorten
          call stretch
       endif
    enddo
    !if (enl)  ! if ring structure
    if (enl) then
       do j=1,3
          do i=1,npart
             x(:,i) = xold(:,i) * (1.0+0.001*j)
          enddo 
          call writexyz('mol_config')
       enddo
       call w_bond
       do j=1,3
          do i=1,npart
             x(:,i) = xold(:,i) * (1.0-0.001*j)
          enddo 
          call writexyz('mol_config')
       enddo
    endif
  endif

  return
end

recursive subroutine track(aa,bb,iring,jring,m)
    use mod_stretch,only: mv
    use mode,only: ctype,types
    use modnei, only: atable,matom
    use modfim, only: cflag
    implicit none
!
    integer      atab(matom)
    integer      at(matom)
    integer      m, mm,i,tm,j
    integer      aa,bb
    integer      iring,jring
!

    atab = 0
    at = 0
    mm=m

    if (cflag(mm) == 0) return

    do i=1,atable(1,mm)+1
       atab(i)=atable(i,mm)
    enddo
!
    cflag(mm)=0
!
    do i=2,atab(1)+1
       if(atab(i)==0) print '(A9)','** error!'
       if(cflag(atab(i))/=0)then
          !cflag(atab(i))=0
          tm=atab(i)
          do j=2,atable(1,aa)+1
             if (tm==atable(j,aa) .and. tm /= bb) then
                print '(A51)','** bonds in circular structure cannot be stretched!'
                print '(A22,2I5,2A4)','** Atoms ID and type: ', &
                       aa,bb,ctype(types(aa)),ctype(types(bb))
                iring = 1
                jring = 1
                return
             endif
          enddo
          mv%num = mv%num + 1
          mv%ind(mv%num) = tm
          call track(aa,bb,iring,jring,tm)
       endif
    enddo
!
    return
end
!

subroutine stretch 
!
!
!          ^
!          |
!          | 
!          |          *    * 
!          |          |   /  
!          |          |  / 
!          |           * ----*   
!          |          /  
!          |         / 
!          |        /      
!          |      
!          |-------------------------->
!

    use mode,only: x
    use mod_stretch
    !use modnei,only: nbond!,bond_value
    use mod_rotate,only: inptype
    implicit none

!   type(move) :: nitr(3)

    integer        ::   i,k,kk
    integer             lent
    character(30)       filename
    !character(5)   ::   c1,c2
    real(8)             rr
    character(5),external::int_to_str
    !character(40),external::str_add

    filename = '                              '
    lent =len_trim(adjustl(stretch_name))
    !print *,stretch_name
    !print *,lent
    filename(1:lent) = stretch_name(1:lent)
    filename(lent+1:lent+1) = trim(adjustl(int_to_str(0)))
    !print *,filename
    filename(lent+2:lent+5) = '.inp'
    !print *,filename

    if(inptype=='nwchem')then
        call write_nwchem(filename)
    elseif(inptype=='gaussian')then
        call write_gaussian(filename)
    elseif(index(inptype,'xyz').NE.0)then
        call writexyz('mol_config')
    endif

!   bond compressing 

!  set moving unit vector 
!  x/r, y/r, z/r


    mv%vec(1) = x(1,mv%b) - x(1,mv%a)
    mv%vec(2) = x(2,mv%b) - x(2,mv%a)
    mv%vec(3) = x(3,mv%b) - x(3,mv%a)
! 归一化
    rr = mv%vec(1)*mv%vec(1) + mv%vec(2)*mv%vec(2) + mv%vec(3)*mv%vec(3)
    rr = dsqrt(rr)
    mv%vec(1) = mv%vec(1)/rr
    mv%vec(2) = mv%vec(2)/rr
    mv%vec(3) = mv%vec(3)/rr

    do k=-1,-nshorten,-1
        do i=1,mv%num
           x(:,mv%ind(i))= x(:,mv%ind(i)) - mv_dr*mv%vec(:)
        enddo

!  evolute filename
       filename = '                    '
       lent =len_trim(stretch_name)
       filename(1:lent) = stretch_name(1:lent)
       if(k>-10) then
          kk = 2
       elseif(k>-100 .and. k<-10)then
          kk = 3
       elseif(k<=-100)then
          kk = 4
       endif

       filename(lent+1:lent+kk) = int_to_str(k)
       filename(lent+kk+1:lent+kk+4) = '.inp'

       if(inptype=='nwchem')then
          call write_nwchem(filename)
       elseif(inptype=='gaussian')then
          call write_gaussian(filename)
       elseif(index(inptype,'xyz').NE.0)then
          call writexyz('mol_config')
       endif

    enddo

    if (comp_var) then
       !open(103,file='bonds.txt',status='old')
       !do i=1,nbond
       !   read(103,*) c1,c2,bond_value(i)
       !enddo
       call w_bond
       !close(103)
    endif

    do i=1,mv%num
       x(:,mv%ind(i))= x(:,mv%ind(i)) + nshorten*mv_dr*mv%vec(:)
    enddo

!   bond enlongation
    do k=1,nenlongation
        do i=1,mv%num
           x(:,mv%ind(i))= x(:,mv%ind(i)) + mv_dr*mv%vec(:)
        enddo

!  evolute filename
       filename = '                    '
       lent =len_trim(stretch_name)
       filename(1:lent) = stretch_name(1:lent)

       if(k<10) then
          kk =1
       elseif(k>=10 .and. k<100)then
          kk = 2
       elseif(k>=100)then
          kk = 3
       endif

       filename(lent+1:lent+kk) = int_to_str(k)
       !print *,filename
       filename(lent+kk+1:lent+kk+4) = '.inp'
       !print *,filename
       if(inptype=='nwchem')then
          call write_nwchem(filename)
       elseif(inptype=='gaussian')then
          call write_gaussian(filename)
       elseif(index(inptype,'xyz').NE.0)then
          call writexyz('mol_config')
       endif
    enddo

    do i=1,mv%num
       x(:,mv%ind(i))= x(:,mv%ind(i)) - nenlongation*mv_dr*mv%vec(:)
    enddo
    
    return
end


!


