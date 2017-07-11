!                         WRITED BY GUO F.                 
!                                                                      
!               Institute of Atomic & Molecular Physics of Si-          
!                     chuan University of China.                       
!                                                                       
!                   email: gfeng.alan@hotmail.com                     
!
!
!-------------------- grid operations -------------------
! 

subroutine grid
    use mode,only:   npart !, x, types,alive
    use mod_grd, only: ngd,gpart
    implicit none

    integer                ig,ix,iline
    real(8)                  dd
    character(10)     tmp
    character(80)     line

     iline = 0
     rewind(5)

     !print *,'* i am here!'
     !open(115,file='grid.txt',status='old')

     do
        read(5,200,end=99,err=99)line

        iline= iline + 1
        line = trim(adjustl(line))
        if(len(line)==0)cycle
        !print *,line
        if(index(line,'%totgrid').NE.0)then
               read(line,*)tmp, ngd
               gpart = npart / ngd
               print *,gpart,'atoms in every grid'
        elseif(index(line,'%movgrid').NE.0)then
               read(line,*)tmp,ig,ix,dd
               call movgrid(ig,ix,dd)
               print *,'moving grid ...'
        elseif(index(line,'%delgrid').NE.0)then
               read(line,*)tmp,ig
               call delgrid(ig)
               print *,'deleting grid ...'
        elseif(index(line,'%fixgrid').NE.0)then
               read(line,*)tmp,ig
               call fixgrid(ig)
               print *,'fixing grid ...'
        endif
     enddo

99  continue

!   close(115)
   !print *,iline
200  format(A80)

   return
end
! 
!
!-------------------- grid operations -------------------
! 

subroutine movgrid(ig,ix,dist)
    use mode,only:  x
    use mod_grd, only: gpart
    implicit none

    integer        ig,gs,ge
    integer        ix,i
    real(8)           dist

    gs= 1+ (ig-1)*gpart
    ge= ig*gpart

    do i=gs,ge                ! atoms in grid
          x(ix,i)=x(ix,i) + dist
    enddo

    return
end
!
!-------------------- grid operations -------------------
! 

subroutine delgrid(ig)
    use mode,only:  alive
    use mod_grd, only: gpart
    implicit none

    integer        ig,gs,ge,i

    gs= 1+ (ig-1)*gpart
    ge= ig*gpart

    do i=gs,ge                ! atoms in grid
          alive(i) = 0
    enddo

    return
end
!
!-------------------- grid operations -------------------
! 

subroutine fixgrid(ig)
    use mode,only:  types!,npart
    use mod_grd, only: gpart,fix
    use modctrl,only: iwr
    implicit none

    integer        ig,gs,ge,i

    fix = 1

    gs= 1+ (ig-1)*gpart
    ge= ig*gpart

    do i=gs,ge         ! atoms in grid
          if(iwr==0) types(i) = types(i) + 4 !fixed
          fix(i) = 0
    enddo

    return
end

