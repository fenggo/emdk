!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@hotmail.com                      =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!------               zoom coordinate acoordining volumn            ------
!-------------------------------------------------------------------------
!
    subroutine zoom
    use mode,   only: a,b,c,r1,r2,r3
    use modctrl,only:izoom
    implicit none

    integer          i

    call proj_coord

    do i=1,3
       a(i)= a(i)*izoom(1)
       b(i)= b(i)*izoom(2)
       c(i)= c(i)*izoom(3)
    enddo

    r1 = r1*izoom(1)
    r2 = r2*izoom(2)
    r3 = r3*izoom(3)

    call cart_coord
    !print *,'debug !!'
    return
    end

