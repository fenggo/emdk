! email: gfeng.alan@foxmail.com

subroutine move
    use mode,only: x,movex,movey,movez,npart
    implicit none

    integer     i

    do i=1,npart
       x(1,i) = x(1,i) + movex
       x(2,i) = x(2,i) + movey
       x(3,i) = x(3,i) + movez
    enddo

    return
end
