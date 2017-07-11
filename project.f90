!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                     email: gfeng.alan@hotmail.com                      =
!=========================================================================
!
!*************************************************************************
!**         Project the coordinate according lattice vectors            **
!*************************************************************************
!
!
subroutine proj_coord
use mode, only: npart, a,b,c,x,xd
implicit none

integer     i


!  应求解三元一次线性方程组

!  x(i) = xd1*a + xd2*b + xd3*c
!  反解出 xd1, xd2, xd3 即：相对坐标 xd(1-3,i)


do i=1,npart
     xd(3,i) = (a(1)*b(2)-a(2)*b(1))*(b(2)*x(3,i)-b(3)*x(2,i)) - (a(3)*b(2)-a(2)*b(3))*(b(2)*x(1,i)-b(1)*x(2,i))
     xd(3,i) = xd(3,i)/((a(1)*b(2)-a(2)*b(1))*(c(3)*b(2)-c(2)*b(3)) - (a(3)*b(2)-a(2)*b(3))*(c(1)*b(2)-c(2)*b(1)))

     xd(1,i) = b(2)*x(1,i)-b(1)*x(2,i) - (c(1)*b(2)-c(2)*b(1))*xd(3,i)
     xd(1,i) = xd(1,i)/(a(1)*b(2)-a(2)*b(1))

     xd(2,i) = x(2,i)- xd(1,i)*a(2) - xd(3,i)*c(2)
     xd(2,i) = xd(2,i)/b(2)
enddo


return
end
!
!*************************************************************************
!** Convert the Direct coordinate to Cartation according lattice vectors**
!*************************************************************************
!
subroutine cart_coord
use mode, only: npart,a,b,c,x,xd
implicit none

integer       i

do i=1,npart
   x(1,i)= xd(1,i)*a(1) + xd(2,i)*b(1) + xd(3,i)*c(1)
   x(2,i)= xd(1,i)*a(2) + xd(2,i)*b(2) + xd(3,i)*c(2)
   x(3,i)= xd(1,i)*a(3) + xd(2,i)*b(3) + xd(3,i)*c(3)
enddo

return
end
!
!
