!============================================================
!                    WRITED BY GUO F.                       =
!                                                           =
!       Institute of Atomic & Molecular Physics of Si-      =
!                chuan University of China.                 =
!                                                           =
!               email: gfeng8410@foxmail.com                =
!============================================================ 
!
!------------------------------------------------------------ 
!
! -----------------------------------------------------------------------
! string and number routines

! -----------------------------------------------------------------------
! returns TRUE if str1 matches 1st chars of str2, FALSE otherwise
! also returns m = loc of next char in str2
! could test if next char is a space or tab

logical function match(str1,str2,m)
      implicit none

      character*(*) str1,str2
      integer m

      match = .FALSE.
      m = len(str1) + 1
      if (len(str1).gt.len(str2)) return
      if (str1.eq.str2(1:len(str1))) match = .TRUE.

      return
end
!
!    When Error Accurred!
! 
subroutine err(str)
     implicit none

     character*(*) str

     Write(*,'(4x,A)')str

     stop
end

function int_to_str(ii)
   implicit none

   character(8) int_to_str
   integer ii
   int_to_str = '          '
   write(int_to_str,'(I8)') ii

   int_to_str = adjustl(int_to_str)
   int_to_str = trim(int_to_str)
   return
end

function str_add(str1,str2)
   implicit none

   character*(*) str1
   character*(*) str2
   !character*(*) str3
   character*(*) str_add
   integer l1,l2,i

   str_add = '                                        '
   str1 = adjustl(str1)
   str2 = adjustl(str2)
   l1 = len(trim(str1))
   l2 = len(trim(str2))
   !print *,l1,l2
   if(l1+l2>40) then
      print *,'str1 and str2 to long: in str_add function!'
      stop
   endif

   do i=1,l1
      str_add(i:i) = str1(i:i)
   enddo

   do i=1,l2
      str_add(l1+i:l1+i) = str2(i:i)
   enddo

   do i= l2+1,40
      str_add(i:i) = ' '
   enddo
   !str_add(1:l1) = adjustl(str1)
   !str_add = adjustl(str_add)
   !str_add=trim(str_add)
   !print *,str_add

   return
end

