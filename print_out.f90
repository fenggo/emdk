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
!
!-------------------------------------------------------------------------
!  ------             MAIN  FUNCTION            ------                   -
!-------------------------------------------------------------------------
!
!
    subroutine write_lp
     use mode,    only: a0,b0,c0,a,b,c,r1,r2,r3,alpha,beta,gamm,sc,mass,npart,types
    implicit none

    integer                   i
    real(8)                    cos_alpha,cos_beta,cos_gamma
    real(8) ::                 pi = 3.1415927
    real(8) ::                 mass_all
    real(8)                    vol, den, den_cm
    !integer                   typ(10)
    !integer                   id

    pi = acos(-1.0d0)
    !print *,pi
!   sc = super cell

    if (sc(1)>1.or.sc(2)>1.or.sc(3)>1) then
       print '(A30)','** Super cell boundary vector:'
       print '(A1)' , ' '
       print '(3F12.5)',(a0(i)*sc(1),i=1,3)
       print '(3F12.5)',(b0(i)*sc(2),i=1,3)
       print '(3F12.5)',(c0(i)*sc(3),i=1,3)    
    endif
    r1 = a0(1)*a0(1)+ a0(2)*a0(2)+ a0(3)*a0(3)
    r1 = dsqrt(r1)
    r2 = b0(1)*b0(1)+ b0(2)*b0(2)+ b0(3)*b0(3)
    r2 = dsqrt(r2)
    r3 = c0(1)*c0(1)+ c0(2)*c0(2)+ c0(3)*c0(3)
    r3 = dsqrt(r3)

    cos_alpha = (b0(1)*c0(1) + b0(2)*c0(2) + b0(3)*c0(3))/(r2*r3)
    cos_beta  = (c0(1)*a0(1) + c0(2)*a0(2) + a0(3)*c0(3))/(r3*r1)
    cos_gamma = (a0(1)*b0(1) + a0(2)*b0(2) + a0(3)*b0(3))/(r1*r2)

    alpha = dacos(cos_alpha)*180.0d0/pi
    beta  = dacos(cos_beta)*180.0d0/pi
    gamm = dacos(cos_gamma)*180.0d0/pi 

    r1= r1*sc(1)
    r2= r2*sc(2)
    r3= r3*sc(3)
    
    !print '(3(A4,F12.4))','a=',r1,'b=',r2,'c=',r3
    !print *,alpha,beta,gamma
    print '(A2)' , '**'

    print '(A39)','** Cell parameters (Angstroms/Degrees):'
    print '(A2)' , '**'
    !stop
    print *,'a=  ', r1 ,'alpha= ', alpha ! '(A4,F12.4,A7,F8.4)'
    print *,'b=  ', r2 ,'beta = ', beta
    print *,'c=  ', r3 ,'gamma=' , gamm
    !stop
    vol   = c(1)*b(2)*a(3) - c(1)*b(3)*a(2) + c(2)*b(3)*a(1) - c(2)*b(1)*a(3) &
             + c(3)*a(2)*b(1) - c(3)*b(2)*a(1)

    !typ = 0

    do i=1,npart
          mass_all = mass_all + mass(types(i))
          !print *,mass(types(i))
          !typ(types(i))=typ(types(i))+1
    enddo
    !stop
    vol = abs(vol)
    den   = mass_all/vol
    den_cm= den*10.0000/6.02253


    !do i=1,4
    !print '(A10,1X,I5,A3,1X,A9,35X)','Totally :',typ(i),ctype(i),'atoms...'
    !enddo
    print *,'*','The total volum of the simulation box:' !'(A1,1X,A37,30X)'
    print *,'*', vol, 'Ang^3'!'(A1,1X,F20.8,2X,A6,39X)'
    print *,'*','The total mass in the simulation box:' !'(A1,1X,A37,30X)'
    print *,'*',mass_all,'g/mole' !'(A1,1X,F20.8,2X,A6,39X)'
    print *,'*','The total density of the system: '!'(A1,1X,A31,36X)'
    print *,'*',den,'g/mole/A^3' !'(A1,1X,F20.8,2X,A10,35X)'
    print *,'*',den_cm,'g/cm^3' !'(A1,1X,F20.8,2X,A10,35X)'

!    print *,r2*cos(gamma),r2*sin(gamma)

    return
    end

