!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                     =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!  ------             WRITE topology file for Gromacs      ------        -
!-------------------------------------------------------------------------
!
subroutine write_uspex(filename)
    use mode,only: x,npart
    use mod_nw,only: atom_type
    !use modnei,only: nbond,bond,n_ang,angle,n_dih,dihedral
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
    call gulp_ff

    open(36,file=filename,status='unknown') 

    write(36,'(I7,3F12.4)')npart
    write(36,'(A10)')'mol_charge'
    do i=1,npart
        write(36,'(A7,F10.5,3F12.4)') atom_type(i),0.0,x(1,i),x(2,i),x(3,i)
    enddo

    close(36)
    return
end

