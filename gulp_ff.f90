!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                      LiaoCheng University of China                     =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!
!-------------------------------------------------------------------------
!Write GULP input file and force  fieldfor Molecular Mechanic Simulations-
!-------------------------------------------------------------------------
!

subroutine gulp_ff
    !use modfim,only: n_mol
    use modnei,only: nbond,bond,n_ang,mang,mdih
    use mode,  only: x, a, b, c,npart,alive,q
    use mod_nw,only: atom_type
    !use mod_grd,only: fix
    use modctrl,only: run_type
    implicit none

    integer         i,iii
    !integer         nw_bonds,nw_angles,nw_dihes
    !integer         bondtype(10,20),ibondtype
    integer,allocatable::id(:)

    real(8),external:: torsion_angle

    allocate(id(npart))

    call  neighbor
    call  findmole

!  write gulp force filed, using NWchem output

!  Determining Angles

    n_ang = 0
    !print *,nbond

    call atm_typ

! ###### write input #######
     !print *,'where ?'
    open(31,file='inp-gulp',status='unknown')
   
    if(run_type==0) then
         write(31,'(A38)')'opti nosymmetry conp molmec #conjugate'
         write(31,*) !opti conv nosymmetry molecule 
         !write(31,'(A9)')'# for MD:'
    elseif(run_type==1) then
         write(31, '(A14)') 'md conv molmec'
         write(31, '(A12)') 'ensemble nvt'

         write(31, '(A23)') 'tau_thermostat  0.05 ps'
         write(31, '(A24)') 'temperature 298.000000 K'

         write(31, '(A20)') 'timestep        1 fs'
         write(31, '(A20)') 'production     20 ps'
         write(31, '(A20)') 'equilibration 0.5 ps'
         write(31, '(A17)') 'write          50'
         write(31, '(A17)') 'sample         50'
    elseif(run_type==2) then

         write(31, '(A14)') 'md conp molmec'
         write(31, '(A26)') 'integrator leapfrog verlet'
         write(31, '(A24)') 'ensemble npt 0.005 0.005'
         write(31, '(A17)') 'temperature 300 K'
         write(31, '(A17)') 'pressure 0.00 GPa'
         write(31, '(A22)') 'tau_barostat    0.1 ps'
         write(31, '(A22)') 'tau_thermostat  0.1 ps'

         write(31, '(A1)')'#'
         write(31, '(A1)')'#'


         write(31, '(A20)') 'timestep    0.001 ps'
         write(31, '(A20)') 'production     20 ps'
         write(31, '(A20)') 'equilibration   0 ps'
         write(31, '(A17)') 'write          50'
         write(31, '(A17)') 'sample         50'
    elseif(run_type==3) then
         write(31,'(A42)')'gradient nosymmetry conp molmec #conjugate'
    endif

    write(31,*)
    write(31,'(A5)') 'title'
    write(31,'(A16)') 'GULP calculation'
    write(31,'(A3)') 'end'
    write(31,*)

     !write(31,'(A10)') 'spacegroup'
     !write(31,'(A4)') ' P 1'
!library
    write(31,'(A16)') 'library mechanic'
    write(31,'(A21)') 'output movie xyz gulp'
    write(31,'(A19)') 'dump 50 restart.grs'
    write(31,*)
    write(31,'(A7)') 'vectors'
    write(31,'(3F16.8)') (a(i),i=1,3)
    write(31,'(3F16.8)') (b(i),i=1,3)
    write(31,'(3F16.8)') (c(i),i=1,3)
    write(31,*)
    write(31,'(A9)') 'cartesian'
    !!! determine atom type
    iii = 0
    !do ii=1,n_mol
       !if(ii>1)then      ! n molecule ii list start
         !n=mlist(ii-1) + 1
       !else
         !n=1
       !endif
       !print *,'molecule ',ii,n_mol
       !do jj=n, mlist(ii)! mlist(ii) molecule ii list end
       do i = 1,npart
          !i = molelist(jj)
          !print *,'atom ',i
           if(alive(i)/=0) then
               write(31,'(A5,A6,3F15.6,3F8.4,3I3)')atom_type(i),' core ',&
               & x(1,i),x(2,i),x(3,i),q(i),1.0,0.0,  1,1,1 ! 1=vary 0=fix
               iii = iii + 1
               id(i) = iii
           endif
       enddo
    !enddo

    do i=1,nbond
       write(31,'(A8,2I8)')'connect ',id(bond(1,i)),id(bond(2,i))
    enddo
!
    deallocate(id)
    close(31)

    return
end

