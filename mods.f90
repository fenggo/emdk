!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!-------------------------------------------------------------------------
!  ------             MODULE                    ------                   -
!-------------------------------------------------------------------------
!

module mode
implicit none
!
     integer, save::              npart  ! total patical number
     integer, save::              cpart  ! patical number in cell
     integer, save::              sc(3)  ! supper cell N*N*N
     integer, save::              nelem  ! Inputted element number

     real(8), save::              cel(3)
! lattice parameter
! a, b, c --> Lattice vector
     real(8), save::              a0(3),a(3)
     real(8), save::              b0(3),b(3)
     real(8), save::              c0(3),c(3)
     real(8)      ::              r1,r2,r3    ! lattice parameter, the lengh of a,b,c
     real(8)      ::              alpha,beta,gamm
     real(8), save::              vel

! cell movation in three directions
     real(8), save::              movex= 0.0
     real(8), save::              movey= 0.0
     real(8), save::              movez= 0.0
     real(8), save::              mass(10)
     real(8), save::              shockv   ! shock velocity for a symmetric impact

     integer, allocatable, save:: types(:)
     integer, allocatable, save:: alive(:)       ! 0 = delete, 1= no change
     integer,              save:: ialive = 0
     real(8), allocatable, save:: x(:,:)
     real(8), allocatable, save:: v(:,:)
     real(8), allocatable, save:: q(:)
     real(8), allocatable, save:: xd(:,:) !Direct coord: projected coordinate according a, b, c 
!
     character(2),         save:: ctype(100)
!
     character(40),        save:: struc_name
     !data ctype /'  ','  ','  ','  ','  '/
     data sc    /1,1,1/
!

end module
!
!-------------------------------------------------------------------------
!                                                                        -
!-------------------------------------------------------------------------

module   modnei
implicit none
!
    integer,  parameter:: mbond=50000
    integer,  parameter:: mang=20000
    integer,  parameter:: mdih=20000
    integer,  parameter:: matom=20
    integer,  save     :: nbond
    integer,  save     :: n_ang
    integer,  save     :: n_dih
!
    integer,  allocatable, save:: atable(:,:)   !  atom-table within the r-cut off
    integer,  save:: bond(2,mbond)
    integer,  save:: angle(3,mang)
    integer,  save:: dihedral(4,mdih)
!
    real(8),  save:: rcut(4,4) ! cut-off for electronic density
    real(8),  save:: rcutsq(4,4)

    real(8),  allocatable, save:: bond_value(:)
    real(8),  allocatable, save:: ang_value(:)
    real(8),  allocatable, save:: dih_value(:)
!
!   contribution of neigbors's density
!   real(8),  allocatable, save:: nlist(mlist) ! neighbor-list
!
    data rcut /1.7892,1.3266,1.7200,2.1200,&
               1.3266,0.8640,1.2574,1.2267,&
               1.7200,1.2574,1.6509,1.6201,&
               2.1200,1.2267,1.6201,1.5894/
!
!

end module

!-------------------------------------------------------------------------
!                                                                        -
!-------------------------------------------------------------------------

module modctrl
implicit none
!
     character(80),        save:: cfile = 'coord.xyz'
     integer,              save:: iwr    !输出类型
     integer,              save:: ilt        = 0
     integer,              save:: run_type   = 0
     real(8),              save:: izoom(3) = (/1.0,1.0,1.0/)
     logical,              save:: fzoom = .false.
     logical,              save:: fremove = .false.
     logical,              save:: ldipole = .false.
     logical,              save:: center=.false.
     logical,              save:: isc=.true.    !whether creat super-cell
     logical,              save:: irec=.false.
     logical,              save:: lmole= .false.
     logical,              save:: istretch = .false. ! stretching the bond used in parameter fitting 
     logical,              save:: istretchfreebond = .false. ! stretching the bond used in parameter fitting 
!     logical,              save:: ipty=.false.
!     logical,              save:: iptz=.false.
!    iptx,ipty,iptz: 是否将分裂的原子再放在一块
     logical,              save:: lmap = .false. ! whether mapping the lattice indices, used by Lammps fix-phonon
     logical,              save:: lcut = .false. ! develop the slap into two part
!

end module

!-------------------------------------------------------------------------
!                                                                        -
!-------------------------------------------------------------------------

module modfim
implicit none

     integer,              save:: n_mol = 0
     integer,allocatable,  save:: cflag(:)
     integer,allocatable,  save:: molecule(:)
!
!   mlist(), molelist() : molecule lookup table!
!   mlist() is the 'table of contents'
!
!   cflag signs which atom has already considered
!   molecue, mlist, and molelist establelish a list
!   that we can find atom that belong to one molecule
!
!   n_mol    :  Total molecule number
!   molecule :  Atom i belong to molecue 'molecule(i)'
!
!   Atoms in molecue i are stored in molelist() from-
!         mlist(i-1)+1 to mlist(i)
!
     integer,allocatable,  save:: molelist(:)
     integer,allocatable,  save:: mlist(:)

     end module

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

module modgr
implicit none
!
     integer,      save:: ntype(4)
     integer,      save:: icom1 = 1
     integer,      save:: icom2 = 1
     real(8),      save:: gcut = 10.0
!
     end module
!
!-------------------------------------------------------------------------
!

module mod_stretch
    implicit none

    !integer   ::   nstretch = 1
    integer        nenlongation
    integer        nshorten
    real(8)        mv_dr
  
    public :: move

    type  move

!   unit vector in stretching
        integer     a,b       ! index of atoms of the two bond end
!       real(8)     dr        ! length of the vector: move vector by dr*vec 
        real(8)     vec(3)   ! unit vector in the bond direction fix-end ---> stretch-end
        integer     num      ! number of atoms to be moved
        integer     ind(100)  ! atom index to be stretched
! Max number of atoms to be moved
    end type move

    type(move)::    mv
    character(30)   stretch_name
    character(6)    species1, species2
    logical   ::    comp_var = .false.
end module

!
!-------------------------------------------------------------------------
!

module lattice
    implicit none

! the lattice indices
! lattice indices: used in FixPhonon coupled in lammps!

    !integer                  natoms= cparts    
    integer,allocatable ::   Natom(:) ! the nth atom in unit-cell
    integer,allocatable ::   LatticeIndex(:,:)

end module
!
!-------------------------------------------------------------------------
!

module mod_def
    implicit none

!   create circular defects 
    integer            ::   idef = 0 ! # of defects
    real,allocatable   ::   Defect(:,:)

end module
!
!-------------------------------------------------------------------------
!

module mod_grd
    implicit none
!
    integer,  save:: ngd, gpart
    integer, allocatable ::  fix(:) ! fixed grid
    logical              ::  grid_op = .false.
!

end module mod_grd
!

module mod_nw
    implicit none

    integer   ,     save::  nw_bonds,nw_angles,nw_dihes
    integer,allocatable::  nw_bond(:,:),nw_angle(:,:),nw_dihe(:,:)

    real(8),allocatable::  nw_bond_val(:)
    real(8),allocatable::  nw_angle_val(:)
    real(8),allocatable::  nw_dihe_val(:)

    character(5),allocatable:: atom_type(:) 
    character(10)       ::  basis = '6-311G**'    
    logical                latmtyp   
end module mod_nw

module mod_rotate
    implicit none
!
    integer,  save:: nrotate, atom1,atom2,atom3,num_rot
    integer          rota(20)
!
    !real(8)          dih(20)
    real(8)          drot

    logical     ::   irotate = .false.
    logical     ::   iswing  = .false.

    character(8)        inptype
end module mod_rotate

module mod_mmfit
    implicit none
!
    integer,  save:: task

end module mod_mmfit
