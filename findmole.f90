!=========================================================================
!                           WRITED BY GUO F.                             =
!                                                                        =
!              Institute of Atomic & Molecular Physics of Si-            =
!                       chuan University of China.                       =
!                                                                        =
!                      email: gfeng.alan@gmail.com                       =
!=========================================================================
!
!
!-------------------------------------------------------------------------
!  ------             FIND MOLECULE             ------                   -
!-------------------------------------------------------------------------
!
!
     subroutine findmole
     use mode,  only: npart
     use modfim,only: molelist,mlist,molecule,cflag,n_mol
     use modnei,only: matom
     implicit none    
!
!
     integer        i
     integer        ii
     integer        j
     !integer        jj
     integer        cc
     !integer        l
!
!
     n_mol=1
!
!   cflag signs which atom has already considered
!   molecue, mlist, and molelist establelish a list
!   that we can find atom that belong to one molecule
!   n_mol    :  Total molecule number
!   molecule :  Atom i belong to molecue 'molecule(i)'
!   Atoms in molecue i are stored in molelist from-
!         mlist(i-1)+1 to mlist(i)
!
     do i=1,npart
        cflag(i)=1
     enddo
! 
!
     do i=1, npart
!
       ii=i
       if(cflag(ii)/=0)then
          cflag(ii)=0
          call trackbond(ii)
          n_mol=n_mol+1
       endif
!
     enddo
!
     n_mol=n_mol-1
     cc=0
!     print 100,'Total Molecule Number:',n_mol
!
      do i=1, n_mol
          do j=1, npart
             if(molecule(j)==i)then
                cc=cc+1
                molelist(cc)=j
             endif
          enddo
          mlist(i)=cc
     enddo
!
!
     return
     end
!
!-------------------------------------------------------------------------
!  ------             RECURSIVE TRACBOND             ------              -
!-------------------------------------------------------------------------
!
   
    recursive subroutine trackbond(m)
    use modnei, only: atable,matom
    use modfim, only: n_mol,cflag,molecule
    implicit none
!
    integer      mmm(matom)
    integer      m, mm,i,tm
!
!
!    if(cflag(m)/=0)return
!
    mm=m
    do i=1,matom
       mmm(i)=atable(i+1,mm)
    enddo
!
    molecule(mm)=n_mol
    cflag(mm)=0
!
    do i=1,matom
!
       if(mmm(i)==0)exit
       if(cflag(mmm(i))/=0)then
           cflag(mmm(i))=0
           molecule(mmm(i))=n_mol
           tm=mmm(i)
           call trackbond(tm)
      endif
!
    enddo
!
    return
    end
!
!-------------------------------------------------------------------------
!

