!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@foxmail.com                     =
!=========================================================================
!
!

subroutine atm_typ_cl20(mstart,mend)
use modfim,only: molelist
use modnei,only: mang,mdih,atable
use mode,  only: types,ctype
use mod_nw,only: atom_type!,latmtyp
!use mod_grd,only: fix
!use modctrl,only: run_type
implicit none

    integer         mstart, mend
    integer         i,j,k,l,m,n,jj
    integer         iatom,jatom,katom,matom,latom,llatom,mmatom
    integer         n_c,n_h,n_n,n_o
    integer      :: dihe_vertical = 0
    integer         dihe_parallel
    integer         dihe_middle
    integer         dihed(2)
    character(7)         tc
    real(8)         dih,dmax,dmin
    !real(8)         radi!,rsq1,rsq2


    real(8),external:: torsion_angle
    !logical           :: dihe_vertical = .false.

    

!@@@@@@@@@@@@@@@@   topology analysis for molecule CL20   @@@@@@@@@@@@@@@@
 

    dmax = 0.0
    dmin = 180.0


    do jj=mstart, mend
       i = molelist(jj)
       if(ctype(types(i))=='H') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          if(n_c==1) atom_type(i) = 'H0'  
       elseif(ctype(types(i))=='C') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          if(n_c==1 .and. n_n==2 .and. n_h==1) atom_type(i) = 'C010'  ! General
       elseif(ctype(types(i))=='N') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  

          if(n_c==2.and.n_n==1) atom_type(i) = 'N010'  ! N sp3 bonded to C, general 
!         N1c3 ! N sp3 bonded to C
          if(n_o==2.and.atable(1,i)==3) then  !######## N in nitro-group
             atom_type(i) = 'N015' ! N1co  ! N in nitro-group general 
          endif                               !########## end analysis N in nitro-group
       elseif(ctype(types(i))=='O') then
          n_c = 0
          n_h = 0
          n_n = 0
          n_o = 0
          do j=2,atable(1,i)+1
             n = atable(j,i)
             tc = ctype(types(n))
             if (tc=='H') n_h = n_h + 1
             if (tc=='C') n_c = n_c + 1
             if (tc=='N') n_n = n_n + 1
             if (tc=='O') n_o = n_o + 1
          enddo  
          if(n_n==1.and.atable(1,i)==1) atom_type(i) = 'O2' ! sp2 Oxygen General
          !!Oxygen in Nitro-group
       endif
    enddo

    ! Find Cc30
    print '(A2)','**'
    print '(A17)','** diheral angle:'
    print '(A2)','**'

    do jj=mstart, mend
       iatom = molelist(jj)
       if ( ctype(types(iatom))=='C') then
          do j=2,atable(1,iatom)+1
             jatom = atable(j,iatom)
             ! dihedral angele
             if (ctype(types(jatom))=='C') then
                dihe_vertical = 0
                dihe_parallel = 0
                dihe_middle = 0
                do l = 2,atable(1,jatom)+1
                   latom = atable(l,jatom)
                   if (ctype(types(latom))=='N') then
                      do m = 2,atable(1,latom)+1
                         matom = atable(m,latom)
                         if (ctype(types(matom))=='N') then
                             dih = torsion_angle(iatom,jatom,latom,matom)
                             !print '(F11.5,4I6)',dih,iatom,jatom,latom,matom
                             dih = abs(dih)
                             !dih = dih - 90.0           
                             if(dih<100.0 .and. dih>60.0) then
                                 dihe_vertical = dihe_vertical + 1
                                 atom_type(matom) ='N016'     ! N in nitro: vertical
                                 atom_type(latom) ='N011'    ! N bond to C: vertical
                                 !print *, matom,latom,atom_type(matom)
                             elseif(dih>165.0) then
                                 dihe_parallel = dihe_parallel + 1
                                 atom_type(matom) ='N014'  !N1co1  N in nitro: parallel
                                 !atom_type(latom) ='Nc31'
                                 !print *,matom,latom,atom_type(latom)
                                 !print *, matom,latom,atom_type(matom)
                             elseif(dih<160.0 .and. dih>100.0) then
                                 dihe_middle = dihe_middle + 1
                                 !atom_type(matom) ='Nco1'
                                 !atom_type(latom) ='Nc31'
                                 !print *, matom,latom,atom_type(matom)
                             endif
                         endif
                      enddo
                   endif
                enddo 

                do l = 2,atable(1,iatom)+1
                   llatom = atable(l,iatom)
                   if (ctype(types(llatom))=='N') then
                      do m = 2,atable(1,llatom)+1
                         mmatom = atable(m,llatom)
                         if (ctype(types(mmatom))=='N') then
                             dih = torsion_angle(jatom,iatom,llatom,mmatom)
                             !print '(F11.5,4I6)',dih,jatom,iatom,llatom,mmatom
                             dih = abs(dih)          
                             if(dih<100.0 .and. dih>60.0) then
                                 dihe_vertical = dihe_vertical + 1
                                 !print *, matom,latom,atom_type(matom)
                             elseif(dih>165.0) then
                                 dihe_parallel = dihe_parallel + 1
                                 !print *, matom,latom,atom_type(matom)
                             elseif(dih<160.0 .and. dih>100.0) then
                                 dihe_middle = dihe_middle + 1
                                 !print *, matom,latom,atom_type(matom)
                             endif
                         endif
                      enddo
                   endif
                enddo 

                ! C-C bond have two vertical dihedral angele
                if (dihe_vertical == 2) then !C010  ! General
                   atom_type(jatom) ='C011'  !C1c30 with two vertical
                   atom_type(iatom) ='C012'  !C1c32
                   !print *,dihe_vertical
                endif
                if (dihe_parallel == 2) then
                   atom_type(jatom) ='C013' !C1c31 with two parallel
                   atom_type(iatom) ='C013'
                endif   
                if (dihe_middle == 4) then
                   atom_type(jatom) ='C014' !C1c31 with four midle
                   atom_type(iatom) ='C014'
                endif   
                !print *,dihe_parallel
             endif
          enddo 
       endif
    enddo

!______________________________________________________________________
!
    do jj=mstart, mend
       iatom = molelist(jj)
       if ( atom_type(iatom)=='C012') then !C012  C with two vertical N-NO2
          do j=2,atable(1,iatom)+1
             jatom = atable(j,iatom)
             ! dihedral angle
             if (ctype(types(jatom))=='N') then
                 atom_type(jatom) = 'N017' !N1c32 vertical N-NO2 (N with-C)
                 do l = 2,atable(1,jatom)+1
                    latom = atable(l,jatom)
                    if (ctype(types(latom))=='N') then
                       atom_type(latom) ='N012' !N of NO2
                    endif
                 enddo
             endif
          enddo
       endif

       if ( atom_type(iatom)=='C014') then !C1c31 with four midle
          !print *,iatom
          do j=2,atable(1,iatom)+1
             jatom = atable(j,iatom)
             ! dihedral angele
             if (ctype(types(jatom))=='N' .and. atom_type(jatom)/='N011') then
                 !print *,atom_type(jatom)
                 atom_type(jatom) ='N018'
                 do l = 2,atable(1,jatom)+1
                    latom = atable(l,jatom)
                    if (ctype(types(latom))=='N') then
                       atom_type(latom) ='N013' !N of NO2
                    endif
                 enddo
             endif
          enddo
       endif
    enddo

    !  rename O atom
    do jj=mstart, mend
       iatom = molelist(jj)
       if ( ctype(types(iatom))=='O') then !C012  C with two vertical N-NO2
          do j=2,atable(1,iatom)+1
             jatom = atable(j,iatom)
             if ( ctype(types(jatom))=='N') then
                do k = 2,atable(1,jatom)+1
                   katom = atable(k,jatom)
                   if ( ctype(types(katom))=='N') then
                      dih = 0.0
                      dihed = 0
                      m = 0
                      do l = 2,atable(1,katom)+1
                         latom = atable(l,katom)
                         if ( ctype(types(latom))=='C') then 
                         	m = m+1
                            dih = torsion_angle(iatom,jatom,katom,latom)
                            if (dih<90.0) then
                               read(atom_type(latom),'(2X,I2)') dihed(1)
                            else
                               read(atom_type(latom),'(2X,I2)') dihed(2)
                            endif
                         endif
                      enddo
                   endif
                   if (dihed(1)<dihed(2)) then
                   	  atom_type(iatom) = 'O21  ' !# O21
                   else
                      atom_type(iatom) = 'O21  ' !# O22
                   endif
                enddo
             endif
          enddo
       endif
    enddo
    return 
end 
