!=========================================================================
!                            WRITED BY GUO F.                            =
!                                                                        =
!                     LiaoCheng University of China                      =
!                                                                        =
!                      email: gfeng.alan@foxmail.com                     =
!=========================================================================
!
!

subroutine atm_typ_hmx(mstart,mend)
    use modfim,only: molelist
    use modnei,only: mang,mdih,atable ! ,dihedral
    use mode,  only: types,ctype!,alive!,q !, a, b, c
    use mod_nw,only: atom_type!,latmtyp
    !use mod_grd,only: fix
    !use modctrl,only: run_type
    implicit none

    integer         i,j,k,l,m,n,jj,o,iii,jjj
    integer         iatom,jatom,katom,latom,jjatom,kkatom,llatom
    integer         n_c,n_h,n_n,n_o
    integer         mstart, mend
    character(7)         tc
    real(8)         dih,dmax,dmin
    real(8),external:: torsion_angle
    logical           :: dihe_vertical = .false.
    

!@@@@@@@@@@@@@@@@       topology analysis      @@@@@@@@@@@@@@@@
    !!! determine atom type
    !print *,n_mol
    !do ii=1,n_mol
       !if(ii>1)then      ! n molecule ii list start
         !n=mlist(ii-1) + 1
       !else
         !n=1
       !endif
       !print *,'molecule ',ii
       do jj=mstart, mend ! mlist(ii) molecule ii list end
          i = molelist(jj)
          !print *,'atom ',i
          if(ctype(types(i))=='C') then
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
             if(n_h==2) atom_type(i) = 'C3' ! Carbon with two H sp3
          endif
       enddo
    !enddo

    dmax = 0.0
    dmin = 180.0

    !do ii=1,n_mol
       !if(ii>1)then
         !n=mlist(ii-1) + 1
       !else
         !n=1
       !endif
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
             if(n_c==1) atom_type(i) = 'H'
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
             dihe_vertical = .false.
             if(n_c==2.and.n_n==1) atom_type(i) = 'N3'  ! N sp3

             if(n_o==2.and.atable(1,i)==3) then  !######## N in nitro-group
                 iatom = i ! N of nitro-group
                 do k = 2,atable(1,iatom)+1    ! find atoms bond to iatom
                    l = atable(k,iatom)
                    if (ctype(types(l))/='O') then
                    ! iatom  jatom katom and latom 
                         jatom = l         ! jatom  have found ...
                      !print *,'jatom',jatom
                    endif
                 enddo
                 do k = 2,atable(1,jatom)+1     !! find katom at jatom
                    l = atable(k,jatom) 
                    if( ctype(types(l))/='H' .and. l/=iatom) then !! katom  have found ...
                       katom = l
                       !print *,'katom',katom    !! find latom at katom
                       do o = 2,atable(1,katom)+1
                          m = atable(o,katom) 
                          if(ctype(types(m))/='H' .and. m/=jatom)then !! latom have found
                              latom = m
                              dih = torsion_angle(iatom,jatom,katom,latom)
                              !print *,dih,iatom,jatom,katom,latom
                              dih = abs(dih)
                              !dih = dih - 90.0

                              if(dih<105.0 ) then
                                  dihe_vertical = .true.

                                  if (dmax < dih) dmax = dih
                                  if (dmin > dih) dmin = dih
                                  !print *,dih,iatom,jatom,katom,latom
                                  !if(ctype(types(katom))=='C') atom_type(katom) = 'C30' ! for HMX  
                              endif
                          endif
                       enddo
                    endif
                 enddo

                 if(dihe_vertical) then
                     atom_type(i) = 'N32' ! N sp3 Nitro-group N_R1
                 else
                     atom_type(i) = 'N31' ! N sp3 Nitro-group N_R2
                 endif
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
             if(n_n==1.and.atable(1,i)==1) atom_type(i) = 'O2' ! sp2 Oxygen
             !!Oxygen in Nitro-group
          endif
       enddo
    !enddo

    ! Find C30
    !do ii=1,n_mol
       !if(ii>1)then
         !n=mlist(ii-1) + 1
       !else
         !n=1
       !endif
       do jj=mstart, mend
          i = molelist(jj)
          if(atom_type(i)=='N32') then
             jjatom = i
             !print *,'N32',jjatom
             do j = 2,atable(1,i)+1
                iatom = atable(j,i)
                if (atom_type(iatom)=='N3') then
                   kkatom = iatom
                   !print *,'N3',kkatom
                   do iii = 2,atable(1,kkatom)+1
                      if (atom_type(atable(iii,kkatom))=='C3') then
                          katom = atable(iii,kkatom)
                          !print *,'C3',katom
                          do jjj = 2,atable(1,katom)+1
                             if (atom_type(atable(jjj,katom))=='N3' .and. atable(jjj,katom) /= kkatom) then
                                llatom = atable(jjj,katom)
                                !print *,'N3',llatom
                                !do kkk = 2,atable(1,llatom)+1
                                   !if (atom_type(atable(kkk,llatom))=='N31') then
                                !#####!iiatom = atable(kkk,llatom)
                                      !print *,'N31',iiatom

                                      dih = torsion_angle(jjatom,kkatom,katom,llatom)
                                      !print *,dih, jjatom,kkatom,katom,llatom
                                      !print *,0.5*(dmax+dmin)
                                      if (dih > 0.5*(dmax+dmin)) atom_type(katom)='C30'
                                   !endif
                                !enddo
                             endif
                          enddo
                      endif
                   enddo
                endif
             enddo
          endif
       enddo
    !enddo
!@@@@@@@@@@@@@@@@       topology analysis      @@@@@@@@@@@@@@@@

    return
end