subroutine atm_typ_nm(mstart,mend)
    use modfim,only: molelist
    use modnei,only: mang,mdih,atable ! ,dihedral
    use mode,  only: types,ctype!,alive!,q !, a, b, c
    use mod_nw,only: atom_type!,latmtyp
    !use mod_grd,only: fix
    !use modctrl,only: run_type
    implicit none

    integer         i,jj
    integer         mstart, mend


    do jj=mstart, mend ! mlist(ii) molecule ii list end
       i = molelist(jj)
       !print *,'atom ',i
       if(ctype(types(i))=='C') then
           atom_type(i) = 'C_3'
       elseif(ctype(types(i))=='H') then
           atom_type(i) = 'H_'
       elseif(ctype(types(i))=='O') then
           atom_type(i) = 'O_2'
       elseif(ctype(types(i))=='N') then
           atom_type(i) = 'N_R'
       endif
    enddo

    return
end


