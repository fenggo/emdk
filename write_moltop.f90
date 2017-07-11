!

!-------------------------------------------------------------------------

!       Write out molecular topology-table for dl_field and dy_poly      -

!-------------------------------------------------------------------------

!

     subroutine write_moltop
     use mode, only: npart,types,ctype,x,q,r1,r2,r3,alpha,beta,gamm
     use modnei,only: matom, atable
     use modfim,only: molecule,n_mol,molelist,mlist
     implicit none


     integer i,j,natom
     !integer,allocatable::       str_num(:)  

!    e.g. n.u. of CHON <==> mol_table(2,2,2,2) because 
!                NO can be mol_table(1,1,2,2)

     !real(8)                     xold(3)   ! for Larger moleculars split by PBC
     character(20)            :: mol_name
     character(4)                str1
     character(4),allocatable :: atm_nm(:)
     character(20)            :: atm_na = '                    '
     character(4)                ntab(8)
     character(3)                mol,atm,res

     integer                     leng,lenn,lens
     integer                     nty(5)
     integer                     k,kk,m,jj,iatom,l
     integer                     nof1,nof2,nof3,nof4,nof5
     real(8)         ::          pi = 3.1415927
!

     allocate(atm_nm(npart))

     open(114,file='moltop.txt',status='unknown')
     open(111,file='moltop.pdb',status='unknown')

     write(114,'(A20)')'Molecular topology'
     write(111,'(A6,3F9.3,3F7.2,A11,4X)')'CRYST1',r1,r2,r3,alpha*180.0d0/pi,beta*180.0d0/pi,gamm*180.0d0/pi,'P1'

     call neighbor
     call findmole
     !print *,'Where is the problem?'

       natom = 0
       do i=1, n_mol    ! cycle from molecules
          if(i>1)then
             jj=mlist(i-1) + 1
          else
             jj=1
          endif

          nof1= 0
          nof2= 0
          nof3= 0
          nof4= 0
          nof5= 0

          do j=jj, mlist(i)
             atm_na = '               '
             iatom = molelist(j)
             if    (types(iatom)==1)then
                  nof1=nof1 + 1
                  lenn = len_trim(ctype(1)) 
                  atm_na(1:lenn)=ctype(1)
                  write(str1,'(I4)')nof1
                  str1=adjustl(str1)
                  lens=len_trim(str1)+lenn
                  atm_na(lenn+1:lens)=str1
             elseif(types(iatom)==2)then
                  nof2=nof2 + 1
                  lenn = len_trim(ctype(2)) 
                  atm_na(1:lenn)=ctype(2)
                  write(str1,'(I4)')nof2
                  str1=adjustl(str1)
                  lens=len_trim(str1)+lenn
                  atm_na(lenn+1:lens)=str1
             elseif(types(iatom)==3)then
                  nof3=nof3 + 1
                  lenn = len_trim(ctype(3)) 
                  atm_na(1:lenn)=ctype(3)
                  write(str1,'(I4)')nof3
                  str1=adjustl(str1)
                  lens=len_trim(str1)+lenn
                  atm_na(lenn+1:lens)=str1
             elseif(types(iatom)==4)then
                  nof4=nof4 + 1
                  lenn = len_trim(ctype(4)) 
                  atm_na(1:lenn)=ctype(4)
                  write(str1,'(I4)')nof4
                  str1=adjustl(str1)
                  lens=len_trim(str1)+lenn
                  atm_na(lenn+1:lens)=str1
             elseif(types(iatom)==5)then
                  nof5=nof5 + 1
                  lenn = len_trim(ctype(5)) 
                  atm_na(1:lenn)=ctype(5)
                  write(str1,'(I4)')nof5
                  str1=adjustl(str1)
                  lens=len_trim(str1)+lenn
                  atm_na(lenn+1:lens)=str1
             endif

             atm_nm(iatom) =adjustl(atm_na)
             !print *,atm_nm(iatom),atm_na
          enddo

!  molecular structure

          nty(1)= nof1
          nty(2)= nof2
          nty(3)= nof3
          nty(4)= nof4  
          nty(5)= nof5   
          !nty(:)= str_name(:)
          mol_name='                    '

          do j=1,5
             if(nty(j)>0)then
                  leng = len_trim(mol_name)
                  lenn = len_trim(ctype(j)) + leng
                  mol_name(leng+1:lenn)=ctype(j)
                  if(nty(j)>1)then
                         write(str1,'(I4)')nty(j)
                         str1=adjustl(str1)
                         lens=len_trim(str1)+lenn
                         mol_name(lenn+1:lens)=str1
                  endif
             endif
          enddo  

          write(114,'(A8,1X,A20,1X,I8,F8.5)')'MOLECULE',mol_name,nof1+nof2+nof3+nof4+nof5,0.0 
          atm_na = '               '

          do j=jj, mlist(i)
             iatom = molelist(j)
             if(types(iatom)==1) then
                   atm_na= 'C_aliphatic1 '
             elseif(types(iatom)==2) then
                   atm_na= 'H_aliphatic'
             elseif(types(iatom)==3) then
                   atm_na= 'O_monoxide'
             elseif(types(iatom)==4) then
                   atm_na= 'N_amide '
             endif
             write(114,'(A5,A20,F8.5)')adjustl(atm_nm(iatom)), adjustl(atm_na), 0.00
          enddo

          do j=jj, mlist(i)
             iatom = molelist(j)
             natom = natom + 1
             kk = atable(1,iatom)
             do k=1,kk
                ntab(k) = adjustl(atm_nm(atable(k+1,iatom)))
             enddo 

             write(114,*)'CONNECT ',atm_nm(iatom), ' > ',(ntab(m),m=1,kk)
             atm_na = '                    '
             if(types(iatom)==1) then
                   atm= 'C1 '
             elseif(types(iatom)==2) then
                   atm= 'H1 '
             elseif(types(iatom)==3) then
                   atm= 'OM '
             elseif(types(iatom)==4) then
                   atm= 'N1 '
             endif
             if(mol_name=='C6H6O12N12') then
                   mol= 'cl2'
                   res= 'EM2'
             elseif(mol_name=='C4H8O8N8') then
                   mol= 'hmx'
                   res= 'EM3'
             endif
             write(111,'(A6,I5,1X,A4,A4,2X,I4,4X,3F8.3,2F6.2,5X,A3,2X,A2)')'ATOM  ',natom,atm,&
                     res,molecule(iatom),(x(l,iatom),l=1,3),1.00,q(i),mol,ctype(types(iatom)) ! HMX-CL20
          enddo
          write(114,'(A12)')'END MOLECULE'

       enddo  ! FindMole Main loop
      
     write(111,'(A3)')'END'

     close(114)
     close(111)
     deallocate(atm_nm)
!
     return
     end

