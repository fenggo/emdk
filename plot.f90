!@@@@@@@@@@@@@@@@@@@@@@ plot function @@@@@@@@@@@@@@@@@@@@@@
subroutine plot(fil,x_lable,y_lable,title,fileplot)
    implicit none
    !integer       i,ii
    character(*) x_lable,y_lable,title,fil,fileplot
!   using gnuplot plot datas
!   asuming you have gnuplot in you system

    !ii = size(x)
    !open(222,file='data.txt',status='unknown')
    !do i=1,ii
       !write(222,*)x(i),y(i)
    !enddo
    !close(222)
    
    open(221,file='p.plt',status='unknown')
    write(221,'(A54)')'set terminal post eps color enhanced solid linewidth 3'
    write(221,'(A9,A,A1)')'set out "',fileplot,'"'
    write(221,'(A12,A,A1)')'set xlabel "',x_lable,'"'
    write(221,'(A12,A,A1)')'set ylabel "',y_lable,'"'
    !write(221,'(A52,A16)')'plot  "data.txt"  using 1:2 with linespoints ps 2 title ',title
    write(221,'(A6,A,A46,A,A1)')'plot "',fil,'" using 1:2 with linespoints ps 2 pt 4 title "',title,'"'
    close(221)
    call system('gnuplot p.plt')
    return
end

