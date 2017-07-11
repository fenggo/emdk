set terminal post eps color enhanced solid linewidth 3
set out "no.eps"
set xlabel "bondNumb"
set ylabel "bondlength"
plot "nw_no.txt" using 1:2 with linespoints ps 2 pt 4 title "NO"
