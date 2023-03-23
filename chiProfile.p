
#set terminal png size 800,500 enhanced font "Helvetica,16"
#set output 'outputT.png' 
set xlabel "Z [-]"
set ylabel "Chi [1/s]"
set xrange [0:1]
set yrange [0:50]

set grid
set key spacing 1 font "Helvetica, 14"

plot "chiProfile.dat" using 1:2  with line lt -1 lw 2    








