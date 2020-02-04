set term pdf
set style arrow 1 head back filled linetype 1 linecolor rgb "red"
set style arrow 2 linetype 5 dt 2 linecolor rgb "blue"

set key left

set key font "Times, 18"
#set key at 11,22

set xrange [0:100]

set output "cumulative_distance.pdf"

set xlabel "{/Times iterations}" font ",22"
set ylabel "{/Times cumulative distance (m)}" font ",22"
set grid
plot 'loc.dat' using 17 title '{/Times O-IST}'  with lp pt 6 pi 5 ps 0.8 dt "_" lw 1.5  lc rgb "#00CC66",\
'' using 18 title '{/Times O-DR}' with lp pt 12 pi 5 ps 0.8 dt ".._" lw 1.5    lc rgb "#FF3300",\
'' using 19 title '{/Times O-DISTA}' with lp pt 8 pi 5 ps 0.8 dt "._" lw 1.5  lc rgb "#3333FF",\


set output "distance.pdf"
set xlabel "{/Times iterations}" font ",22"
set ylabel "{/Times  distance (m)}" font ",22"
set grid
set grid
plot 'loc.dat' using 20 title '{/Times O-IST}'  with lp pt 6 pi 5 ps 0.8 dt "_" lw 1.5  lc rgb "#00CC66",\
'' using 21 title '{/Times O-DR}' with lp pt 12 pi 5 ps 0.8 dt ".._" lw 1.5    lc rgb "#FF3300",\
'' using 22 title '{/Times O-DISTA}' with lp pt 8 pi 5 ps 0.8 dt "._" lw 1.5  lc rgb "#3333FF",\
