set term pdf enhanced 
set key font ",16" 
set key spacing 1.1
set grid
set xlabel '{/Times time (s)}' font ",22"
unset label
unset border

set mytics 2
set mxtics 2
set grid xtics ytics mxtics mytics

set xrange [0:1]

set output 'a.pdf'
set key top samplen 2.5
set ylabel '{/Times-Italic a_1}' font ",22"
plot 'patch_a1.txt' using ($1)*0.001:2 notitle  with l dt "." lw 5 lc rgb "#606060",\
'csi.dat'  using ($1)*0.001:2 title "{/Times true}" with l dt "." lw 5 lc rgb "#606060",\
'' using ($1)*0.001:3 title "{/Times O-IST - Alg. 1}" with lp pt 6 pi 50 ps 0.8  lw 1.5  lc rgb "#00CC66",\
'' using ($1)*0.001:27 title "{/Times O-DR - Alg. 3}" with lp pt 2 pi 50 ps 0.8  lw 1.5    lc rgb "#FF3300",\
'' using ($1)*0.001:51 title "{/Times O-DISTA - Alg. 5}" with lp pt 9 pi 50 ps 0.8  lw 1.5  lc rgb "#3333FF",\

set output 'b.pdf'
set key bottom samplen 2.5
set ylabel '{/Times-Italic b_1}' font ",22"
plot  'patch_a2.txt' using ($1)*0.001:2 notitle with l dt "." lw 5 lc rgb "#606060",\
'csi.dat'  using ($1)*0.001:4 title "{/Times true}" with l dt "." lw 5 lc rgb "#606060",\
'' using ($1)*0.001:5 title "{/Times O-IST - Alg. 1}" with lp pt 6 pi 50 ps 0.8  lw 1.5  lc rgb "#00CC66",\
'' using ($1)*0.001:29 title "{/Times O-DR - Alg. 3}" with lp pt 2 pi 50 ps 0.8  lw 1.5    lc rgb "#FF3300",\
'' using ($1)*0.001:53 title "{/Times O-DISTA - Alg. 5}" with lp pt 9 pi 50 ps 0.8  lw 1.5  lc rgb "#3333FF",\


set key top samplen 2.5
set output 'zero.pdf'
set ylabel '{/Times null parameters (mean)}' font ",22"
plot 'csi.dat'  using ($1)*0.001:(sum [col=6:23] column(col))/18 title "{/Times O-IST - Alg. 1}"  with lp pt 6 pi 50 ps 0.7  lw 1.5  lc rgb "#00CC66",\
'' using ($1)*0.001:(sum [col=30:47] column(col))/18 title "{/Times O-DR - Alg. 3}" with lp  pt 2 pi 50 ps 0.7 lw 1.5    lc rgb "#FF3300",\
'' using ($1)*0.001:(sum [col=54:71] column(col))/18 title "{/Times O-DISTA - Alg. 5}" with lp  lw 1.5  pt 9 pi 50 ps 0.7 lc rgb "#3333FF",\
'' using ($1)*0.001:($4)*0 title "{/Times true}" with l dt "." lw 5 lc rgb "#606060",\


set output 'mse.pdf'
set yrange [0.0002:0.2]
set ytics 0.0002,10,0.2
set logscale y
set key top samplen 2.5
set ylabel "{/Times MSE}" font ",22"
plot 'csi.dat'  using ($1)*0.001:24 title  "{/Times O-IST - Alg. 1}" with lp pt 6 pi 50 ps 0.8  lw 1.5  lc rgb "#00CC66",\
'' using ($1)*0.001:48 title "{/Times O-DR - Alg. 3}" with lp    pt 2 pi 50 ps 0.8 lw 1.5    lc rgb "#FF3300",\
'' using ($1)*0.001:72 title "{/Times O-DISTA - Alg. 5}" with lp   pt 9 pi 50 ps 0.8 lw 1.5  lc rgb "#3333FF",\

set output 'dr.pdf'
unset logscale y
set yrange [0:0.1]
set ytics 0.01

set key top samplen 2.5
set ylabel '{/Times Reg^d_t / t}' font ",22"
plot 'csi.dat' using ($1)*0.001:($74/$1) title  "{/Times O-IST - Alg. 1}" with lp pt 6 pi 50 ps 0.8  lw 1.5  lc rgb "#00CC66",\
'' using ($1)*0.001:($75/$1) title "{/Times O-DR - Alg. 3}" with lp    pt 2 pi 50 ps 0.8 lw 1.5    lc rgb "#FF3300",\
'' using ($1)*0.001:($76/$1) title "{/Times O-DISTA - Alg. 5}" with lp   pt 9 pi 50 ps 0.8 lw 1.5  lc rgb "#3333FF"



