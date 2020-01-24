set term pdf
#set style data lines
set style arrow 1 head back empty size screen 0.013,15,130  linetype 1  linecolor rgb "black"
set style arrow 2 head  filled size screen 0.01,15,45  linetype 1  dt "_" linecolor rgb "#00CC66"
set style arrow 3 head  filled size screen 0.015,15,45  linetype 1  dt "_" linecolor rgb "#FF3300"
set style arrow 4 head filled size screen 0.015,15,45  linetype 1  dt "_" linecolor rgb "#3333FF"
set style line 1  linetype 1  linecolor rgb "#E9D92A"  dt "." lw 2.5


set xrange [0:24]
set yrange [0:24]
set xtics -0.5,25,24.5
set ytics -0.5,25,24.5
set mxtics 25
set mytics 25

set size ratio -1
set grid mxtics mytics

set key font "Times, 14"
set key at 11,22


set output "locist.pdf"
plot 'loc.dat' using 1:2:3:4 title '{/Times moving target}'  w vectors arrowstyle 1,\
'' using 5:6:7:8 title '{/Times O-IST }' w vectors arrowstyle 2,\
'sensors.dat' u ($1)/100:($2)/100 w p linestyle 1 pt 7 ps 0.2 notitle 

set output "locdr.pdf"
plot 'loc.dat' using 1:2:3:4 title '{/Times moving target}'  w vectors arrowstyle 1,\
'' using 9:10:11:12 title '{/Times O-DR}' w vectors arrowstyle 3,\
'sensors.dat' u ($1)/100:($2)/100 w p linestyle 1 pt 7 ps 0.2 notitle 

set output "locdista.pdf"
plot 'loc.dat' using 1:2:3:4 title '{/Times moving target}'  w vectors arrowstyle 1,\
'' using 13:14:15:16 title '{/Times O-DISTA}' w vectors arrowstyle 4,\
'sensors.dat' u ($1)/100:($2)/100 w p linestyle 1 pt 7 ps 0.4 notitle,\
'connections.dat' u ($1)/100:($2)/100 w l linestyle 1 notitle 


