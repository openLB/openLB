reset
set grid
set terminal pdf enhanced
set xrange [0:6]
set xtics  0,1,6
set xlabel "Time [s]"

set output "x_coordinates.pdf"
set ylabel "horizontal position [m]"
set yrange [0:0.02]
set ytics  0,0.002,0.02
set key left bottom
plot "gnuplot.dat" using 1:4 with lines title "Leading Particle", \
     "gnuplot.dat" using 1:5 with lines title "Trailing Particle"

set output "y_coordinates.pdf"
set ylabel "vertical position [m]"
set yrange [0:0.08]
set ytics  0,0.005,0.08
set key left bottom
plot "gnuplot.dat" using 1:2 with lines title "Leading Particle", \
     "gnuplot.dat" using 1:3 with lines title "Trailing Particle"
