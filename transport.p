reset
set terminal png
set xrange [0: 2*pi]
set yrange [-2:2]
set output 'transport.png'
plot 'fort.100', 'fort.200'
replot
reset

