#This is a simple program that will plot data from the 
#file file1.dat

Plot "file1.dat"
#This will plot the graph with the key/legend as file1.dat. To remove this use
unset key

#To use a custom legend:
set key

plot "file1.dat" title "my plot" lt 7 lc 0 w lp

#We can plot multiple plots together. Just type the data informations separated by a comma. Let's plot a sine fuction along with the same data.
p "file1.dat" t "My first plot" lc 0 lt 7, 5*sin(x) lt 7 lc 1 t "sine function"

#x and y range
set xrange [-2: 12]
set yrange [-5:10]
replot

#Axis title
set xlabel "x-axis"
set ylabel "y-axis"
replot 
