#Contour plot
#In gnu-lot type: load 'contour.p'
reset
set xrange [0:2*pi]
set yrange [0:2*pi]
set contour
set isosamples 30
set view map scale 0.5
unset key
unset surface
f(x,y)=sin(x)*sin(y)+cos(y)
splot f(x,y)
