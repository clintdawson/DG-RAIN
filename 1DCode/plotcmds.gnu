# set terminal  
set terminal jpeg

# Line width of the axes
set border linewidth 1.5

# Line styles
set style line 1 linecolor rgb "red" linetype 1 linewidth 2
set style line 2 linecolor rgb "blue" linetype 1 linewidth 2
set style line 3 linecolor rgb "magenta" linetype 1 linewidth 2

# specify the format of the input data and the format of the label on the x-axis

# format the axes
set key right top box 1
set xlabel "x"
set ylabel "Water Surface Elevation"
#set yrange[0:0.45]

set output "MacDonaldSmoothTransitionAndShockZeta1D200El.jpg"

plot "039_Zeta.dat" using 1:3 title "Bottom Elevation" with lines linestyle 1, \
"SmallAnalyticSol" using 1:6 title "Analytical Solution" with lines linestyle 2, \
"039_Zeta.dat" using 1:5 title "DG solution" with lines linestyle 3

#plot "analyticSol" using 1:2 title 'Bottom Elevation' with lines linestyle 1,\
#"analyticSol" using 1:2 title 'Analytical solution' with lines linestyle 2,\

#set key right top box 1
set key left top box 1
set xlabel "x"
set ylabel "Discharge"
#set yrange[0.15:0.2]
set output "MacDonaldSmoothTransitionAndShockQ1D200El.jpg"

plot "SmallAnalyticSol" using 1:5 title 'Analytical solution' with lines linestyle 2,\
"039_Q.dat" using 1:4 title "DG solution" with lines linestyle 3


