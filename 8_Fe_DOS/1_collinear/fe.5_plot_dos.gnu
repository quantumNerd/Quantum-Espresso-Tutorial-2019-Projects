set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Density of states (DOS) of Fe crystal"
set xlabel "Energy (eV)"
set ylabel "DOS"
set arrow 1 from 22.090,-2.5 to 22.090,2.5 nohead ls 10 dt 2
#set key 0.01,100
set xr [10:35]
#set yr [0:325]
plot    "fe.dos.dat" using 1:2 title 'DOS of spin-up' with linespoints,\
	"fe.dos.dat" using 1:(-$3) title 'DOS of spin-down' with linespoints
pause -1 "Hit any key to continue\n"    #so that the code doesn't exit automatically
