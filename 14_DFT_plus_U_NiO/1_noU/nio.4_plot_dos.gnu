set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Density of states (DOS) of NiO with standard DFT"
set xlabel "Energy (eV)"
set ylabel "DOS"
set arrow 1 from 14.432,-15 to 14.432,15 nohead ls 10 dt 2
#set xr [5:25]
#set yr [0:325]
plot    "nio.dos.dat" using 1:2 title 'DOS of spin-up' with lines,\
	"nio.dos.dat" using 1:(-$3) title 'DOS of spin-down' with lines
pause -1 "Hit any key to continue\n"    #so that the code doesn't exit automatically
