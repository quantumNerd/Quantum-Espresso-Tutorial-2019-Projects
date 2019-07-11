set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "absorption spectrum of benzene molecule"
set xlabel "Energy (eV)"
set ylabel "Absoprtion"
plot 'C6H6.plot_S.dat' u 1:2 title 'light absorption: C6H6' with linespoints
pause -1 "Hit any key to continue\n"    #so that the code doesn't exit automatically
