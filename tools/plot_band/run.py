#!/usr/bin/env python
from Bands import *

datafile='si.bands.dat.gnu'
fermi = 6.1330
symmetryfile='si.3_bands.pp.out'
bool_shift_efermi= True
fig, ax = plt.subplots()

#bndplot(datafile,fermi,symmetryfile,ax)
bndplot(datafile,fermi,symmetryfile,ax,shift_fermi=0,\
color='black',linestyle='solid',name_k_points=['L','G','X','U','G'],legend='Si, PBE')


#fig.savefig("test.png")
plt.show()
