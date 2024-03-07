#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

FigTitle = "AdiabaticCollapse_XCFC"

ta  = np.loadtxt( 't_UniGrid_dr1.0km.dat' , skiprows = 1 ) / 2.99792458e8 * 1.0e3
dta = np.loadtxt( 'dt_UniGrid_dr1.0km.dat', skiprows = 1 ) / 2.99792458e8 * 1.0e3



tb  = np.loadtxt( 't_amr_dr1.0km.dat' , skiprows = 1 ) / 2.99792458e8 * 1.0e3
dtb = np.loadtxt( 'dt_amr_dr1.0km.dat', skiprows = 1 ) / 2.99792458e8 * 1.0e3

fig, ax = plt.subplots( 1, 1 )

ax.plot( ta, dta * 1.0e3,       \
         color  = 'red',        \
         marker = '',           \
         label  = 'Unigrid'     )

ax.plot( tb, dtb * 1.0e3,       \
         color  = 'blue',        \
         marker = '',           \
         label  = 'AMReX'     )

ax.legend(  loc = "upper right" )

ax.set_xlabel( r'$t\,\left[\mathrm{ms}\right]$' )
ax.set_ylabel( r'$dt\,\left[\mu\mathrm{s}\right]$' )

ax.set_yscale('log')

ax.set_title( '{:} - Time Step over Time'.format(FigTitle))

SaveFig = True
if SaveFig:

    FigName = 'fig.{:s}_dtvst.png'.format( FigTitle )
    plt.savefig(FigName, dpi = 300 )

else:

    plt.show()
