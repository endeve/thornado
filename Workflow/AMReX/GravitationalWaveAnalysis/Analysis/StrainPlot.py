from math import pi

import numpy as np
import matplotlib.pyplot as plt
import subprocess
import math as m
from UtilitiesModule import ChoosePlotFile


# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

Root = 'NR2D_M2.8_Mdot0.3_Rs120'
WritingDirectory = HOME + 'Research/Data/'
F = open( '{:}{:}{:}.txt'.format( WritingDirectory, Root, '_Data' ), 'r' )

Length = 0
#for l in F:
#    Length = Length + 1

T = []
rC = []
rP = []
#Times = np.zeros( shape = (Length) ) 
#Rh_c  = np.zeros( shape = (Length) )
#Rh_p  = np.zeros( shape = (Length) )
lines = F.readlines()

for line in lines:
    print(line)
    Time, rh_c, rh_p = line.split(" ")
    
    Time = float(Time)
    rh_c = float(rh_c)
    rh_p = float(rh_p)

    T.append(Time)
    rC.append(rh_c)
    rP.append(rh_p)
    
    Length = Length + 1

F.close()

Times = np.asarray( T ) 
Rh_c  = np.asarray( rC )
Rh_p  = np.asarray( rP )

plt.plot( Times * 1000, Rh_c, color = 'black' )
plt.plot( Times * 1000, Rh_p, color = 'red' )
plt.xlabel("Postbounce Time (ms)")
plt.ylabel("D h (cm)")
plt.legend( ['$D h_x$','$D h_+$'], loc = "upper left" )
#plt.plot( X1_C[:], rh_r[nFiles-2,0,:] )
#plt.plot( X1_C[:], rh_r[nFiles-2,1,:] )
plt.show()
