import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


Dinitial=np.loadtxt('Dinitial.dat');

I1initial=np.loadtxt('I1initial.dat');

I2initial=np.loadtxt('I2initial.dat');

I3initial=np.loadtxt('I3initial.dat');

Dfinal=np.loadtxt('Dfinal.dat');

I1final=np.loadtxt('I1final.dat');

I2final=np.loadtxt('I2final.dat');

I3final=np.loadtxt('I3final.dat');

nIter=np.loadtxt('nIterations.dat');

V=np.loadtxt('vel.txt');


l=len(Dinitial)

ErrorD=np.zeros([l-1], dtype=np.float64 )
ErrorI1=np.zeros([l-1], dtype=np.float64 )
ErrorI2=np.zeros([l-1], dtype=np.float64 )
ErrorI3=np.zeros([l-1], dtype=np.float64 )
N=np.zeros([l-1], dtype=np.float64 )

for i in range(0,l-1):
    ErrorD[i]=100*abs((Dinitial[i+1]-Dfinal[i+1])/Dinitial[i+1])
    ErrorI1[i]=100*abs((I1initial[i+1]-I1final[i+1])/I1initial[i+1])
    ErrorI2[i]=100*abs((I2initial[i+1]-I2final[i+1])/I2initial[i+1])
    ErrorI3[i]=100*abs((I3initial[i+1]-I3final[i+1])/I3initial[i+1])


for i in range(0,l-1):
    if(ErrorD[i]>10.0):
        print(ErrorD[i], Dinitial[i+1],I1initial[i+1],I2initial[i+1],I3initial[i+1])

print(V)
 
for i in range(0,l-1):
    N[i]=nIter[i+1]


Davg=0
for i in range(0,l-1):
    Davg=Davg+ErrorD[i]
Davg=Davg/(l-1)

I1avg=0
for i in range(0,l-1):
    I1avg=I1avg+ErrorI1[i]
I1avg=I1avg/(l-1)

I2avg=0
for i in range(0,l-1):
    I2avg=I2avg+ErrorI2[i]
I2avg=I2avg/(l-1)

I3avg=0
for i in range(0,l-1):
    I3avg=I3avg+ErrorI3[i]
I3avg=I3avg/(l-1)


Navg=0
for i in range(0,l-1):
    Navg=Navg+N[i]
Navg=Navg/(l-1)




Avg=np.zeros([5], dtype=np.float64 )

Avg[0]=Navg
Avg[1]=Davg
Avg[2]=I1avg
Avg[3]=I2avg
Avg[4]=I3avg

np.savetxt("092M0.txt",Avg, delimiter=',')

Std=np.zeros([5], dtype=np.float64 )

Dstd=np.std(ErrorD)
I1std=np.std(ErrorI1)
I2std=np.std(ErrorI2)
I3std=np.std(ErrorI3)
Nstd=np.std(N)


Std[0]=Nstd
Std[1]=Dstd
Std[2]=I1std
Std[3]=I2std
Std[4]=I3std

np.savetxt("092M0std.txt",Std, delimiter=',')



minD=np.log10(np.min(ErrorD[np.nonzero(ErrorD)]))
maxD=np.log10(np.max(ErrorD[np.nonzero(ErrorD)]))

minI1=np.log10(np.min(ErrorI1[np.nonzero(ErrorI1)]))
maxI1=np.log10(np.max(ErrorI1[np.nonzero(ErrorI1)]))

minI2=np.log10(np.min(ErrorI2[np.nonzero(ErrorI2)]))
maxI2=np.log10(np.max(ErrorI2[np.nonzero(ErrorI2)]))

minI3=np.log10(np.min(ErrorI3[np.nonzero(ErrorI3)]))
maxI3=np.log10(np.max(ErrorI3[np.nonzero(ErrorI3)]))


N=N.astype(int)
minN=np.min(N)
maxN=np.max(N)




fig1=plt.figure(1)
plt.hist(ErrorD, bins=10**np.linspace(minD, maxD, 100), weights=np.ones(len(ErrorD)) / len(ErrorD))
plt.gca().set_xscale("log")
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Percent Error D', fontsize='large')
plt.title("|V| =  " +str(V[0]))

fig2=plt.figure(2)
plt.hist(ErrorI1, bins=10**np.linspace(minI1, maxI1, 100), weights=np.ones(len(ErrorI1)) / len(ErrorI1))
plt.gca().set_xscale("log")
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Percent Error I1', fontsize='large')
plt.title("|V| =  " +str(V[0]))

fig3=plt.figure(3)
plt.hist(ErrorI2, bins=10**np.linspace(minI2, maxI2, 100), weights=np.ones(len(ErrorI2)) / len(ErrorI2))
plt.gca().set_xscale("log")
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Percent Error I2', fontsize='large')
plt.title("|V| =  " +str(V[0]))

fig4=plt.figure(4)
plt.hist(ErrorI3, bins=10**np.linspace(minI3, maxI3, 100), weights=np.ones(len(ErrorI3)) / len(ErrorI3))
plt.gca().set_xscale("log")
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Percent Error I3', fontsize='large')
plt.title("|V| =  " +str(V[0]))


fig5=plt.figure(5)
plt.hist(N,bins=range(minN, maxN+1, 1), weights=np.ones(len(N)) / len(N))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Number of iterations', fontsize='large')
plt.title("|V| =  " +str(V[0]))

