import matplotlib.pyplot as plt
import sys
import numpy as np

ProblemName = sys.argv[1]
FileDir     = sys.argv[2]
PlotVar     = sys.argv[3]

TallyFiles = {
    "CM_D" : str(ProblemName) + ".Tally_BaryonicMass.dat",
    "CM_S1": str(ProblemName) + ".Tally_MHDMomentumX1.dat",
    "CM_S2": str(ProblemName) + ".Tally_MHDMomentumX2.dat",
    "CM_S3": str(ProblemName) + ".Tally_MHDMomentumX3.dat",
    "CM_E" : str(ProblemName) + ".Tally_MHDEnergy.dat",
    "Tem11": str(ProblemName) + ".Tally_Tem11.dat",
    "Tem33": str(ProblemName) + ".Tally_Tem33.dat"
}

TallyData = np.loadtxt(FileDir + "/" + TallyFiles[PlotVar], skiprows = 1)

RelDiff = TallyData[:,4] / TallyData[:,3] 

plt.plot(TallyData[:,0], RelDiff[:], linewidth = 2.5)
plt.ylim( -1e0, 1e0 )
plt.yscale("symlog", linthreshy = 1.0e-8)
plt.xlabel("Time [ms]")
plt.ylabel("Relative Difference in " + str(PlotVar))
plt.savefig(FileDir + "/" + TallyFiles[PlotVar] + ".png")
