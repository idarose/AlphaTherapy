import matplotlib.pyplot as plt
import numpy as np

analysisRegion = "cells"
activity = 10
cellLine = "C4-2"

filePath = f"../Mathematica/Output/{cellLine}/{analysisRegion}/Activity{activity}_"


dataN212PbTotal = np.loadtxt(f"{filePath}N212PbTotal.csv", delimiter=",")

dataN212Pb = np.loadtxt(f"{filePath}N212Pb.csv", delimiter=",")

dataN212Bi = np.loadtxt(f"{filePath}N212Bi.csv", delimiter=",")

dataN208Tl = np.loadtxt(f"{filePath}N208Tl.csv", delimiter=",")

dataN212Po = np.loadtxt(f"{filePath}N212Po.csv", delimiter=",")

dataN208Pb = np.loadtxt(f"{filePath}N208Pb.csv", delimiter=",")


plt.title(f"Number of radionuclides in {analysisRegion}")
plt.plot(dataN212Pb[:,0],dataN212Pb[:,1],label="N212Pb")
plt.plot(dataN212Po[:,0],dataN212Po[:,1],label="N212Po")
plt.plot(dataN212Bi[:,0],dataN212Bi[:,1],label="N212Bi")
plt.plot(dataN208Tl[:,0],dataN208Tl[:,1],label="N208Tl")
plt.plot(dataN208Pb[:,0],dataN208Pb[:,1],label="N208Pb")
plt.semilogy()
plt.legend()
plt.savefig(f"Figures/A{activity}_N_{cellLine}_{analysisRegion}.png")
plt.clf()


dataDecays212Pb212Bi = np.loadtxt(f"{filePath}TotalDecays212Pb212Bi.csv", delimiter=",")

dataDecays212Bi212Po = np.loadtxt(f"{filePath}TotalDecays212Bi212Po.csv", delimiter=",")

dataDecays212Bi208Tl = np.loadtxt(f"{filePath}TotalDecays212Pb212Bi.csv", delimiter=",")

dataDecays208Tl208Pb = np.loadtxt(f"{filePath}TotalDecays208Tl208Pb.csv", delimiter=",")

dataDecays212Po208Pb = np.loadtxt(f"{filePath}TotalDecays212Po208Pb.csv", delimiter=",")


plt.title(f"Total Number of Decays in {analysisRegion}")
plt.plot(dataDecays212Pb212Bi[:,0],dataDecays212Pb212Bi[:,1],label="212Pb->212Bi")
plt.plot(dataDecays212Bi212Po[:,0],dataDecays212Bi212Po[:,1],label="212Bi->212Po")
plt.plot(dataDecays212Bi208Tl[:,0],dataDecays212Bi208Tl[:,1],label="212Bi->208Tl")
plt.plot(dataDecays208Tl208Pb[:,0],dataDecays208Tl208Pb[:,1],label="208Tl->208Pb")
plt.plot(dataDecays212Po208Pb[:,0],dataDecays212Po208Pb[:,1],label="212Po->208Pb")
plt.semilogy()
plt.legend()
plt.savefig(f"Figures/A{activity}_D_{cellLine}_{analysisRegion}.png")
plt.clf()



analysisRegion = "solution"
activity = 10
cellLine = "C4-2"

filePath = f"../Mathematica/Output/{cellLine}/{analysisRegion}/Activity{activity}_"


dataN212PbTotal = np.loadtxt(f"{filePath}N212PbTotal.csv", delimiter=",")

dataN212Pb = np.loadtxt(f"{filePath}N212Pb.csv", delimiter=",")

dataN212Bi = np.loadtxt(f"{filePath}N212Bi.csv", delimiter=",")

dataN208Tl = np.loadtxt(f"{filePath}N208Tl.csv", delimiter=",")

dataN212Po = np.loadtxt(f"{filePath}N212Po.csv", delimiter=",")

dataN208Pb = np.loadtxt(f"{filePath}N208Pb.csv", delimiter=",")

ratio = dataN212Pb[:,1]/dataN212PbTotal[:,1]

plt.title(f"PDF, radionuclide in solution")
plt.plot(dataN212Pb[:,0], ratio)
plt.savefig("Figures/RatioPlot.png")
plt.clf()




plt.title(f"Number of radionuclides in {analysisRegion}")
plt.plot(dataN212Pb[:,0],dataN212Pb[:,1],label="N212Pb")
plt.plot(dataN212Po[:,0],dataN212Po[:,1],label="N212Po")
plt.plot(dataN212Bi[:,0],dataN212Bi[:,1],label="N212Bi")
plt.plot(dataN208Tl[:,0],dataN208Tl[:,1],label="N208Tl")
plt.plot(dataN208Pb[:,0],dataN208Pb[:,1],label="N208Pb")
plt.semilogy()
plt.legend()
plt.savefig(f"Figures/A{activity}_N_{cellLine}_{analysisRegion}.png")
plt.clf()



dataDecays212Pb212Bi = np.loadtxt(f"{filePath}TotalDecays212Pb212Bi.csv", delimiter=",")

dataDecays212Bi212Po = np.loadtxt(f"{filePath}TotalDecays212Bi212Po.csv", delimiter=",")

dataDecays212Bi208Tl = np.loadtxt(f"{filePath}TotalDecays212Pb212Bi.csv", delimiter=",")

dataDecays208Tl208Pb = np.loadtxt(f"{filePath}TotalDecays208Tl208Pb.csv", delimiter=",")

dataDecays212Po208Pb = np.loadtxt(f"{filePath}TotalDecays212Po208Pb.csv", delimiter=",")


plt.title(f"Total Number of Decays in {analysisRegion}")
plt.plot(dataDecays212Pb212Bi[:,0],dataDecays212Pb212Bi[:,1],label="212Pb->212Bi")
plt.plot(dataDecays212Bi212Po[:,0],dataDecays212Bi212Po[:,1],label="212Bi->212Po")
plt.plot(dataDecays212Bi208Tl[:,0],dataDecays212Bi208Tl[:,1],label="212Bi->208Tl")
plt.plot(dataDecays208Tl208Pb[:,0],dataDecays208Tl208Pb[:,1],label="208Tl->208Pb")
plt.plot(dataDecays212Po208Pb[:,0],dataDecays212Po208Pb[:,1],label="212Po->208Pb")
plt.semilogy()
plt.legend()
plt.savefig(f"Figures/A{activity}_D_{cellLine}_{analysisRegion}.png")
plt.clf()


