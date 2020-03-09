import numpy as np
import matplotlib.pyplot as plt

y=0
x = 8
V = 0.5
directory = "./OutputV05/"
means = np.zeros([6], dtype=np.float64 )
nx = np.zeros([6], dtype=np.float64 )
while x <= 256 :
    dataList = [[float(s) for s in item.strip().split(",")] for item in open(directory + "thornado" + str(x) + "init.txt").readlines()]
    rows = len(dataList)
    cols = len(dataList[0])
    arrinit = np.zeros([rows,cols], dtype=np.float64)
    for i, row in enumerate(dataList):
        for j, number in enumerate(row):
            arrinit[i][j] = number
    dataList = [[float(s) for s in item.strip().split(",")] for item in open(directory + "thornado" + str(x) + "fin.txt").readlines()]
    rows = len(dataList)
    cols = len(dataList[0])
    arrfin = np.zeros([rows,cols], dtype=np.float64)
    for i, row in enumerate(dataList):
        for j, number in enumerate(row):
            arrfin[i][j] = number

    means[y] = np.mean(abs(arrfin-arrinit))
    nx[y] = x
    y = y + 1
    x = x * 2 

plt.loglog(nx,means, label = 'Data')
plt.loglog(nx,1/(nx*nx), label = 'Second Order')
plt.legend()
plt.xlabel('nX')
plt.ylabel('Error')
plt.title('V = ' + str(V))
plt.show()