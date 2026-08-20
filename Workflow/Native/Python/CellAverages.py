import numpy as np

def CellAverage( dim, nNodes, quant ):

	W1D = np.zeros(nNodes)
	nX  = np.array( np.shape( quant ) ) // nNodes

	if(dim == 1):
		nX[0] = 1
		nX[1] = 1
	elif(dim == 2):
		nX[0] = 1

	Avg = np.zeros(shape = (nX[0],nX[1],nX[2]))

	if(nNodes == 1):
		W1D = [1.0]
	elif(nNodes == 2):
		W1D = [0.5, 0.5]
	elif(nNodes == 3):
		W1D = [5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0]

	if( dim == 1):
		WnD = np.zeros( shape = ( 1, 1, nNodes ) )
		WnD[0,0,:] = W1D
	elif( dim == 2 ):
		WnD = np.zeros( shape = ( 1, nNodes, nNodes ) )
		WnD[0,:,:] = np.multiply.outer( W1D, W1D )
	elif( dim == 3 ):
		WnD = np.zeros( shape = ( nNodes, nNodes, nNodes ) )
		WnD[:,:,:] = np.multiply.outer ( np.multiply.outer ( W1D, W1D ), W1D )


	for l in range(nX[2]):
		for m in range(nX[1]):
			for n in range(nX[0]):
				Avg[n,m,l] = 0
				for k in range(np.size(WnD[:,0,0])):
					for j in range(np.size(WnD[0,:,0])):
						for i in range(np.size(WnD[0,0,:])):
							Avg[n,m,l] = Avg[n,m,l] + WnD[k, j, i] * quant[k + nNodes * n,j + nNodes * m,i + nNodes * l]

	return Avg
