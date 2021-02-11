import csv
import time
import os
import numpy as np
import glbls
from trilateration_utils import *
from calibration_utils import *



#-----------------------TRILATERATION------------------------

# initialize the state initial guess X_guess 
# by iterative trilateration starting from anchors
# 1) build the beacons connectivity graph
# 2) visit the graph starting from the 1st anchor
# 3) while visiting survey the unknown beacons 
#    (available modes: GEO - LS - NLLS)

MODE = 'NLLS'
connectivity_graph = {}

# build the graph
for IDs,z in glbls.Z.items():
	i,j = IDs[0],IDs[1]
	
	if i not in connectivity_graph:
		connectivity_graph[i] = {'adj':[j],'vi':False,'su':False,'coord':[]}
		
		if i in glbls.A:
			connectivity_graph[i]['coord'] = glbls.A[i]
			connectivity_graph[i]['su'] = True
	else:
		connectivity_graph[i]['adj'].append(j)

# visit the graph and survey the beacons
print('Trilateration started')
tic = time.time()
visit(0, connectivity_graph, glbls.Z, mode=MODE)
tac = tac = time.time()
print('Trilateration ended')
print('Elapsed time %.2f s' % (tac-tic))

# update X_guess
for ID,data in connectivity_graph.items():
	glbls.X_guess[:,ID] = data['coord']

# compute Euc error per beacon
x,y = [],[]
for c in range(glbls.X_guess.shape[1]):
	x_pred = glbls.X_guess[:,c].reshape(3,1)
	x_true = glbls.X_true[:,c].reshape(3,1)
	euc = np.linalg.norm(x_true - x_pred)
	x.append(c)
	y.append(euc)


print('Max Euc error: %.4f' %  max(y))
print('Avg Euc error: %.4f' %  (sum(y)/len(y)))

# plot the results
tri_y = list(y)
plot_tri(x,y) 
print()

# --------------------CALIBRATION---------------------------

# refine the state initial guess exploiting NLLS
N_iter = 20
chi_log = []

# execute refinement iteration
print('Calibration started')
tic = time.time()
for _ in range(N_iter):	
	chi = calibration_NLLS()
	chi_log.append(chi)
	print('Chi square: %.5f' % chi)
tac = time.time()
print('Calibration ended')
print('Elapsed time: %.4f s' % (tac-tic))

# compute Euc error per beacon
x,y = [],[]
for c in range(glbls.X_guess.shape[1]):
	x_pred = glbls.X_guess[:,c].reshape(3,1)
	x_true = glbls.X_true[:,c].reshape(3,1)
	euc = np.linalg.norm(x_true - x_pred)
	x.append(c)
	y.append(euc)


print('Max Euc error: %.4f' %  max(y))
print('Avg Euc error: %.4f' %  (sum(y)/len(y)))

# plot the results
plot_cal(chi_log, x, y, tri_y)
print()

