import os
import csv
import numpy as np






def load(data, ground_truth):

	# X_true
	X_true = []
	with open(ground_truth, 'r') as file:
	    reader = csv.reader(file, skipinitialspace=True)
	    for row in reader:
	    	X_true.append(row[1:])
	X_true = np.array(X_true, dtype=np.float64).T


	# Z and Anchored beacons
	Z, A, U = {},{},{}
	with open(data, 'r') as file:
	    reader = csv.reader(file, skipinitialspace=True)
	    for index, row in enumerate(reader):
	    	if index < 3:
	    		ID = int(row[0])-1 
	    		A[ID] = np.array(row[2:], dtype=np.float64)
	    	else:
		    	i, j = int(row[0])-1, int(row[1])-1
		    	d = float(row[2])
		    	Z[(i, j)] = d
		    	
		    	if i not in U:
		    		U[i] = {'out':1, 'in':0}
		    	else:
		    		U[i]['out'] += 1
		    	
		    	if j not in U:
		    		U[j] = {'out':0, 'in':1}
		    	else:
		    		U[j]['in'] += 1


	return X_true, Z, A, U


# load the data (.csv to memory)------------------------------------------
# X_true: true state
# Z		: measurements
# A		: anchored beacons

os.system("octave-cli dat2csv.m")
files = os.listdir('.')
assert('beacon_data.csv' in files and 'beacon_gt.csv' in files)

X_true, Z, A, U = load('beacon_data.csv', 'beacon_gt.csv')
X_guess = np.zeros(X_true.shape)
state_var_dim = X_true[:,0].shape[0]
state_dim = X_true.shape[1]
omega_z = 3

nnz = len(Z)*state_var_dim**2 *2 + len(U)*state_var_dim**2

print('Beacons num: %d' % X_true.shape[1])
print('Measurements num: %d' % len(Z))
print('H shape: ', (state_dim*state_var_dim, state_dim*state_var_dim))
print('H elements: %d' % (state_var_dim**2*state_dim**2))
print('H nnz: %d\n' % nnz)