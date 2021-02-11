from scipy.sparse import csc_matrix
from sksparse.cholmod import cholesky
import matplotlib.pyplot as plt
import numpy as np
import glbls




def x_is2H_is(x_is):
	H_is = []
	for i in x_is:
		H_is += [i*3, (i*3)+1, (i*3)+2]
	return H_is

# compute the error ij
def error(state_i, state_j, measurement):
	pred_measurement =  np.linalg.norm_ij(state_i-state_j)
	e = pred_measurement - measurement
	return e

# compute the jacobian block ij
def jacobianBlocks(state_i, state_j):
	J = 1.0/np.linalg.norm_ij(state_i-state_j)* \
			np.array([[state_i[0]-state_j[0], state_i[1]-state_j[1], state_i[2]-state_j[2]]])

	return J, -J

# refine the initial state guess 
# exploiting NLLS
def calibration_NLLS():


	anchor_num = len(glbls.A)
	anchor_ID = list(glbls.A.keys())

	H = np.zeros((glbls.state_var_dim*(glbls.state_dim-anchor_num),
			  glbls.state_var_dim*(glbls.state_dim-anchor_num)))
	b = np.zeros((glbls.state_var_dim*(glbls.state_dim-anchor_num),1))
	c = 0	

	# for each (i,j) measurement
	for ids,z_ij in glbls.Z.items():
		if ids[0] not in anchor_ID and ids[1] not in anchor_ID:

			i,j = ids[0],ids[1]

			x_i = glbls.X_guess[:,i]
			x_j = glbls.X_guess[:,j]
			
			# -----------original-------------
			#
			# e_ij = error(x_i, x_j, z)
			# J_i, J_j = jacobianBlocks(x_i, x_j)

			# H_ii = J_i.T@omega_z@J_i
			# H_ij = J_i.T@omega_z@J_j
			# H_jj = J_j.T@omega_z@J_j

			# b_i = J_i.T@omega_z@e_ij
			# b_j = J_j.T@omega_z@e_ij

			# c += e_ij@omega_z@e_ij
			#-----------------------------------


			# compute (i,j) error and error jacobian block
			diff = (x_i-x_j).reshape(3,1)
			norm2_ij = np.dot(diff.T, diff)[0,0]
			norm_ij = np.sqrt(norm2_ij)
			f = 1/norm2_ij

			H_ii = f*glbls.omega_z*(diff@diff.T)
			H_ij = -H_ii
			H_jj = H_ii

			e_ij = norm_ij-z_ij
			chi_ij = e_ij*e_ij

			# if chi_ij>0.007:
			# 	print('OUTLIER')
			# 	e_ij = e_ij * np.sqrt(1.0/chi_ij)
			# 	chi_ij = 0.007
			
			b_i = (e_ij/norm_ij)*glbls.omega_z*diff
			b_j = -b_i

			# hard coded index shift (anchors)
			if i<50: i-=2
			else: i-=3

			if j<50: j-=2
			else: j-=3

			# populate H matrix, b vector
			# update c the chi square error
			H[glbls.state_var_dim*i:glbls.state_var_dim*i+glbls.state_var_dim,
				glbls.state_var_dim*i:glbls.state_var_dim*i+glbls.state_var_dim] += H_ii

			H[glbls.state_var_dim*i:glbls.state_var_dim*i+glbls.state_var_dim,
				glbls.state_var_dim*j:glbls.state_var_dim*j+glbls.state_var_dim] += H_ij

			H[glbls.state_var_dim*j:glbls.state_var_dim*j+glbls.state_var_dim,
				glbls.state_var_dim*j:glbls.state_var_dim*j+glbls.state_var_dim] += H_jj

			H[glbls.state_var_dim*j:glbls.state_var_dim*j+glbls.state_var_dim,
				glbls.state_var_dim*i:glbls.state_var_dim*i+glbls.state_var_dim] += H_ij.T

			b[glbls.state_var_dim*i:glbls.state_var_dim*i+glbls.state_var_dim] += b_i

			b[glbls.state_var_dim*j:glbls.state_var_dim*j+glbls.state_var_dim] += b_j

			c += chi_ij

	

	# solve the square sparse and damped linear system with cholesky
	f = cholesky(csc_matrix(H + np.eye(glbls.state_var_dim*(glbls.state_dim-anchor_num))*0.1))
	dx = f(-b)

	# copy paste dx in dx_ 
	dx_ = np.zeros((glbls.state_var_dim*glbls.state_dim,1))

	for i in range(glbls.state_dim-anchor_num):
		
		# hard coded index unshift (anchors)
		if i < 48: j=i+2
		else: j=i+3

		dx_[j*3] = dx[i*3]
		dx_[(j*3)+1] = dx[(i*3)+1]
		dx_[(j*3)+2] = dx[(i*3)+2]

	# update the state
	glbls.X_guess = glbls.X_guess + np.reshape(dx_, glbls.X_guess.shape, 'F')


	return c

# produce plots
def plot_cal(chi_log, x, y, tri_y):
	f,a = plt.subplots(1,3, figsize=(14,4))
	a[0].scatter(glbls.X_true[0,:], glbls.X_true[1,:], s=5)
	a[0].scatter(glbls.X_guess[0,:], glbls.X_guess[1,:], s=5)
	f.suptitle('Calibration')
	a[0].set_title('Ground truth vs. estimated position')
	a[0].axis('off')

	a[1].bar(x,y)
	a[1].set_title('Euclidean error per beacon')
	a[1].set_ylabel('euclidean error')
	a[1].set_ylim([0.0,max(tri_y)])
	a[1].set_xlabel('beacon ID')


	a[2].scatter(glbls.X_true[0,:], glbls.X_true[1,:], c=np.array([[e/max(y),e/max(y),e/max(y)] for e in y]), s=5)
	a[2].set_title('Euclidean error gray map')
	a[2].axis('off')
	plt.show()


	plt.plot(chi_log)
	plt.title('Chi square error vs. iteration')
	plt.ylabel('chi square error')
	plt.xlabel('iteration')
	plt.show()
