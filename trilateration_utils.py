import matplotlib.pyplot as plt
import itertools
import numpy as np
import glbls




# given 2 circles data (C, r)
# return the intersection points 
def intersect(adj_data):

	w_c1 = np.array([[adj_data[0]['coord'][0]],
					 [adj_data[0]['coord'][1]]])
	
	w_c2 = np.array([[adj_data[1]['coord'][0]],
					 [adj_data[1]['coord'][1]]])
	
	r1 = adj_data[0]['range'] 
	r2 = adj_data[1]['range'] 

	diff = w_c2-w_c1
	theta = np.arctan2(diff[1,0],diff[0,0]);
	w_R_f1 = np.array([[np.cos(theta),-np.sin(theta)],
	          		   [np.sin(theta),np.cos(theta)]])
	     
	f1_c1 = np.zeros((2,1))
	f1_c2 = w_R_f1.T@(w_c2-w_c1)

	f1_pa = np.zeros((2,1))
	f1_pb = np.zeros((2,1))
	     
	f1_pa[0,0] = 0.5*((np.power(f1_c2[0,0],2)+(np.power(r1,2)-np.power(r2,2)))/f1_c2[0,0])

	d = np.power(r1,2)-np.power(f1_pa[0,0],2)
	if d < 0:
		return -1 # SOS gestire nel chiamante?
		
	f1_pa[1,0] = np.sqrt(d);
	f1_pb[0,0] = f1_pa[0,0];
	f1_pb[1,0] = -f1_pa[1,0];

	w_pa = w_R_f1@f1_pa + w_c1;
	w_pb = w_R_f1@f1_pb + w_c1;

	return [w_pa, w_pb]

# given <3 already surveyed beacons
# returns the coordinates of a third 
# adjacent beacon based on geometry chriteria
def lt3(adj_data):

	candidates = intersect(adj_data)

	if candidates != -1:

		w_pa = candidates[0]
		w_pb = candidates[1]

		d_a = np.linalg.norm(w_pa)
		d_b = np.linalg.norm(w_pb)

		if d_a <= d_b:
			return w_pb
		else:
			return w_pa
	else:
		print('WARNING: NO CANDIDATES lt3') # SOS
		avg = 0
		for i in adj_data:
			avg += i['coord'][:2].reshape(2,1)
		avg /= 2.0
		
		return avg

# given >=3 already surveyed beacons
# return the coordinates of a third 
# adjacent beacon based on geometry 
# and K-NN density chriteria
def gte3(adj_data):

	candidates_set = []
	
	combi_2 = list(itertools.combinations(adj_data, 2))

	for anchors in combi_2:
		candidates = intersect(anchors)
		if candidates != -1:
			candidates_set += candidates

	if len(candidates_set) == 0:
		print('WARNING: NO CANDIDATES gte3') # SOS
		avg = 0
		for i in adj_data:
			avg += i['coord'][:2].reshape(2,1)
		avg /= len(adj_data)
		
		return avg
	
	#pprint(candidates_set)

	K = (len(candidates_set)//2)-1

	D = np.zeros((len(candidates_set),len(candidates_set)))
	i=0
	while i < len(candidates_set):
		for j in range(i+1,len(candidates_set)):
			D[i,j] = np.linalg.norm(candidates_set[i] - candidates_set[j])
		i+=1

	D = D.T + D

	#print(D)


	D_dict = {} # row:[[col,d],[col,d],...,[col,d]]
	for r in range(D.shape[0]):
		D_dict[r] = []
		for c in range(D.shape[1]):
			D_dict[r].append([c,D[r,c]])
		D_dict[r].sort(key = lambda x: x[1])


	#pprint(D_dict)
	min_Kth_d, max_density_candidateID = 999, -1
	for ID, d_list in D_dict.items():
		Kth_d = d_list[K][1]
		if Kth_d < min_Kth_d:
			min_Kth_d = Kth_d
			max_density_candidate_ID = ID

	
	K_best_candidates_ID = D_dict[max_density_candidate_ID][1:K+1]
	#pprint(K_best_candidates_ID)
	centroid = np.zeros((2,1))
	for Kth_data in K_best_candidates_ID:
		Kth_ID = Kth_data[0]
		centroid += candidates_set[Kth_ID]

	centroid /= K

	return centroid

# switch between the possible surveying cases and modes
def survey(n_ID, graph, measurements, mode):

	surveyed_adj = []
	
	adj = graph[n_ID]['adj']
	for a in adj:
		if graph[a]['su']:
			surveyed_adj.append({'ID':a, 
								 'coord':graph[a]['coord'], 
								 'range':(measurements[(a,n_ID)]+measurements[(n_ID,a)])/2.0}) # SOS ho usato la media

	n_coord = np.zeros((2,1))

	#print(n_ID)
	#pprint(surveyed_adj)

	if mode == 'GEO':
		if len(surveyed_adj) >= 3:
			n_coord = gte3(surveyed_adj)
		else: 
			n_coord = lt3(surveyed_adj)

	elif mode == 'LS':
		if len(surveyed_adj) >= 3:
			n_coord = trilateration_LS(surveyed_adj)
		else:
			n_coord = lt3(surveyed_adj)

	elif mode == 'NLLS':
		if len(surveyed_adj) >= 3:
			n_coord = gte3(surveyed_adj)
		else: 
			n_coord = lt3(surveyed_adj)

		for _ in range(5):
			n_coord, chi = trilateration_NLLS(n_coord, surveyed_adj)


	#print(n_coord)

	# update graph
	graph[n_ID]['coord'] = np.array([n_coord[0,0], n_coord[1,0], 1.0])
	graph[n_ID]['euc'] = np.linalg.norm(np.array([n_coord[0,0], n_coord[1,0], 1.0]) - glbls.X_true[:,n_ID])
	graph[n_ID]['su'] = True

# visit the graph and survey the beacons (while visiting)
def visit(root_ID, graph, measurements, mode='GEO'):

	graph[root_ID]['vi'] = True
	frontier = []
	frontier.append(root_ID)

	while len(frontier) != 0:
		n = frontier.pop(0)
		if not graph[n]['su']:
			survey(n, graph, measurements, mode)
			graph[n]['su'] = True

		adj = graph[n]['adj']
		for a in adj:
			if not graph[a]['vi']:
				graph[a]['vi'] = True
				frontier.append(a)

# given some already surveyed beacons
# return the coordinate of a third 
# adjacent beacon exploiting NLLS
# (init state guess = geometric)
def trilateration_NLLS(state, measurements):
	state_var_dim = state.shape[0]
	state_dim = state.shape[1]
	
	H = np.zeros((state_var_dim*state_dim,
				  state_var_dim*state_dim))
	b = np.zeros((state_var_dim*state_dim,1))
	c = 0

	for su_adj in measurements:

		e = np.linalg.norm(su_adj['coord'][:2].reshape((2,1)) - state) - su_adj['range']

		J = (1/np.linalg.norm(su_adj['coord'][:2].reshape((2,1)) - state)) \
		*np.array([[state[0,0] - su_adj['coord'][0], state[1,0] - su_adj['coord'][1]]])

		H += J.T@J
		b += J.T*e
		c += e*e

	#det = np.linalg.det(H)
	#print(det)

	dx = np.linalg.solve(H + np.eye(state_var_dim*state_dim)*0.1, -b)
	new_state = state + dx

	return new_state, c

# given some already surveyed beacons
# return the coordinate of a third 
# adjacent beacon exploiting LS
def trilateration_LS(measurements):
	# A dimension
	m = len(measurements)-1
	n = 2

	# ATA, ATb initialization
	ATA = np.zeros((n,n))
	ATb = np.zeros((n,1))

	p_last = measurements[-1]
	for su_adj in measurements[:-1]:
		
		Ai = 2*np.array([[p_last['coord'][0] - su_adj['coord'][0],p_last['coord'][1] - su_adj['coord'][1]]])
		
		p_last_vector = p_last['coord'][:2].reshape((2,1))
		su_adj_vector = su_adj['coord'][:2].reshape((2,1))
		bi = (p_last_vector.T@p_last_vector)[0,0] - (su_adj_vector.T@su_adj_vector)[0,0] + (np.power(su_adj['range'],2) - np.power(p_last['range'],2))
		
		ATA += Ai.T@Ai
		ATb += Ai.T*bi

	x = np.linalg.solve(ATA + np.eye(n)*0.1, ATb)

	# A = []
	# b = []

	# p_last = measurements[-1]
	# for su_adj in measurements[:-1]:
	# 	A.append([(p_last['coord'][0] - su_adj['coord'][0])*2,(p_last['coord'][1] - su_adj['coord'][1])*2])
		
	# 	p_last_vector = p_last['coord'][:2].reshape((2,1))
	# 	su_adj_vector = su_adj['coord'][:2].reshape((2,1))
	# 	b.append((p_last_vector.T@p_last_vector)[0,0] - (su_adj_vector.T@su_adj_vector)[0,0] + (np.power(su_adj['range'],2) - np.power(p_last['range'],2)))

	# A = np.array(A)
	# b = np.asarray(b).reshape(len(b),1)

	
	#x = np.linalg.lstsq(A, b, rcond=None)[0]

	return x

# produce plots
def plot_tri(x,y):
	f,a = plt.subplots(1,3, figsize=(14,4))
	a[0].scatter(glbls.X_true[0,:], glbls.X_true[1,:], s=5)
	a[0].scatter(glbls.X_guess[0,:], glbls.X_guess[1,:], s=5)
	f.suptitle('Trilateration')
	a[0].set_title('Ground truth vs. estimated position')
	a[0].axis('off')


	a[1].bar(x,y)
	a[1].set_title('Euclidean error per beacon')
	a[1].set_ylabel('euclidean error')
	a[1].set_xlabel('beacon ID')


	a[2].scatter(glbls.X_true[0,:], glbls.X_true[1,:], c=np.array([[e/max(y),e/max(y),e/max(y)] for e in y]), s=5)
	a[2].set_title('Euclidean error gray map')
	a[2].axis('off')

	plt.show()