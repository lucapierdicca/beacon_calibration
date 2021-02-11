# Beacon Calibration
Calibrate beacon positions in 2D (3D with fixed z) starting from 3 known anchors and beacon 1-hop range measurements.  
The calibration is formulated as a NLLS problem - GN algorithm is used for the solution.

### requirements
scikit-sparse  
https://scikit-sparse.readthedocs.io/en/latest/overview.html

### requirements installation
<code>sudo apt-get install libsuitesparse-dev</code>  
<code>pip install --user scikit-sparse</code>

### run
'''MODE = {'GEO','LS','NLLS'}'''
'''python beacon_calibration_opt.py'''

### trilateration
The GN algorithm initial state guess is obtained by iterative trilateration of the beacons starting from 3 known anchors.  
In more details, iterative trilateration is realized leveraging on the construction of a connectivity graph G = (V,E) where each node b in V is a beacon_ID and each edge (bi_ID, bj_ID) in E represents the existence of a measurement between bi_ID and bj_ID. The graph is then visited and surveyed: for each selected beacon b it is defined a set S_b that contains all its adjacent and already surveyed beacons relevant data (position and range measurements). These data will finally be processed in differrent ways in order to localize b.  

The surveying of a node can be performed using 3 different modalities:
   * GEO: geometry and (k-NN)density based heuristics are used in order to select among the possible candidate solutions
   * LS: least squares - the surveying is casted into a linear least squares problem (A.T * A + a * I)x = A.T * b 
   * NLLS: non-linear least squares - the surveying is casted into a non-linear least squares problem (J.T * J + a * I)x = -J.T * b. The initial node position guess is computed by GEO. 

### calibration
Once the initial state guess is obtained it is refined using the whole set of range measurements.

### results
The results of the best performing trilateration method (NLLS) are reported and commented here.    
In the first row we can see the results after the trilateration is performed, in the second row are the results after the calibration and in the third row the chi square error evolution is reported.   
Looking at the central column graphs it is possible to see how a refinement based on the whole set of measurements helps to uniformly reduce the euclidean error. This is particularly evident for the beacons on the borders of the grid - those beacons are surveyed using only 2 adjacent beacons. In fact, from the first column graphs, we can clearly see that the mislocalization of the beacons is mainly caused by the mislocalization of the bottom grid row and first grid column beacons which, given the iterative nature of the algorithm, makes the error to accumulate over and over.   
It is not possible to compare the two gray maps because they are normalized wrt to their own maximum error value, but from them we can see the zones where the error is concentrated (the darker the better).   
Quantitatively the maximum euclidean error is halved (from 0.25 to 0.11) and the average euclidean error is also reduced (from 0.09 to 0.06)  
![alt text](https://github.com/lucapierdicca/beacon_calibration/blob/main/Tri_NL.png)
![alt text](https://github.com/lucapierdicca/beacon_calibration/blob/main/Cal_NL.png)
<p align="center">
  <img src="https://github.com/lucapierdicca/beacon_calibration/blob/main/Chi_NLsmall.png">
</p>


### known limitations
There are several known limitations that should be adressed in the future in order to make the whole procedure more robust and efficient:  

- memory and sparse matrices: the final H matrix dimension is (7500 * 7500) and it is mainly sparse. The total amount of memory consumption is ~ 500MB (64 bits per element) therefore this suggests smart sparse memory storage. In the current version this is not realized, the transformation into sparse format is just realized a posteriori (augmenting the computational complexity even further) in order to be able to exploit the sparse cholesky routine from scikit-sparse.  
- precomputation of state components order: to speed up the solution of the sparse linear system a precomputation of the permutation matrix should be done and the consequent reordering of the state components should be performed once at the beginning.   
- trilateration in case of 2 anchor beacons: the surveying method currently implemented for this kind of cases is not really robust and unfortunately leverages on some a priori knowledge of the geometry of the problem.
- iterative trilateration and error accumulation: iterative trilateration can perform poorly because of the error accumulation, this fact and the geometry of the problem (the beacons are placed in a grid) might heavily concurr to the reasons why the LS method for trilateration performs particularly bad. A better scheme might help reducing the error in general.
