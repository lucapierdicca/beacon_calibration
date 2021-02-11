# Beacon Calibration
Calibrate beacon positions in 2D (3D with fixed z) starting from 3 known anchors and beacon 1-hop range measurements.  
The calibration is formulated as a NLLS problem - GN algorithm is used for the solution.

### requirements
scikit-sparse  
https://scikit-sparse.readthedocs.io/en/latest/overview.html

### requirements installation
<code>sudo apt-get install libsuitesparse-dev</code>  
<code>pip install --user scikit-sparse</code>

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
Results of the best performing method (NLLS) are reported and commented here  

![alt text](https://github.com/lucapierdicca/beacon_calibration/blob/main/Tri_NL.png)
![alt text](https://github.com/lucapierdicca/beacon_calibration/blob/main/Cal_NL.png)
<p align="center">
  <img src="https://github.com/lucapierdicca/beacon_calibration/blob/main/Chi_NLsmall.png">
</p>


### known limitations
sparse matrix in memory  
precomputation of state components permutation 
difficult trilateration in case of 2 surveyed beacons  
iterative trilateration accumulates error 
