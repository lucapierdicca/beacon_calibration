dataset = load('beacon_calibration_data.dat');
csv = 'beacon_data.csv';
fid = fopen(csv,'w');

anchors_num = size(dataset.data.anchors)(2);
measurements_num = size(dataset.data.measurements)(2);


for(i=1:anchors_num)
	
	fprintf(fid, "%d, %d, %.17g, %.17g, %.17g\n", 
		dataset.data.anchors(i).id,
		0,
		dataset.data.anchors(i).coordinates(1),
		dataset.data.anchors(i).coordinates(2),
		dataset.data.anchors(i).coordinates(3));

endfor;

for(i=1:measurements_num)
	
	fprintf(fid, "%d, %d, %.17g, %.17g, %.17g\n", 
		dataset.data.measurements(i).variables_ids(1),
		dataset.data.measurements(i).variables_ids(2),
		dataset.data.measurements(i).range,
		0.0,
		0.0);

endfor;

fclose(fid);



dataset = load('beacon_calibration_gt.dat');
csv = 'beacon_gt.csv';
fid = fopen(csv,'w');

beacons_num = size(dataset.gt.beacons)(2);

for(i=1:beacons_num)
	
	fprintf(fid, "%d, %.17g, %.17g, %.17g\n", 
		dataset.gt.beacons(i).id,
		dataset.gt.beacons(i).coordinates(1),
		dataset.gt.beacons(i).coordinates(2),
		dataset.gt.beacons(i).coordinates(3));

endfor;

fclose(fid);

