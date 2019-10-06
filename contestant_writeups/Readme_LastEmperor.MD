A high-level description.

1) --- Delete all zero values and delete outliers in at the beginning and the end dataset by calculating the ratio of two consecutive standing values.
2) --- Check ratio between mean of values and median, if ratio more k than:
 	a) Split all values on n = shape(values) / 2 bins and 		check is there a gap where there are no values in the bin 	in last ? values.
		-  if yes that gap will be point of split values on 		the cluster, then we just count the medians of the 		clusters
 		- else use kmeans
else:
	b)Split all values on n = shape(values) / 2 bins and 		check is there a gap where there are no values in the bin 	in interval between 15 and 85 percentile.
		-  if yes that gap will be point of split values on 		the cluster, then we just count the medians of the 		clusters
		 - else use kmeans

All computation split between processors.
After this calculation, for gens with the same barcode we calculate medians and remember which value is greater. And for all values for this gens we swap values if they stand in the wrong order. After this we clip all values for all gens that more median in 2 times.