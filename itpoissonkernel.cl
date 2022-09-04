__kernel void itpoisson(__global float* data, __global float* mapdata, __global float* rhomap, __global float* outdata) 
{
	size_t gidx = get_global_id(0);
	size_t gidy = get_global_id(1);
	size_t gxsize = get_global_size(0);
	size_t gysize = get_global_size(1);

	size_t gid = (get_global_id(1) * get_global_size(0)) + get_global_id(0);
	
	if(mapdata[gid] == 0)
	{
		outdata[gid] =  data[gid] + 1.0 * 	(- data[gid] * rhomap[gid]		 
												+	(
														2.0 * 	(	rhomap[gid + gxsize] 	* data[gid + gxsize] 	+ rhomap[gid - gxsize] 		* data[gid - gxsize] 
																+ 	rhomap[gid + 1] 		* data[gid + 1] 	 	+ rhomap[gid - 1] 			* data[gid - 1]
																)
																+ 	rhomap[gid-1 + gxsize] 	* data[gid-1 + gxsize] 	+ rhomap[gid+1 + gxsize] 	* data[gid+1 + gxsize] 
																+ 	rhomap[gid-1 - gxsize] 	* data[gid-1 - gxsize] 	+ rhomap[gid+1 - gxsize] 	* data[gid+1 - gxsize]
													)/12.0
																						);
	} else {
		outdata[gid] = data[gid];
	}
}