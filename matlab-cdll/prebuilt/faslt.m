% faslt.m Help file for faslt.mexw64
%
%	Author: Harald Bârzan (Transylvanian Institute of Neuroscience)
% 	Date:	June 2020
%
%	Performs fractional adaptive superlet transform (FASLT)	on input data.
%	The input matrix contains the input trials, row by row.
%	If multiple trials are provided, the result is the average spectrum.
%
%
%	Parameter Name				- Parameter Type							- Description
%
%  		input_data				- Scalar Matrix								- Input buffers (row major) - each row is a trial
%		sampling_rate			- Scalar Number								- The sampling frequency in Hz
%		frequency_interval		- Scalar Vector								- tuple (vector of size 2) containing the lowest and highest frequency
%		frequency_count			- Integer Number							- the number of frequency bins in the interval
%		cycle_count				- Scalar Number								- the number of cycles of the shortest wavelets
%		superresolution_order	- Scalar Vector								- tuple containing the lowest and the highest superresolution orders
%		multiplicative*			- Scalar Number								- 0 to use additive superresolution, multiplicative otherwise (default: true)
%		fractional*				- Scalar Number								- 0 to use integral ASLT, uses fractional (FASLT) otherwise (default: true)
%
%	Return values
%		S						- Scalar Matrix								- A matrix of size length(frequencies) x size(input_data, 2) containing the result
%
%	Parameters marked with asterisk (*) are optional.
%
%
%	Example usage:
%
%		n 	= 1000;											% signal size	
%		fs 	= 1000;											% sampling rate
% 		t 	= (1 : n) / fs;									% time dimension
%		sig = sin(2 * pi * 50 * t);							% sine at 50 Hz
%				
%		S1 	= faslt(sig, [10 90], 81, 3, [10 10]);			% perform SLT in range 10-90Hz (81 points - 1 point per Hz) with c1=3 and order 10 
%		S2	= faslt(sig, [10 90], 161, 5, [5 20]);			% perform ASLT in range 10-90Hz (161 points - 2 points per Hz) with c1=5 and order between 5 and 20
%
%	WARNING: EACH ROW IS A TRIAL. THIS METHOD WILL FAIL FOR COLUMN VECTORS.
%
