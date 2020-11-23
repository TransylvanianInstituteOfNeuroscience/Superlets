# Superlets
Superlet Transform (SLT) code for MATLAB.

Two versions are supplied:
1. Pure MATLAB implementation (`aslt.m`, `faslt.m`)
	- Can be run completely in MATLAB, no further steps required.
	- The `matlab-pure` folder also contains a script `Superlet_Toy_Data.m` which produces a test signal that proves the utility of SLT.

2. MATLAB interface for C++ backend (`faslt.mex64`)
	- The folder contains a Visual Studio 2019 (v16.0) solution with the necessary files to build the mex64 file
	- The mex64 file will only work with a 64-bit MATLAB (who uses 32-bit for data analysis anyway?)
	- The project file must be modified to include two sets of libraries: 
		- the MATLAB C libs - usually in your MATLAB folder `($(MATLABPath)\extern\lib\win64\microsoft)` - the last directory might be different for other operating systems
		- the FFTW libs - a `.lib` file is included in 'matlab-cdll/lib' folder. The required `.dll` file may be downloaded from the FFTW website. It must be put in the same folder with the `.mex64` file in order for the script to work.