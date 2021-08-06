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

3. Pyure Python implementation ([superlet.py](./python/superlet.py))
   - multiplicative and fractional adaptive SLT
   - needs only minimal dependencies as defined in [environment.yml](./python/environment.yml), standard for any scientific Python environment
   - when run as a script from the command line via `python superlet.py` will produce the example output [synthetic_example.png](./python/synthetic_example.png)
   - can be imported as stand-alone module via `from superlet import superlet`