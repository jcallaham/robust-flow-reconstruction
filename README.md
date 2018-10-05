# Robust flow field reconstruction from limited measurements

Code for data analysis and figures from J. Callaham, K. Maeda, and S. L. Brunton (2018)

### Obtaining data

With the exception of the mixing layer data (which was generated using a DoD resource), all data is publicly available.
* Re=100 cylinder flow data available from [DMDbook.com](www.DMDbook.com/DATA.zip)
* NOAA OISST v2 data available [here](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html) or can be downloaded automatically with the script in noaa-sst\remote_read.m
* Same for the HYCOM Gulf of Mexico data - it's online [here](https://hycom.org/hycom), but can be downloaded and formatted automatically with the script in hycom\remote_read.m

### Required packages

All code available here was written for MATLAB R2018a and relies heavily on the freely available [cvx](http://cvxr.com/cvx/) package. Some parts of the code also use Ron Rubenstein's implementations of the K-SVD and OMP algorithms, available [here](http://www.cs.technion.ac.il/~ronrubin/software.html)
