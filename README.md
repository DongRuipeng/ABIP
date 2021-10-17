## R example 

The R examples are in the folder "ABIP-beta/test_CLIME.R", which requires the abip.dll file attached in the folder "ABIP-beta/src/abip.dll". 

Although there has attached the dll file, we maybe need to recompile it when downloading the code because of the difference of the system/devices. 

### Compiling method 

The compiling needs to install Rtools (https://cran.r-project.org/bin/windows/Rtools/) on Windows system. In the terminal of the local system under the path "ABIP-beta/src", running "make main" will generate the abip.dll file. 

## Matlab example 

The matlab examples are in the folder "Original-ABIP/abip/test_CLIME.m", which requires the abip_direct.mexw64 file attached in the same folder. 

Similar to the dll file, we maybe need to recompile abip_direct.mex file because of the difference of the system/devices. 

### Compiling method 

Similar to Rtools, we need to install the minGW64 plug-in tools in Matlab (https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler). 

After installing it, directly running the script "Original-ABIP/abip/make_abip.m" will generate the abip_direct.mex file. 