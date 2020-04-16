@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2014a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2014a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=triangular_algorithm_TBM_mex
set MEX_NAME=triangular_algorithm_TBM_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for triangular_algorithm_TBM > triangular_algorithm_TBM_mex.mki
echo COMPILER=%COMPILER%>> triangular_algorithm_TBM_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> triangular_algorithm_TBM_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> triangular_algorithm_TBM_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> triangular_algorithm_TBM_mex.mki
echo LINKER=%LINKER%>> triangular_algorithm_TBM_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> triangular_algorithm_TBM_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> triangular_algorithm_TBM_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> triangular_algorithm_TBM_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> triangular_algorithm_TBM_mex.mki
echo BORLAND=%BORLAND%>> triangular_algorithm_TBM_mex.mki
echo OMPFLAGS= >> triangular_algorithm_TBM_mex.mki
echo OMPLINKFLAGS= >> triangular_algorithm_TBM_mex.mki
echo EMC_COMPILER=msvcsdk>> triangular_algorithm_TBM_mex.mki
echo EMC_CONFIG=optim>> triangular_algorithm_TBM_mex.mki
"C:\Program Files\MATLAB\R2014a\bin\win64\gmake" -B -f triangular_algorithm_TBM_mex.mk
