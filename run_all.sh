#!/bin/bash

gcc -o coefficientsWDM_time coefficientsWDM_time.c csubs.c -lm -lgsl
gcc -o coefficientsWDM_freq coefficientsWDM_freq.c csubs.c -lm -lgsl
gcc -o Chirp_WDM Chirp_WDM.c -lm -lgsl
gcc -o transform_time transform_time.c csubs.c-lm -lgsl
gcc -o transform_freq transform_freq.c csubs.c -lm -lgsl
gcc -o match match.c csubs.c -lm -lgsl

mkdir coeffs

#compute the fast frequency Taylor expansion coefficients
./coefficientsWDM_freq
#compute the fast time Taylor expansion coefficients
./coefficientsWDM_time
#generate a time domain chirplet and its fast WDM transform four ways
./Chirp_WDM
#compute the WDM transform of the chirplet directly in the time domain
./transform_time chrp_time.dat
#compute the WDM transform of the chirplet by FFTing first
./transform_freq chrp_time.dat 0
#compute the match between the dirtect time and frequency domain transforms
./match BinaryT.dat BinaryF.dat
#compute the match between the dirtect frequency domain transform and each of
#the fast transforms
./match BinaryF.dat BinaryTaylorT.dat
./match BinaryF.dat BinarySparseT.dat
./match BinaryF.dat BinaryTaylorF.dat
./match BinaryF.dat BinarySparseF.dat

gnuplot tran.gnu
open tran.png

