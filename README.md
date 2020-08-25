# WDM_Transform
Codes to compute the WDM wavelet transform

The directory contains several codes and two run scripts. The scripts were developed for OS X and gcc, and may need to be modified for other hardwares.

If you want to run everything, including the generation of lookup tables, type

source run_all.sh

If the lookup tables have already been generated, type

source run.sh

The run_all.sh script compiles all the codes then runs the codes coefficientsWDM_time.c and coefficientsWDM_freq.c to generate the lookup tables used for the fast Taylor expansion based WDM transforms. The time domain code takes a lot longer to run than the frequency domain code. These only need to be run once for a given set of WDM parameters.
The WDM parameters are set in the header file wdm.h. 

If the lookup tables have already been generated, that step can be skipped. The run.sh script is used when the lookup tables already exits.

The scripts next call the code Chirp_WDM.c. This code illustrates the four fast methods for generating the WDM transform for a binary signal using a chirplet waveform with a linear chirp. The code first computes the chirplet waveform in the time domain and saves it to the file chrp_time.dat. It then computes the WDM transform using the the TaylorT, TaylorF, and sparse T and F transforms. These are saved to the files BinaryTaylorT.dat, BinarySparseT.dat, BinaryTaylorF.dat, BinarySparseF.dat. Here the term "Binary" refers to the time-frequency grid being rectangular, as opposed to "Dyadic", which is more commonly used for discrete wavelet transforms.

The scripts next call transform_time.c and transform_freq.c which are used to transform the time domain chirplet into the WDM basis using either the direct time domain transform or the direct frequency domain transform rolling an FFT of the data (with a Tukey window applied to limit leakage in the FFT).

Next the scripts call the code match.c, which is used to compute the match between the time domain and frequency domain versions of the transform. This is followed by calling the match code to compare the direct frequency domain transform to the various fast transforms.

The direct transform codes also write gnuplot scripts, tran.gnu and tranf.gnu, that can be used to produce images of the WDM transform.


