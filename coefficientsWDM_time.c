/*******************************************************************************************

Copyright (c) 2020 Neil Cornish

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Constants.h"
#include "chead.h"
#include "wdm.h"

// gcc -o coefficientsWDM_time coefficientsWDM_time.c csubs.c -lm -lgsl

int main()
{
  int i, j, k, n, K;
  int Nstep;
  char filename[1024];
  double Tfilt, Tchunk, df;
  double DT, DF;
  double x, y, z, xa, xb, dx;
    int N, NC, M, NK, up;
    int Nx, NT;
    int n0, n1, n2, n3;
    double f, t;
    double c, s, delf;
    double phase, freq;
    double m1, m2, Mt, Mc, eta, tc, Theta, fisco, fdot;
    double tx, tmid;
    double fcent, delt;
    double odc, ods, evc, evs;
    double *fd;
    double OM, DOM, A, B, BW;
    double insDOM;
    double om, fac, dom, nrm;
    double *phi, *DX;
    
    int *Nfsam;
    
    int ii, jj;
    
    double *wave;
    
    int prime;
    char ch;
    
  FILE *in;
  FILE *fout;
  FILE *out;
    
     DT = dt*(double)(Nf);           // width of wavelet pixel in time
     DF = 1.0/(2.0*dt*(double)(Nf));   // width of wavelet pixel in frequency
     
     OM = PI/dt;
     
     DOM = OM/(double)(Nf);
     
     insDOM = 1.0/sqrt(DOM);
     
     B = OM/(double)(2*Nf);
     
     A = (DOM-B)/2.0;
     
     K = mult*2*Nf;
     
     Tfilt = dt*(double)(K);
     
     printf("Filter length (seconds) %e\n", Tfilt);
     
     dom = TPI/Tfilt;  // max frequency is K/2*dom = pi/dt = OM
     
     printf("full filter bandwidth %e  samples %d\n", (A+B)/PI, (int)(((A+B)/PI)*Tfilt));

     DX = (double*)malloc(sizeof(double)* (2*K));

     wave = (double*)malloc(sizeof(double)* (K));
     
     //zero frequency
     REAL(DX,0) =  insDOM;
     IMAG(DX,0) =  0.0;
     
     // postive frequencies
     for(i=1; i<= K/2; i++)
     {
         om = (double)(i)*dom;
         z = phitilde(om, insDOM, A, B);
         REAL(DX,i) =  z;
         IMAG(DX,i) =  0.0;
     }
     
     // negative frequencies
     for(i=1; i< K/2; i++)
     {
         om = -(double)(i)*dom;
         z = phitilde(om, insDOM, A, B);
         REAL(DX,K-i) =  z;
         IMAG(DX,K-i) =  0.0;
     }
     
     gsl_fft_complex_radix2_backward(DX, 1, K);
     
     
     phi = (double*)malloc(sizeof(double)* (K));
     
            for(i=0; i < K/2; i++)
             {
              phi[i] = REAL(DX,K/2+i);
             }
            for(i=0; i< K/2; i++)
              {
              phi[K/2+i] = REAL(DX,i);
              }

     
     nrm = 0.0;
     for(i=0; i < K; i++) nrm += phi[i]*phi[i];
     nrm = sqrt(nrm);
    // printf("phi norm %e\n", nrm);
    
    // it turns out that all the wavelet layers are the same modulo a
    // shift in the reference frequency. Just have to do a single layer
    // we pick one far from the boundaries to avoid edge effects
    
    k = Nf/16;
    
    wavelet(k, wave, K, nrm, dom, DOM, A, B, insDOM);
    
    // total width of wavelet in frequency
    BW = (A+B)/PI;
    
    fd = double_vector(Nfd);
    
    fd[0] = 0.0;
    fd[1] = DF/Tfilt*dfdot;  // sets the f-dot increment
    for(j=2; j< Nfd; j++) fd[j] = (double)(j)*fd[1];
    
    printf("%e %.14e %.14e %e %e\n", DT, DF, DOM/TPI, fd[1], fd[Nfd-1]);
    
    // frequency spacing
    df = BW/(double)(Nsf);
    
    Nfsam = int_vector(Nfd);
    
    for(j=0; j< Nfd; j++)
    {
     Nfsam[j] = (int)((BW+fd[j]*Tfilt)/df);
     if(Nfsam[j]%2 != 0) Nfsam[j]++; // makes sure it is an even number
    }

        // The odd wavelets coefficienst can be obtained from the even.
        // odd cosine = -even sine, odd sine = even cosine

        // each wavelet covers a frequency band of width DW
        // execept for the first and last wavelets
        // there is some overlap. The wavelet pixels are of width
        // DOM/PI, except for the first and last which have width
        // half that
    
        fcent = (double)(k)*DF;
            
            for(jj=0; jj< Nfd; jj++)  // loop over f-dot slices
            {
                
            printf("%d %d\n", jj, Nfsam[jj]);
                
            sprintf(filename, "coeffs/WDMcoeffs%d.dat", jj);
            out = fopen(filename,"w");
                
            for(j=0; j< Nfsam[jj]; j++)  // loop over frequency slices
            {
                f = fcent+((double)(j-Nfsam[jj]/2)+0.5)*df;

                evc = 0.0;
                evs = 0.0;
                
                for(i=0; i< K; i++)
                {
                    t = ((double)(i-K/2))*dt;
                    z = TPI*(f*t+0.5*fd[jj]*t*t);
                    c = cos(z);
                    s = sin(z);
                    evc += wave[i]*c;
                    evs += wave[i]*s;
                }
                
                fprintf(out,"%d %.14e %.14e\n", j, evc, evs);
                
            }
                
                fclose(out);
                
            }
    
   
       
   free_int_vector(Nfsam);
   free_double_vector(wave);
   free_double_vector(fd);
   free_double_vector(phi);
   free_double_vector(DX);
    
    return 0;

}

