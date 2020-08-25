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

// gcc -o coefficientsWDM_freq coefficientsWDM_freq.c csubs.c -lm -lgsl

int main()
{
  int i, j, k, n, m, K;
  int Nstep, Np;
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
    double fcent, delt, Tobs;
    double odc, ods, evc, evs;
    double *td;
    double OM, DOM, A, B, BW;
    double insDOM;
    double om, fac, dom, nrm;
    double *phi;
    
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
    
     Tobs = dt*(double)(Nt*Nf);
     
     printf("Filter length (seconds) %e\n", Tfilt);
    
    // total width of wavelet in frequency
    BW = (A+B)/PI;
    
    // The wavelet is computed at the sample frequency 1/Tobs
    Np = (int)(BW*Tobs);
    
    dom = 2.0*(A+B)/(double)(Np);
    
    phi = (double*)malloc(sizeof(double)* (Np));
    
    for(i=0; i < Np; i++)
    {
     f = -BW/2.0+(double)(i)/Tobs;
     phi[i] = phitilde(TPI*f, insDOM, A, B);
    }
    
       nrm = 0.0;
       for(i=0; i < Np; i++) nrm += phi[i]*phi[i];
       nrm = sqrt(nrm);
       nrm *= (sqrt(2.0)*(double)(Nf)*dt);
    
    // time spacing
    delt = (Tfilt/2.0)/(double)(Nst);
   
    td = double_vector(Ntd);
    
    //printf("%e %e\n", Tfilt/DF, 1.0/(DF*DF));
    
    td[0] = 0.0;
    td[1] = dtdf/(DF*DF);
    for(j=2; j< Ntd; j++) td[j] = (double)(j)*td[1];
    
    Nfsam = int_vector(Ntd);
    
    for(j=0; j< Ntd; j++)
    {
     Nfsam[j] = (int)((Tfilt/2.0+0.5*td[j]*DF)/delt);
     
    }
    
    for(k=0; k< Ntd; k++)
    {
        printf("%d %d\n", k, Nfsam[k]);
        sprintf(filename, "coeffs/WDMcoeffsf%d.dat", k);
        out = fopen(filename,"w");
    for(j=-Nfsam[k]; j < Nfsam[k]; j++)
    {
        t = delt*(double)(j);
        x = 0.0;
        y = 0.0;
           for(i=0; i < Np; i++)
           {
               f = -BW/2.0+(double)(i)/Tobs;
               c = cos(2.0*PI*(t*f+0.5*td[k]*f*f));
               s = sin(2.0*PI*(t*f+0.5*td[k]*f*f));
               x += c*phi[i]/nrm;
               y += s*phi[i]/nrm;
           }
      fprintf(out,"%d %e %e\n", j, x, y);
    }
        fclose(out);
    }
       
   free_int_vector(Nfsam);
   free_double_vector(td);
   free_double_vector(phi);
  

    
    return 0;

}


