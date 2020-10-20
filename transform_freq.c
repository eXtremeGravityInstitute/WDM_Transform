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

#include <time.h>

// gcc -o transform_freq transform_freq.c csubs.c -lm -lgsl

int main(int argc, char *argv[])
{
  int i, j, k, l, jj, tf;
    int n, m;
  char filename[1024];
  double Tobs, Tchunk, df;
    double DT, DF;
    double x, y, z, alpha;
    double c, s;
    double xx, yy;
    double A, B, DOM, OM;
    double insDOM;
    double om, fac;
    double f0;
    double *R;
    double *DX;
    double *data, *time, *freq;
    double *wdata;
    double **wave;
    
    int K, ND, NC, M,  L, up;
    int Nx, Mx;
    double f, t, T, dom, nrm;
    double **data1, **data2;
    double *hist;
    double *phif;
    int *tim;
    
    clock_t start, end;
    double cpu_time_used;

  FILE *in;
  FILE *ifp;
  FILE *out;
    
    if(argc<3)
    {
        printf("./transform_freq filename time/freq\n");
        return 1;
    }
    
    in = fopen(argv[1],"r");
    
    tf = atoi(argv[2]); // 0 for TD, 1 for FD
    
    ND = Nt*Nf;
    Tobs = dt*(double)(ND);
    
    if(tf == 0)
    {
        // the data stream
    data = (double*)malloc(sizeof(double)* (ND));
    time = (double*)malloc(sizeof(double)* (ND));
    for(i=0; i< ND; i++) fscanf(in,"%lf%lf", &time[i], &data[i]);
    }
    else
    {
        // the data stream
        data = (double*)malloc(sizeof(double)* (ND));
        freq = (double*)malloc(sizeof(double)* (ND/2));
        for(i=1; i< ND/2; i++) fscanf(in,"%lf%lf%lf", &freq[i], &data[i], &data[ND-i]);
    }
  
    DT = dt*(double)(Nf);           // width of wavelet pixel in time
    DF = 1.0/(2.0*dt*(double)(Nf)); // width of wavelet pixel in frequency
    
    OM = PI/dt;
    
    L = 2*Nf;
    
    DOM = OM/(double)(Nf);
    
    insDOM = 1.0/sqrt(DOM);
    
    B = OM/(double)(L);
    
    A = (DOM-B)/2.0;
    
    Tobs = dt*(double)(ND);
    
    dom = TPI/Tobs;
    
    printf("Pixel size DT (seconds) %e DF (Hz) %e\n", DT, DF);
    printf("full filter bandwidth %e\n", (A+B)/PI);
    
    phif = (double*)malloc(sizeof(double)* (Nt/2+1));

    for(i=0; i<= Nt/2; i++)
    {
        om = (double)(i)*dom;
        phif[i] = phitilde(om, insDOM, A, B);
    }

    
     nrm = 0.0;
     for(l=-Nt/2; l<= Nt/2; l++) nrm += phif[abs(l)]*phif[abs(l)];
     nrm = sqrt(nrm/2.0);
     nrm *= (double)(Nt);
    
     for(j=0; j<= Nt/2; j++) phif[j] /= nrm;

    wave = double_matrix(Nt,Nf);  // wavelet wavepacket transform of the signal
    
    DX = double_vector(2*Nt);
    
    // Window the data and FFT
    // Tukey window parameter. Flat for (1-alpha) of data
    alpha = (2.0*(4.0*DT)/Tobs);
    
    start = clock();
    
    if(tf == 0)
    {
    tukey(data, alpha, ND);
    gsl_fft_real_radix2_transform(data, 1, ND);
    }

    for(m=0; m< Nf; m++)
     {
         for(i=0; i< Nt; i++)
         {
            REAL(DX,i) = 0.0;
            IMAG(DX,i) = 0.0;
         }
         
        for(j=-Nt/2+1; j< Nt/2; j++)
        {
            jj = j + m*Nt/2;
            
            if(jj > 0 && jj < ND/2)
            {
                if(j >= 0)
                {
                REAL(DX,j) = data[jj]*phif[abs(j)];
                IMAG(DX,j) = data[ND-jj]*phif[abs(j)];
                }
                else
                {
                 REAL(DX,Nt+j) = data[jj]*phif[abs(j)];
                 IMAG(DX,Nt+j) = data[ND-jj]*phif[abs(j)];
                }
            }
        }
         
         gsl_fft_complex_radix2_backward(DX, 1, Nt);
         
        
         for(n=0; n < Nt; n++)
         {
             x = -1.0;
             if(m%2==0 || n%2==0) x = 1.0;
             
          
                 if((n+m)%2 ==0)
                 {
                      wave[n][m] = x*REAL(DX,n);
                 }
                 else
                 {
                    wave[n][m] = -x*IMAG(DX,n);
                 }
            
         }
         
       }
    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("The transform took %f seconds\n", cpu_time_used);
    
    out = fopen("BinaryF.dat","w");

    xx = 0.0;
    yy = 0.0;
    
      x = 0.0;
    y = 0.0;
       for(i=0; i< Nf; i++)
       {
          for(j=0; j< Nt; j++)
           {
           //fprintf(out, "%e %e %.14e\n", (double)(j)*DT, (double)(i)*DF, wave[j][i]);
            fprintf(out, "%d %d %.14e\n", j, i, wave[j][i]);
               z = fabs(wave[j][i]);
               if(z > y)
               {
                   n = j;
                   m = i;
                   y = z;
               }
               
               xx += wave[j][i];
               yy += wave[j][i]*wave[j][i];
               
           }
               
       fprintf(out, "\n");
               
     }
    fclose(out);
    
    xx /= (double)(Nf*Nt);
    yy /= (double)(Nf*Nt);
    
    printf("mean = %e sigma = %e sq = %e\n", xx, sqrt(yy-xx*xx), yy);
    
    printf("max %e %d %d\n", y, n, m);
    
    x = pow(10.0,floor(log10(y)));
    
    z = y/x;
    z = ceil(z)*x;
    
    printf("%e\n", z);
    
    
    out = fopen("tranf.gnu","w");
    fprintf(out,"set term png enhanced truecolor crop font Helvetica 18  size 1200,800\n");
    fprintf(out,"set output 'tranf.png'\n");
    fprintf(out,"set pm3d map corners2color c1\n");
   // fprintf(out,"set ylabel 'f (Hz)'\n");
    //fprintf(out,"set xlabel 't (s)'\n");
     fprintf(out,"set yrange [0:200]\n");
      fprintf(out,"set ylabel 'frequency'\n");
      fprintf(out,"set xlabel 'time'\n");
    fprintf(out,"set cbrange [%e:%e]\n", -z, z);
    fprintf(out,"set palette defined (0 '#b2182b', 1 '#ef8a62', 2 '#fddbc7', 3 '#ffffff', 4 '#d1e5f0', 5 '#67a9cf', 6 '#2166ac')\n");
    fprintf(out,"splot 'BinaryF.dat' using 1:2:3 notitle\n");
    fclose(out);
    
    return 0;

}


