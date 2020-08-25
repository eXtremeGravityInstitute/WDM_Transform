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

// gcc -o transform_time transform_time.c csubs.c -lm -lgsl

int main(int argc, char *argv[])
{
   int i, j, k, jj;
   int n, m;
   char filename[1024];
   double Tobs, Tchunk, df;
   double DT, DF;
   double x, y, z;
    double c, s;
    double xx, yy;
    double A, B, DOM, OM;
    double insDOM;
    double om, fac;
    double f0;
    double *R;
    double *DX;
    double *data, *time;
    double *wdata;
    double **wave;
    
    int K, ND, NC, M,  L, up;
    int Nx, Mx;
    double f, t, T, dom, nrm;
    double **data1, **data2;
    double *hist;
    double *phi;
    int *tim;
    
    clock_t start, end;
    double cpu_time_used;

  FILE *in;
  FILE *ifp;
  FILE *out;

    
    if(argc<2)
    {
        printf("./transform_time filename\n");
        return 1;
    }
    
    in = fopen(argv[1],"r");
    
    ND = Nt*Nf;
    
        
    // the time domain data stream
    data = (double*)malloc(sizeof(double)* (ND));
    time = (double*)malloc(sizeof(double)* (ND));
    
    for(i=0; i< ND; i++) fscanf(in,"%lf%lf", &time[i], &data[i]);
  
    DT = dt*(double)(Nf);           // width of wavelet pixel in time
    DF = 1.0/(2.0*dt*(double)(Nf)); // width of wavelet pixel in frequency
    
    OM = PI/dt;
    
    M = Nf;
    
    L = 2*M;
    
    DOM = OM/(double)(M);
    
    insDOM = 1.0/sqrt(DOM);
    
    B = OM/(double)(L);
    
    A = (DOM-B)/2.0;
    
    K = mult*2*M;
    
    T = dt*(double)(K);
    
    printf("filter length = %d\n", K);
    printf("Pixel size DT (seconds) %e DF (Hz) %e\n", DT, DF);
    printf("Filter length (seconds) %e full filter bandwidth %e\n", T, (A+B)/PI);
    
    dom = TPI/T;  // max frequency is K/2*dom = pi/dt = OM

    DX = (double*)malloc(sizeof(double)* (2*K));

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
          fclose(out);
    
    nrm = sqrt((double)(K)/dom);
  
    
    // windowed data packets
    wdata = (double*)malloc(sizeof(double)* (K));
    
    wave = double_matrix(Nt,Nf);  // wavelet wavepacket transform of the signal
    
    fac = sqrt(2.0)/nrm;
    
    for(i=0; i < K; i++) phi[i] *= fac;
    
     start = clock();
    for(i=0; i< Nt; i++)
     {
         
        for(j=0; j< K; j++)
        {
            jj = i*Nf-K/2+j;
            if(jj < 0) jj += ND;  // periodically wrap the data
            if(jj >= ND) jj -= ND; // periodically wrap the data
            wdata[j] = data[jj]*phi[j];  // apply the window
        }
         
        gsl_fft_real_radix2_transform(wdata, 1, K);
         
         wave[i][0] = sqrt(2.0)*wdata[0];
        
         for(j=1; j< Nf; j++)
         {
            // printf("%d %d\n", i, j);
             if((i+j)%2 ==0)
             {
                 wave[i][j] = wdata[j*mult];
             }
             else
             {
                wave[i][j] = -wdata[K-j*mult];
             }
         }
         
       }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("The transform took %f seconds\n", cpu_time_used);
    
    out = fopen("BinaryT.dat","w");

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
    
    printf("mean = %e sigma = %e\n", xx, sqrt(yy-xx*xx));
    
    printf("max %e %d %d\n", y, n, m);
    
    x = pow(10.0,floor(log10(y)));
    
    z = y/x;
    z = ceil(z)*x;
    
    printf("%e\n", z);
    
    
    out = fopen("tran.gnu","w");
    fprintf(out,"set term png enhanced truecolor crop font Helvetica 18  size 1200,800\n");
    fprintf(out,"set output 'tran.png'\n");
    fprintf(out,"set pm3d map corners2color c1\n");
   // fprintf(out,"set ylabel 'f (Hz)'\n");
    //fprintf(out,"set xlabel 't (s)'\n");
      fprintf(out,"set yrange [0:200]\n");
      fprintf(out,"set ylabel 'frequency'\n");
      fprintf(out,"set xlabel 'time'\n");
    fprintf(out,"set cbrange [%e:%e]\n", -z, z);
    fprintf(out,"set palette defined (0 '#b2182b', 1 '#ef8a62', 2 '#fddbc7', 3 '#ffffff', 4 '#d1e5f0', 5 '#67a9cf', 6 '#2166ac')\n");
    fprintf(out,"splot 'BinaryT.dat' using 1:2:3 notitle\n");
    fclose(out);
    
    return 0;

}

