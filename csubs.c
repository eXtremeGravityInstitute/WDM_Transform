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

#include "Constants.h"
#include "chead.h"
#include "wdm.h"

void wavelet(int m, double *wave, int N, double nrm, double dom, double DOM, double A, double B, double insDOM)
{
    
    int i;
    double om;
    double x, y, z;
    double *DE;
    
     DE = (double*)malloc(sizeof(double)* (2*N));

         // zero and postive frequencies
          for(i=0; i<= N/2; i++)
           {
            om = (double)(i)*dom;
            
            y = phitilde(om+(double)(m)*DOM, insDOM, A, B);
            z = phitilde(om-(double)(m)*DOM, insDOM, A, B);
               
               x = y+z;
            
               REAL(DE,i) = IRTWO*x;
               IMAG(DE,i) = 0.0;

           }
     
           // negative frequencies
            for(i=1; i< N/2; i++)
            {
             om = -(double)(i)*dom;
                         
              y = phitilde(om+(double)(m)*DOM, insDOM, A, B);
              z = phitilde(om-(double)(m)*DOM, insDOM, A, B);
                         
              x = y+z;
    
                REAL(DE,N-i) = IRTWO*x;
                IMAG(DE,N-i) = 0.0;
             
             }
    
          gsl_fft_complex_radix2_backward(DE, 1, N);
        
        
                 for(i=0; i < N/2; i++)
                    {
                        wave[i] = REAL(DE,N/2+i)/nrm;
                    }
                 for(i=0; i< N/2; i++)
                 {
                      wave[i+N/2] = REAL(DE,i)/nrm;
                 }

      free(DE);
    
}


double phitilde(double om, double insDOM, double A, double B)
{
    double x, y, z;
    
       z = 0.0;
       
       if(fabs(om) >= A && fabs(om) < A+B)
       {
           x = (fabs(om)-A)/B;
           y = gsl_sf_beta_inc(nx, nx, x);
           z = insDOM*cos(y*PI/2.0);
       }
       
       if(fabs(om) < A) z = insDOM;
    
    return(z);
    
}

void tukey(double *data, double alpha, int N)
{
  int i, imin, imax;
  double filter;
  
  imin = (int)(alpha*(double)(N-1)/2.0);
  imax = (int)((double)(N-1)*(1.0-alpha/2.0));
  
    int Nwin = N-imax;

    for(i=0; i< N; i++)
  {
    filter = 1.0;
    if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
    if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
    data[i] *= filter;
  }
  
}


int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}


double **double_matrix(int N, int M)
{
    int i;
    double **m = malloc( (N+1) * sizeof(double *));
    
    for(i=0; i<N+1; i++)
    {
        m[i] = malloc( (M+1) * sizeof(double));
    }
    
    return m;
}

void free_double_matrix(double **m, int N)
{
    int i;
    for(i=0; i<N+1; i++) free_double_vector(m[i]);
    free(m);
}

double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
{
    free(v);
}



