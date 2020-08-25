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
#include "wdm.h"
#include "chead.h"

// gcc -o match match.c csubs.c -lm -lgsl

int main(int argc, char *argv[])
{
  int i, j, k;
  char filename[1024];
  double Tobs, Tchunk, df;
  double x, y, z;
    double *D;
  int N, NC, M, up;
    int Nx, Mx;
    double f, t;
   double **data1, **data2;
    double DT, DF;
    int *tim;
    

  FILE *in;
  FILE *ifp;
  FILE *out;
    
    if(argc<3)
    {
        printf("./match data1 data2\n");
        return 1;
    }
    
    
    N = Nt*Nf;
    
    Tobs = dt*(double)(N);   // duration (~ 1 year)
    df = 1.0/Tobs;
    
    printf("Tobs = %e  %d\n", Tobs, N);
    
    DT = dt*(double)(Nf);           // width of wavelet pixel in time
    DF = 1.0/(2.0*dt*(double)(Nf));   // width of wavelet pixel in frequency

    
    data1 = double_matrix(Nt,Nf);
    data2 = double_matrix(Nt,Nf);

    in = fopen(argv[1],"r");
    for(i=0; i< Nf; i++)
        {
        for(j=0; j< Nt; j++)
         {
            fscanf(in,"%lf%lf%lf", &x, &y, &data1[j][i]);
         }
           fscanf(in,"\n");
        }
    fclose(in);
    
    in = fopen(argv[2],"r");
    for(i=0; i< Nf; i++)
        {
        for(j=0; j< Nt; j++)
         {
            fscanf(in,"%lf%lf%lf", &x, &y, &data2[j][i]);
         }
           fscanf(in,"\n");
        }
    fclose(in);
    
    out = fopen("matchT.dat","w");
     for(j=mult; j< (Nt-mult); j++)
      {
          x = 0.0;
          y = 0.0;
          z = 0.0;
      for(i=0; i< Nf; i++)
      {
         x += data1[j][i]*data1[j][i];
         y += data2[j][i]*data2[j][i];
         z += data1[j][i]*data2[j][i];
      }
          fprintf(out,"%d %f %e\n", j, z/sqrt(x*y), x);
      }
    fclose(out);

    
    x = 0.0;
    y = 0.0;
    z = 0.0;

     for(j=mult; j< (Nt-mult); j++)
     {
     for(i=0; i< Nf; i++)
     {
        x += data1[j][i]*data1[j][i];
        y += data2[j][i]*data2[j][i];
        z += data1[j][i]*data2[j][i];
     }
     }
    
    printf("%e %e %.12f %e\n", sqrt(x), sqrt(y), z/sqrt(x*y), 1.0-z/sqrt(x*y));
    

    return 0;

}

