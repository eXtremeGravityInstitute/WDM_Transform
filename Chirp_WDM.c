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
#include <string.h>
#include <time.h>

#include "Constants.h"
#include "Header.h"
#include "wdm.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>


// gcc -o Chirp_WDM Chirp_WDM.c -lm -lgsl


int main()
{
    int i, j, k, ii, jj, kk, mc, q,  n, m, mq, l, dn, kx;
    char filename[1024];
    double df, tmid;
    double x, y, z, dx;
    double *D;
    int NC, M, NK, NM, NU, NMF, up;
    int Nx, Nstep, prime;
    int n0, n1, n2, n3;
    double mx, f, ff, t, alpha;
    double *c0, *c1, *c2, *c3;
    double *s0, *s1, *s2, *s3;
    double *fa;
    double c, s, delf;
    double *data0, *data1, *data2, *data3;
    double phase, freq, Tw;
    double *wv, *wvs;
    double m1, m2, Mt, Mc, eta, tc, Theta, fisco, fdot, fddot, DL;
    double ***lookup;
    double **wave;
    double *waveT, *waveF;
    double *params;
    double *h;
    double fac;
    int *Tlist, *Flist;
    int *list;
    double *Amp, *Phase;
    int *Nfsam;
    double *fd;
    double tau, gamma, tp, fp, Phi0, A0;
    double *phif, *phi, *DX;
    
    double *flayer;
    
    clock_t start, end;
    double cpu_time_used;
    
    WDM *wdmpars = malloc(sizeof(WDM));

    FILE *in;
    FILE *infile[Nfd];
    FILE *out;
    
    // compute all the quantities we need to describe the WDM wavelets
    create_WDM(wdmpars);

    // Want gamma*tau large so that the SPA is accurate
    // Pick fdot so that both the Taylor expanded time and frequency domain transforms are valid
    // => fdot < 8 DF/T_w and fdot > DF^2/8 = DF/(16 DT)
    // Also need to ensure that the sparse time domain transform is valid
    // => fdot > DF/T_w
    // fdot = gamma/tau
    
    fdot = 3.105*wdmpars->DF/wdmpars->Tw; // used an irrational fraction to ensure the fdot lands between samples
    tp = wdmpars->Tobs/2.0;
    fp = fdot*tp+4.0*wdmpars->DF;  // ensures that frequency starts positive
    printf("%f %f\n", fp/wdmpars->DF, tp/wdmpars->DT);
    gamma = fp/4.0;
    tau = gamma/fdot;
    Phi0 = 0.0;
    A0 = 10000.0;
    
    params = double_vector(10);
    params[0] = tau;     // time spread
    params[1] = 0.2;      // costh
    params[2] = 1.0;      // phi
    params[3] = A0;  // amplitude
    params[4] = -0.3;    // cosi
    params[5] = 0.8;   // psi
    params[6] = Phi0;  // phi0
    params[7] = gamma;  // frequency spread
    params[8] = tp;     // central time
    params[9] = fp;     // central frequency

    
    printf("computing time domain waveforms\n");
    h = double_vector(wdmpars->N);
    start = clock();
    ChirpTime(wdmpars, params, h);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Direct waveform calculation took %f seconds\n", cpu_time_used);
    out = fopen("chrp_time.dat","w");
    for(i=0; i< wdmpars->N; i++)
    {
        t = (double)(i)*dt;
        fprintf(out, "%e %.16e\n", t, h[i]);
    }
    fclose(out);
    free(h);
    printf("finished time domain waveforms\n");
    
    // Time domain Taylor expansion transform
    wavelet_TaylorT(wdmpars,params,"BinaryTaylorT.dat");
    
    // Time domain sparse transform
    wavelet_SparseT(wdmpars,params,"BinarySparseT.dat");
    
    // Frequency domain Taylor expansion transform
    wavelet_TaylorF(wdmpars,params,"BinaryTaylorF.dat");

    // Frequency domain sparse transform
    wavelet_SparseF(wdmpars,params,"BinarySparseF.dat");
    
    return 0;

}

void amp_phase_f(double f, double *params, double *phase, double *amp)
{
    double x;
    x = (f-params[9])/params[7];
    *phase = params[6] + 2.0*PI*params[8]*f + PI*params[0]*params[7]*x*x;
    *amp = params[3]*exp(-x*x/2.0);
}

void amp_phase_t(double t, double *params, double *phase, double *amp)
{
    double f, x, fac;
    x = (t-params[8])/params[0];
    f = params[7]*x+params[9];
    fac = sqrt(params[7]/params[0]);
    *phase = params[6] + 2.0*PI*(params[8]-t)*f + PI*params[7]*params[0]*x*x + PI/4.0;
    *amp = fac*params[3]*exp(-x*x/2.0);
}

double toff(double f, double *params)
{
    return(params[0]*(f-params[9])/params[7] + params[8]);
}

double foft(double t, double *params)
{
    return(params[7]*(t-params[8])/params[0]+params[9]);
}

void wavelet_SparseF(WDM *wdmpars, double *params, char *file)
{
    int i, j, k, l, m, mc;
    int Nsam, NM, NU;
    double **wave;
    int *Flist;
    double *waveF;
    char filename[100];
    double om, nrm, x, y, z;
    double *DX, *phif;
    
    FILE *out;
    
    clock_t start, end;
    double cpu_time_used;

    // Compute the frequency domain wavelet filter and normalize
    
      phif = double_vector(mult+1);
    
      for(j=0; j<= mult; j++)
      {
           om = wdmpars->dom*(double)(j);
           phif[j] = phitilde(om, wdmpars);
      }
      
      nrm = 0.0;
      for(l=-mult; l<= mult; l++) nrm += phif[abs(l)]*phif[abs(l)];
      nrm = sqrt(nrm/2.0);
      nrm *= 0.25*dt*(double)(Nt);
     
      for(j=0; j<= mult; j++) phif[j] /= nrm;
     
     NM = Nt*2*mult;
     
     Flist = int_vector(NM);
     waveF = double_vector(NM);
    
     sparse_wavelet_freq(wdmpars, params, phif, Flist, waveF, &NU);
     
     start = clock();
     for(mc=0; mc<1000; mc++)
     {
     sparse_wavelet_freq(wdmpars, params, phif, Flist, waveF, &NU);
     }
     end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("Sparse FD waveform calculation took %f seconds\n", cpu_time_used/1000.0);
    
     wave = double_matrix(Nt,Nf);  // wavelet wavepacket transform of the signal
    
    printf("Sparse frequency wavelet pixels used %d\n", NU);
     
     
     // initialize the array
     for(i=0; i< Nf; i++)
      {
       for(j=0; j< Nt; j++)
        {
            wave[j][i] = 0.0;
        }
     }
     
     // unpack the signal
     for(k=0; k< NU; k++)
     {
             j = (Flist[k]%Nt);
             i = (Flist[k]-j)/Nt;
             wave[j][i] = waveF[k];
     }

      out = fopen(file,"w");
        for(i=0; i< Nf; i++)
         {
          for(j=0; j< Nt; j++)
            {
            fprintf(out, "%e %e %.14e\n", (double)(j)*wdmpars->DT, (double)(i)*wdmpars->DF, wave[j][i]);
             // fprintf(out, "%d %d %.14e\n", j, i, wave[j][i]);
            }
          fprintf(out, "\n");
        }
       fclose(out);
    
    free_double_matrix(wave,Nt);
    free_int_vector(Flist);
    free_double_vector(waveF);
    free_double_vector(phif);
    
}


void wavelet_SparseT(WDM *wdmpars, double *params, char *file)
{
    int i, j, k, l, m, mc;
    int Nsam, NM, NU;
    double **wave;
    int *Tlist, *Nfsam;
    double *waveT;
    char filename[100];
    double om, nrm, x, y, z;
    double *DX, *phi;
    double **cH, **sH;
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *in;
    FILE *out;
    
    // Compute the time domain wavelet filter and normalize
    
     DX = (double*)malloc(sizeof(double)* (2*wdmpars->L));
     
        // negative frequencies
        for(i=1; i< wdmpars->L/2; i++)
        {
            om = -(double)(i)*wdmpars->dom;
            z = phitilde(om, wdmpars);
            REAL(DX,wdmpars->L-i) =  z;
            IMAG(DX,wdmpars->L-i) =  0.0;
        }
     
     //zero frequency
     REAL(DX,0) =  wdmpars->insDOM;
     IMAG(DX,0) =  0.0;
     
     // postive frequencies
     for(i=1; i<= wdmpars->L/2; i++)
     {
         om = (double)(i)*wdmpars->dom;
         z = phitilde(om, wdmpars);
         REAL(DX,i) =  z;
         IMAG(DX,i) =  0.0;
     }
    
     gsl_fft_complex_radix2_backward(DX, 1, wdmpars->L);
     
     phi = (double*)malloc(sizeof(double)* (wdmpars->L));
     
        for(i=0; i < wdmpars->L/2; i++)
         {
            phi[i] = REAL(DX,wdmpars->L/2+i);
         }
        for(i=0; i< wdmpars->L/2; i++)
         {
           phi[wdmpars->L/2+i] = REAL(DX,i);
        }
     
     free(DX);
     
         nrm = 0.0;
         for(l=0; l< wdmpars->L; l++) nrm += phi[l]*phi[l];
         nrm = sqrt(2.0*nrm/(double)(wdmpars->p));
         for(j=0; j< wdmpars->L; j++) phi[j] /= nrm;
    

     // pre-computed phase coefficients
     cH = double_matrix(Nf,wdmpars->L);
     sH = double_matrix(Nf,wdmpars->L);
     
      for(m=0; m< Nf; m++)
        {
         for(j=-wdmpars->L/2; j< wdmpars->L/2; j++)
         {
          i = j+wdmpars->L/2;
          x = TPI*(double)(j*m*mult)/(double)(wdmpars->L);
          cH[m][i] = cos(x);
          sH[m][i] = sin(x);
         }
        }
     
     NM = wdmpars->kx*Nt;
     
     // list of wavelet pixels that hold the signal
     Tlist = int_vector(NM);
     waveT = double_vector(NM);  // wavelet wavepacket transform of the signal

     sparse_wavelet_time(wdmpars, params, phi, Tlist, waveT, cH, sH, &NU);
     
     start = clock();
     for(mc=0; mc<1000; mc++)
     {
     sparse_wavelet_time(wdmpars, params, phi, Tlist, waveT, cH, sH, &NU);
     }
     end = clock();
     cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     printf("Sparse TD waveform calculation took %f seconds\n", cpu_time_used/1000.0);
    
    printf("Sparse time wavelet pixels used %d\n", NU);
    
     wave = double_matrix(Nt,Nf);  // wavelet wavepacket transform of the signal
     
      // initialize the array
      for(i=0; i< Nf; i++)
       {
        for(j=0; j< Nt; j++)
         {
             wave[j][i] = 0.0;
         }
      }
      
      // unpack the signal
      for(k=0; k< NU; k++)
      {
              j = (Tlist[k]%Nt);
              i = (Tlist[k]-j)/Nt;
              wave[j][i] = waveT[k];
      }
     
      out = fopen(file,"w");
        for(i=0; i< Nf; i++)
         {
          for(j=0; j< Nt; j++)
            {
            fprintf(out, "%e %e %.14e\n", (double)(j)*wdmpars->DT, (double)(i)*wdmpars->DF, wave[j][i]);
            }
          fprintf(out, "\n");
        }
       fclose(out);
    
    free_double_matrix(wave,Nt);
    free_int_vector(Tlist);
    free_double_vector(waveT);
    free_double_matrix(cH,Nf);
    free_double_matrix(sH,Nf);
    free_double_vector(phi);
}

void wavelet_TaylorT(WDM *wdmpars, double *params, char *file)
{
    int i, j, k, l, mc;
    int Nsam, NM, NU;
    double df, fdot;
    double *fd;
    double ***lookup;
    double **wave;
    int *Tlist, *Nfsam;
    double *waveT;
    char filename[100];
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *in;
    FILE *out;
    
    // frequency spacing
    df = (wdmpars->BW)/(double)(Nsf);
    
    fd = double_vector(Nfd);
    fd[0] = 0.0;
    fd[1] = wdmpars->DF/wdmpars->Tw*dfdot;  // sets the f-dot increment
    for(j=2; j< Nfd; j++) fd[j] = (double)(j)*fd[1];
    
    fdot = params[7]/params[0];
    
    printf("fdot %e fdot_max %e\n", fdot, fd[Nfd-1]);
    
    // number of samples for each frequency derivative layer (grows with increasing fdot)
    Nfsam = int_vector(Nfd);

    for(j=0; j< Nfd; j++)
    {
       Nfsam[j] = (int)((wdmpars->BW+fd[j]*wdmpars->Tw)/df);
       if(Nfsam[j]%2 != 0) Nfsam[j]++; // makes sure it is an even number
    }
    
    lookup = double_tensor(Nfd,Nfsam[Nfd-1],2);  // records the cos and sin coefficients
        
    for(j=0; j< Nfd; j++)
    {
        sprintf(filename, "coeffs/WDMcoeffs%d.dat", j);
        in = fopen(filename,"r");
        for(k=0; k< Nfsam[j]; k++)
        {
        fscanf(in,"%d%lf%lf\n", &i, &lookup[j][k][0], &lookup[j][k][1]);
        }
        fclose(in);
    }
    
    // max number of wavelet layers times the number of time pixels
    NM = (int)(ceil((wdmpars->BW+fd[Nfd-1]*wdmpars->Tw)/wdmpars->DF))*Nt;
    
    printf("max fast time samples %d\n", NM);
    
    // list of wavelet pixels that hold the signal
    Tlist = int_vector(NM);
    waveT = double_vector(NM);  // wavelet wavepacket transform of the signal
    
    TaylorTime(wdmpars, params, NM, df, Nfsam, fd, lookup, Tlist, waveT, &NU);

    printf("fast time samples used %d\n", NU);
    
    start = clock();
    for(mc=0; mc<1000; mc++)
    {
    TaylorTime(wdmpars, params, NM, df, Nfsam, fd, lookup, Tlist, waveT, &NU);
    }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Fast TD waveform calculation took %f seconds\n", cpu_time_used/1000.0);
    
    
    wave = double_matrix(Nt,Nf);  // wavelet wavepacket transform of the signal
    
    // initialize the array
    for(i=0; i< Nf; i++)
     {
      for(j=0; j< Nt; j++)
       {
           wave[j][i] = 0.0;
       }
    }
    
    // unpack the signal
    for(k=0; k< NU; k++)
    {
            j = (Tlist[k]%Nt);
            i = (Tlist[k]-j)/Nt;
            wave[j][i] = waveT[k];
    }
    
      out = fopen(file,"w");
        for(i=0; i< Nf; i++)
        {
           for(j=0; j< Nt; j++)
            {
            fprintf(out, "%e %e %.14e\n", (double)(j)*wdmpars->DT, (double)(i)*wdmpars->DF, wave[j][i]);
            }
        fprintf(out, "\n");
      }
     fclose(out);
    
    free_int_vector(Tlist);
    free_double_vector(waveT);
    free_double_matrix(wave,Nt);
    free_double_tensor(lookup,Nfd,Nfsam[Nfd-1]);
    free_int_vector(Nfsam);
     
}

// Fast frequency domain transform using lookup table
void wavelet_TaylorF(WDM *wdmpars, double *params, char *file)
{
    int i, j, k, l, mc;
    int Nsam, NMF, NU;
    double delt, fdot;
    double *td;
    double ***lookup;
    double **wave;
    int *Flist, *Nfsam;
    double *waveF;
    char filename[100];
    
    clock_t start, end;
    double cpu_time_used;
    
    FILE *in;
    FILE *out;

    // time spacing
    delt = (wdmpars->Tw/2.0)/(double)(Nst);
    
    td = double_vector(Ntd);
          
    td[0] = 0.0;
    td[1] = dtdf/(wdmpars->DF*wdmpars->DF);
    for(j=2; j< Ntd; j++) td[j] = (double)(j)*td[1];
    
    fdot = params[7]/params[0];
       
    printf("tprime %e tprime_max %e\n", 1.0/fdot, td[Ntd-1]);
          
    Nfsam = int_vector(Ntd);
          
          for(j=0; j< Ntd; j++)
          {
           Nfsam[j] = (int)((wdmpars->Tw/2.0+0.5*td[j]*wdmpars->DF)/delt);
          }
       
       lookup = double_tensor(Ntd,2*Nfsam[Ntd-1],2);  // records the coefficients
           
       for(j=0; j< Ntd; j++)
       {
           //printf("%d %d\n", j, Nfsam[j]);
           sprintf(filename, "coeffs/WDMcoeffsf%d.dat", j);
           in = fopen(filename,"r");
           for(k=0; k< 2*Nfsam[j]; k++)
           {
           fscanf(in,"%d%lf%lf\n", &i, &lookup[j][k][0], &lookup[j][k][1]);
           }
           fclose(in);
       }
       
       
       // list of wavelet pixels that hold the signal
          NMF = (int)((double)(2*Nfsam[Ntd-1])*delt/wdmpars->DT)*Nf;
    
          Flist = int_vector(NMF);
          waveF = double_vector(NMF);  // wavelet wavepacket transform of the signal
          printf("max fast freq samples %d\n", NMF);
       
          TaylorFreq(wdmpars, params, NMF, delt, Nfsam, td, lookup, Flist, waveF, &NU);
    
          printf("fast freq samples used %d\n", NU);

           
          start = clock();
          for(mc=0; mc<1000; mc++)
          {
           TaylorFreq(wdmpars, params, NMF, delt, Nfsam, td, lookup, Flist, waveF, &NU);
          }
          end = clock();
          cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
          printf("Fast FD waveform calculation took %f seconds\n", cpu_time_used/1000.0);
           
    
       wave = double_matrix(Nt,Nf);  // wavelet wavepacket transform of the signal
       
       // initialize the array
          for(i=0; i< Nf; i++)
           {
            for(j=0; j< Nt; j++)
             {
                 wave[j][i] = 0.0;
             }
          }
       
       
       // unpack the signal
          for(k=0; k< NU; k++)
          {
             
                  j = (Flist[k]%Nt);
                  i = (Flist[k]-j)/Nt;
                  wave[j][i] = waveF[k];
          }
       
       out = fopen(file,"w");
        for(i=0; i< Nf; i++)
            {
             for(j=0; j< Nt; j++)
                {
                fprintf(out, "%e %e %.14e\n", (double)(j)*wdmpars->DT, (double)(i)*wdmpars->DF, wave[j][i]);
                }
              fprintf(out, "\n");
            }
        fclose(out);
       
       free_double_matrix(wave,Nt);
       free_double_tensor(lookup,Ntd,2*Nfsam[Ntd-1]);
       free_double_vector(td);
       free(waveF);
       free(Flist);
       free_int_vector(Nfsam);
    
    
}

void sparse_wavelet_time(WDM *wdmpars, double *params, double *phi, int *list, double *wave, double **cH, double **sH, int *NU)
{
    double BW, Tw, Tobs, dom, om, nrm;
    double f, t;
    int i, j, k, l, p, is;
    int m, mc, mq, n, dn, mm, nn, nu;
    int K, kx;
    double tau, A0, Phi0, gamma, tp, fp, x, y, dtds;
    double s1, s2, fg, fac;
    double *Amp, *Phase, *tlayer, *DX;
    double *cP, *sP;
    double Ax, px;
    
    FILE *out;
    
    // indicates this pixel not used
       for(i=0; i< (wdmpars->kx*Nt); i++)
       {
           list[i] = -1;
           wave[i] = 0.0;
       }
    
    //nn = Nf/p;
    nn = (Nf*wdmpars->L)/wdmpars->K;

    dtds = dt*(double)(wdmpars->p);  // downsampled data spacing in time

    Amp = double_vector(wdmpars->M);
    //Phase = double_vector(wdmpars->M);
    cP = double_vector(wdmpars->M);
    sP = double_vector(wdmpars->M);
    
    
    for(j=0; j< wdmpars->M; j++)
    {
        t = dtds*(double)(j);
        amp_phase_t(t, params, &x, &Amp[j]);
        //Phase[j] = -x;
        cP[j] = cos(x);
        sP[j] = -sin(x);
    }
    
    DX = double_vector(2*wdmpars->L);
    
    nu = 0;
    
    for(n=0; n< Nt; n++)
    {
        t = (double)(n)*wdmpars->DT;
        is = nn*n;
        f = foft(t,params);
        mc = (int)(f/wdmpars->DF);
        //printf("central m layer %d\n", mc);
        mq = mc*mult;
        
        if(mc > -1 && mc < Nf)
        {
        
        // have sped this up by computing cos(Phase[k]), sin(Phase[k])
        // and passing in pre-computed arrays for cos(TPI*(double)(j*mq)/(double)(L))
        // and sin(TPI*(double)(j*mq)/(double)(L)). The array is L*Nf in size
        for(j=-wdmpars->L/2; j< wdmpars->L/2; j++)
        {
            i = j+wdmpars->L/2;
            k = j+is;
            //printf("%d %d %d\n", n, k, M);
            if(k > -1 && k < wdmpars->M)
            {
            Ax =  Amp[k]*phi[i];
            //px =  Phase[k]-TPI*(double)(j*mq)/(double)(L);
           // REAL(DX,i) = Ax*cos(px);
            //IMAG(DX,i) = Ax*sin(px);
             REAL(DX,i) = Ax*(cP[k]*cH[mc][i]+sP[k]*sH[mc][i]);
             IMAG(DX,i) = Ax*(sP[k]*cH[mc][i]-cP[k]*sH[mc][i]);
            }
            else
            {
                REAL(DX,i) = 0.0;
                IMAG(DX,i) = 0.0;
            }
            
        }
        
        gsl_fft_complex_radix2_forward(DX, 1, wdmpars->L);
       
        
        //negative frequencies
        for(j=wdmpars->kx/2; j< wdmpars->kx; j++)
        {
          m = j-wdmpars->kx+mc;
          if(m > -1 && m < Nf)
          {
            if((n+m)%2 ==0)
            {
              wave[nu] = REAL(DX,j*mult);
            }
            else
            {
              wave[nu] = -IMAG(DX,j*mult);
            }
              list[nu] = n+m*Nt;
              nu++;
          }
        }
        
        //postive frequencies
        for(j=0; j< wdmpars->kx/2; j++)
        {
            m = j+mc;
            if(m > -1 && m < Nf)
            {
              if((n+m)%2 ==0)
              {
                wave[nu] = REAL(DX,j*mult);
              }
              else
              {
                wave[nu] = -IMAG(DX,j*mult);
              }
                list[nu] = n+m*Nt;
                nu++;
            }
        }
            
        }
    
    }
    
    free(DX);
    free(Amp);
    //free(Phase);
    free(cP);
    free(sP);
    
    *NU = nu;
    
}

void sparse_wavelet_freq(WDM *wdmpars, double *params, double *phif, int *list, double *wave, int *NU)
{
    double f, t;
    int i, j, k, l, is;
    int m, mq, n, dn, mm, nn;
    double x, y;
    double A, P;
    double s1, s2, sx;
    double *Amp, *Phase, *flayer, *DX;
    
    Amp = double_vector(wdmpars->K);
    Phase = double_vector(wdmpars->K);

    for(j=0; j< wdmpars->K; j++)
    {
        f= (double)(j)/wdmpars->Tw;
        amp_phase_f(f, params, &P, &A);
        Amp[j] = A;
        Phase[j] = -P;
    }
    
    flayer = double_vector(2*mult);
    DX = double_vector(4*mult);
    
    mm = 0;
    
     for(m=0; m< Nf; m++)
      {
          
        mq = m*mult;
        f = (double)(m)*wdmpars->DF;
        t = toff(f,params);
        n = (int)(t/wdmpars->DT);
          
          sx = -1.0;
          if(n%2 == 0 || m%2 == 0)
          {
          sx = 1.0;
          }
          
          
        // These are used to unpack the iFFT.
        if(m%2 == 0)
        {
           if(n%2 == 0)
           {
               s1 = 1.0;
               s2 = 1.0;
           }
           else
           {
               s1 = -1.0;
               s2 = -1.0;
           }
        }
        else
        {
                s1 = 1.0;
                s2 = -1.0;
        }
        
          
         if(t > 0.0 && t < wdmpars->Tobs)
         {
             //printf("%e %e\n", f, t);
             // Using a lookup table here is actually slower!
             // The reason is that both require K=2*mult*Nf sines and cosines
             // But the lookup also requires additional multiplies and adds
             for(l=-mult; l< mult; l++)
             {
                 REAL(DX,l+mult) = 0.0;
                 IMAG(DX,l+mult) = 0.0;
                 i = mq+l;
                 if(i >= 0 && i < wdmpars->K)
                 {
                   y = Phase[i]+ TPI*(double)(n*Nf)*((double)(i))/(double)(wdmpars->K);
                   REAL(DX,l+mult) = phif[abs(l)]*Amp[i]*cos(y);
                   IMAG(DX,l+mult) = phif[abs(l)]*Amp[i]*sin(y);
                 }
             }
              
            gsl_fft_complex_radix2_backward(DX, 1, (2*mult));

             
        //Unpacking the complex Fourier transform. This took some trial and
        //error to get right!
             for(j=0; j < mult; j++)
             {
                 if((j+n+m)%2 == 0)
                 {
                     x = s1*REAL(DX,j+mult);
                     y = s1*REAL(DX,j);
                 }
                 else
                 {
                    x = s2*IMAG(DX,j+mult);
                    y = s2*IMAG(DX,j);
                 }
                 
                 flayer[j] = x;
                 flayer[j+mult] = y;
                }
          
           for(j=-mult; j< mult; j++)
           {
               nn = n+j;
               if(nn > 0 && nn < Nt)
               {
               // which wavelet is lit up
               list[mm] = nn+m*Nt;
               // value of said wavelet
               wave[mm] = flayer[j+mult];
               mm++;
               }
            }
         }
          
       }
    
    // number of pixels used
    *NU = mm;
    
    free_double_vector(DX);
    free_double_vector(flayer);
    free_double_vector(Amp);
    free_double_vector(Phase);
    
}



double phitilde(double om, WDM *wdmpars)
{
    double x, y, z;
    
       z = 0.0;
       
       if(fabs(om) >= wdmpars->A && fabs(om) < wdmpars->A+wdmpars->B)
       {
           x = (fabs(om)-wdmpars->A)/wdmpars->B;
           y = gsl_sf_beta_inc(nx, nx, x);
           z = wdmpars->insDOM*cos(y*PI/2.0);
       }
       
       if(fabs(om) < wdmpars->A) z = wdmpars->insDOM;
    
    return(z);
    
}



void TaylorFreq(WDM *wdmpars, double *params, int NMF, double delt, int *Nfsam, double *td, double ***lookup, int *Flist, double *waveF, int *NU)
{
    
    int i, j, k;
    double *Amp, *Phase;
    double *ta, *tpa;
    double *FT;

    // indicates this pixel not used
    for(i=0; i< NMF; i++)
    {
        Flist[i] = -1;
        waveF[i] = 0.0;
    }
    
    Amp = (double*)malloc(sizeof(double)* (Nf));
    Phase = (double*)malloc(sizeof(double)* (Nf));
    ta = (double*)malloc(sizeof(double)* (Nf));  // t(f)
    tpa = (double*)malloc(sizeof(double)* (Nf)); // t'(f)
    FT = (double*)malloc(sizeof(double)* (Nf));
    
    for (i=0; i< Nf; i++)
    {
        FT[i] = ((double)(i))*wdmpars->DF;  // center of the pixel
    }
    
    // Note: with the simple analytic model we don't have to evaluate the
    // times and time derivatives numerically. Only doing so
    // here to illustrate the general method that can be used with any
    // waveform model
    ChirpWaveletF(wdmpars, params, FT, Amp, Phase, ta, tpa);
    
    wavemakef(wdmpars, delt, Nfsam, td, Phase, ta, tpa, Amp, lookup, Flist, waveF, NU);
    
    free(Amp);
    free(Phase);
    free(ta);
    free(tpa);
    free(FT);

    
}


void TaylorTime(WDM *wdmpars, double *params, int NM, double df, int *Nfsam, double *fd, double ***lookup, int *Tlist, double *waveT, int *NU)
{
    
    int i, j, k;
    double *Amp, *Phase;
    double *fa, *fda;
    double *TF;

    // indicates this pixel not used
    for(i=0; i< NM; i++)
    {
        Tlist[i] = -1;
        waveT[i] = 0.0;
    }
    
    Amp = (double*)malloc(sizeof(double)* (Nt));
    Phase = (double*)malloc(sizeof(double)* (Nt));
    fa = (double*)malloc(sizeof(double)* (Nt));  // f(t)
    fda = (double*)malloc(sizeof(double)* (Nt)); // fdot(t)
    TF = (double*)malloc(sizeof(double)* (Nt));
    
    for (i=0; i< Nt; i++)
    {
        TF[i] = ((double)(i))*wdmpars->DT;  // center of the pixel
    }
    
    // Note: with the simple analytic model we don't have to evaluate the
    // frequencies and frequency derivatives numerically. Only doing so
    // here to illustrate the general method that can be used with any
    // waveform model
    ChirpWaveletT(wdmpars, params, TF, Amp, Phase, fa, fda);
    
    wavemaket(wdmpars, df, Nfsam, fd, Phase, fa, fda, Amp, lookup, Tlist, waveT, NU);
    
    free(Amp);
    free(Phase);
    free(fa);
    free(fda);
    free(TF);

    
}

void wavemakef(WDM *wdmpars, double delt, int *Nfsam, double *td, double *Phase, double *ta, double *tpa, double *Amp, double ***lookup, int *list, double *wave, int *NU)
{
    
    int i, j, k, kmin, kmax, is, ii, jj, kk, kx, n, dn, mm, mx;
    int k1, k2;
    double dx1, dx2, dy;
    double phase, f, t, fmid, fsam, A;
    double c, s, x, y, z, yy, zz, tprime;
    double tdmx, tdt, t0, xt, tx, t1, t2, dtd;
    int flag;
    
    // maximum second phase derivative
    tdmx = td[Ntd-1];
    dtd = td[1]; //  increment
    
    mm = 0;

        //is = (int)(fp/DF)+1;
        //for(j=is; j< is+1; j++)
        for(j=0; j< Nf; j++)
        {
            
           // central frequency of layer
            f = (double)(j)*wdmpars->DF;
            
            A = Amp[j];
            phase = Phase[j];
            t = ta[j];
            n = (int)(t/wdmpars->DT);
            
            c = A*cos(phase);
            s = A*sin(phase);
            
            tprime = tpa[j];
            // lower t-prime layer
            ii = (int)(floor(tprime/dtd));
            // linearly interpolate over t-prime
            dy = tprime/dtd-(double)(ii);
            
            if(ii < Ntd-2) // only proceed if t-prime is covered by lookup table
            {
            
            // time pixel range at lower t-prime
            jj = (int)((double)(Nfsam[ii])*delt/wdmpars->DT);
            
            // printf("%f %f %f %f\n", A, phase, c, s);
            
            // central time of pixel nearest to track
            t0 = (double)(n)*wdmpars->DT;
           
            // offset from central time
            xt = t-t0;
            
           // printf("%d %f %e\n", n, xt, xt/delt);
            
             // other pixels will be offset by (double)(k-n)*DT+xt
                   
          if(t > 0.0 && t < wdmpars->Tobs)
          {
              
    
           // k = 0 is central pixel in time
           for(k=-jj; k< jj; k++)
           {
               
            //kx = k+n;
            
            // time direction is reversed. Maybe be due to sign issue in phase
            kx = n - k;
               
           // printf("%d\n", kx);
               
            if(kx > -1 && kx < Nt)
            {
            
            tx = (double)(k)*wdmpars->DT+xt;
            
            // time interpolation for lower beta layer
            t1 = (double)(Nfsam[ii])*delt + tx;
            k1 = (int)(floor(t1/delt));
            dx1 = t1/delt-(double)(k1);
                
            // time interpolation for upper beta layer
            t2 = (double)(Nfsam[ii+1])*delt + tx;
            k2 = (int)(floor(t2/delt));
            dx2 = t2/delt-(double)(k2);
                
           // printf("%d %d %d %f %f\n", 2*Nfsam[ii], k1, k2, dx1, dx2);
            
            if(k1 > -1 && k2 > -1 && k1 < 2*Nfsam[ii]-1 && k2 < 2*Nfsam[ii+1]-1)
            {
            // interpolate over time
           
            y = (1.0-dx1)*lookup[ii][k1][0]+dx1*lookup[ii][k1+1][0];
            z = (1.0-dx1)*lookup[ii][k1][1]+dx1*lookup[ii][k1+1][1];
            
            yy = (1.0-dx2)*lookup[ii+1][k2][0]+dx2*lookup[ii+1][k2+1][0];
            zz = (1.0-dx2)*lookup[ii+1][k2][1]+dx2*lookup[ii+1][k2+1][1];
                
            // interpolate over beta
            y = (1.0-dy)*y+dy*yy;
            z = (1.0-dy)*z+dy*zz;
                
               // printf("y z %e %e\n", y, z);
                
            if(j%2 == 0)
            {
             if((j+kx)%2 == 0)
             {
                wave[mm] = (c*y-s*z);
             }
             else
             {
                wave[mm] = (c*z+s*y);
             }
            }
            else
            {
              if((j+kx)%2 == 0)
                {
                 wave[mm] = -(c*y-s*z);
                }
                else
                {
                 wave[mm] = (c*z+s*y);
                }
            }
                
                
            list[mm] = kx+j*Nt;
            mm++;
                
            }
                
            }
            
           }
        
          }
              
          }
            
        }
    
    *NU = mm;
    
   // printf("%d\n", mm);
    
}



void wavemaket(WDM *wdmpars, double df, int *Nfsam, double *fd, double *Phase, double *freq, double *freqd, double *Amp, double ***lookup, int *list, double *wave, int *NU)
{
    
    int i, j, k, kmin, kmax, ii, jj, kk, n, mm, mx;
    double dx, dy;
    double phase, f, fdot, fmid, fsam, A;
    double c, s, x, y, z, yy, zz;
    double fmx, fdmx, dfd, HBW;
    int flag;
    
    // maximum frequency and frequency derivative
    fmx = (double)(Nf-1)*wdmpars->DF;
    fdmx = fd[Nfd-1];
    dfd = fd[1]; // f-dot increment
    
    mm = 0;

        for(j=0; j< Nt; j++)
        {
    
           phase = Phase[j];
           f = freq[j];
           fdot = freqd[j];
           A = Amp[j];
            
           if(f < fmx && fdot < fdmx)
           {
            
           // lower f-dot layer
           n = (int)(floor(fdot/dfd));
               
           dy = fdot/dfd-(double)(n);
               
            HBW  = 0.5*(double)(Nfsam[n]-1)*df;
           
           // lowest frequency layer
           kmin = (int)(ceil((f-HBW)/wdmpars->DF));
               
            // highest frequency layer
           kmax = (int)(floor((f+HBW)/wdmpars->DF));
               
               
              // printf("%d %d\n", kmin, kmax);
            
         
           c = A*cos(phase);
           s = A*sin(phase);
               
            for(k=kmin; k<= kmax; k++)
            {
            
            // central frequency
            fmid = (double)(k)*wdmpars->DF;
                
            kk = (int)(floor((f-(fmid+0.5*df))/df));
            fsam = fmid+((double)(kk)+0.5)*df;
            dx = (f-fsam)/df; // used for linear interpolation
            //printf("%d %d %d %d %d %e %f %f\n", j, k, kk, kk+Nfsam[n]/2, Nfsam[n], fsam, dx, dy);
                
            // interpolate over frequency
            jj = kk+Nfsam[n]/2;
            y = (1.0-dx)*lookup[n][jj][0]+dx*lookup[n][jj+1][0];
            z = (1.0-dx)*lookup[n][jj][1]+dx*lookup[n][jj+1][1];
            jj = kk+Nfsam[n+1]/2;
            yy = (1.0-dx)*lookup[n+1][jj][0]+dx*lookup[n+1][jj+1][0];
            zz = (1.0-dx)*lookup[n+1][jj][1]+dx*lookup[n+1][jj+1][1];
                
            // interpolate over fdot
            y = (1.0-dy)*y+dy*yy;
            z = (1.0-dy)*z+dy*zz;
                
                
            if((j+k)%2 == 0)
            {
                wave[mm] = (c*y-s*z);
            }
            else
            {
                wave[mm] = -(c*z+s*y);
            }
            list[mm] = j+k*Nt;
            mm++;
                
            }  // end loop over frequency layers
            
        }
            
        }
    
    *NU = mm;
    
   // printf("%d\n", mm);
    
}


void ChirpTime(WDM *wdmpars, double *params, double *h)
{
    double *Amp, *Phase;
    double af, fr, df, fac, deltaF, t, f, x;
    int i, j, n;
    double AA, P;
    double px, fnew, fonfs;
    double f0, fdot0, fddot0;
    double fA, fE;
    double fdA, fdE;
    
    int NF = (int)(wdmpars->Tobs/wdmpars->DT)+1;
    
    double *TF;

    FILE *out;
    
    TF = (double*)malloc(sizeof(double)* (NF));
    
    for (i=0; i< NF; i++)
    {
        t = (double)(i)*wdmpars->DT;
        TF[i] = t;
    }
    
    Amp = (double*)malloc(sizeof(double)* (NF));
    Phase = (double*)malloc(sizeof(double)* (NF));
    

    //  compute the phase and amplitude
    for(n=0; n< NF; n++) amp_phase_t(TF[n], params, &Phase[n], &Amp[n]);
       
    // flip phase
    for(n=0; n< NF; n++) Phase[n] = -Phase[n];
      
    gsl_interp_accel *Pacc = gsl_interp_accel_alloc();
    gsl_spline *Pspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Pspline, TF, Phase, NF);
    
    gsl_interp_accel *Aacc = gsl_interp_accel_alloc();
    gsl_spline *Aspline = gsl_spline_alloc (gsl_interp_cspline, NF);
    gsl_spline_init(Aspline, TF, Amp, NF);
    
    out = fopen("fastAP_time.dat","w");
    
    for (i=0; i< wdmpars->N; i++)
    {
        t = dt*(double)(i);
        
        // The instaneous frequency is given by 1/(2 pi) dPhase/dt
        // Includes the doppler shift and the polarization phase
        
        
            P = gsl_spline_eval (Pspline, t, Pacc);
            AA = gsl_spline_eval (Aspline, t, Aacc);
            h[i] = AA*cos(P);
        
            // instaneous frequency of the A channel signal
            fA = gsl_spline_eval_deriv(Pspline, t, Pacc)/TPI;
            
            // fdot
            fdA = gsl_spline_eval_deriv2(Pspline, t, Pacc)/TPI;
        
        fprintf(out,"%e %e %e %e\n", t, AA, fA, fdA);
 
    }
    
    fclose(out);
    
    free(TF);
    free(Amp);
    free(Phase);
  
    gsl_spline_free(Pspline);
    gsl_spline_free(Aspline);
    gsl_interp_accel_free(Pacc);
    gsl_interp_accel_free(Aacc);

}

void ChirpWaveletF(WDM *wdmpars, double *params, double *FT, double *Amp, double *Phase, double *ta, double *tpa)
{
    double af, fr, df, fac, deltaF, t, f, x;
    int i, j, n;
    int Nfs;
    double A, P, fny;
    double dfs;
    double *AS, *PS, *FS;
    
    fny = (1.0/(2.0*dt));
    
    // sample spacing for spline
    dfs = fny/200.0;
    
    // pad out past ends to avoid edge effects with the spline
    Nfs = (int)(fny/dfs)+11;
    
    FS = double_vector(Nfs);
    AS = double_vector(Nfs);
    PS = double_vector(Nfs);
    
    for (i=0; i< Nfs; i++) FS[i] = (double)(i-5)*dfs;

    //  compute the phase and amplitude
    for(n=0; n< Nfs; n++) amp_phase_f(FS[n], params, &PS[n], &AS[n]);

    gsl_interp_accel *Aacc = gsl_interp_accel_alloc();
    gsl_spline *Aspline = gsl_spline_alloc (gsl_interp_cspline, Nfs);
    gsl_spline_init(Aspline, FS, AS, Nfs);
      
    gsl_interp_accel *Pacc = gsl_interp_accel_alloc();
    gsl_spline *Pspline = gsl_spline_alloc (gsl_interp_cspline, Nfs);
    gsl_spline_init(Pspline, FS, PS, Nfs);
    
    free_double_vector(FS);
    free_double_vector(AS);
    free_double_vector(PS);
    
    for (i=0; i< Nf; i++)
    {
        f = FT[i];
        
        Amp[i] = gsl_spline_eval(Aspline, f, Aacc);
        // instaneous phase, time and dt/df of the signal
        Phase[i] = gsl_spline_eval(Pspline, f, Pacc);
        ta[i] = gsl_spline_eval_deriv(Pspline, f, Pacc)/TPI;
        tpa[i] = gsl_spline_eval_deriv2(Pspline, f, Pacc)/TPI;

    }
    
    gsl_spline_free(Aspline);
    gsl_interp_accel_free(Aacc);
    gsl_spline_free(Pspline);
    gsl_interp_accel_free(Pacc);
    
}


void ChirpWaveletT(WDM *wdmpars, double *params, double *TF, double *Amp, double *Phase, double *fa, double *fda)
{
    double af, fr, df, fac, deltaF, t, f, x;
    int i, j, n;
    int Nts;
    double A, P;
    double px, fnew, fonfs;
    double f0, fdot0, fddot0;
    double fA, fE, dts;
    double *TS, *AS, *PS;
    
    FILE *out;
    
    // sample spacing for spline
    dts = wdmpars->Tobs/300.0;
    
    // pad out past ends to avoid edge effects with the spline
    Nts = (int)(wdmpars->Tobs/dts)+11;
    
    TS = double_vector(Nts);
    AS = double_vector(Nts);
    PS = double_vector(Nts);
    
    for (i=0; i< Nts; i++) TS[i] = (double)(i-5)*dts;
    
       //  compute the phase and amplitude
    for(n=0; n< Nts; n++) amp_phase_t(TS[n], params, &PS[n], &AS[n]);
     
    // need sign flip on phase
    for(n=0; n< Nts; n++) PS[n] = -PS[n];

    gsl_interp_accel *Aacc = gsl_interp_accel_alloc();
    gsl_spline *Aspline = gsl_spline_alloc (gsl_interp_cspline, Nts);
    gsl_spline_init(Aspline, TS, AS, Nts);
      
    gsl_interp_accel *Pacc = gsl_interp_accel_alloc();
    gsl_spline *Pspline = gsl_spline_alloc (gsl_interp_cspline, Nts);
    gsl_spline_init(Pspline, TS, PS, Nts);
    
    free_double_vector(TS);
    free_double_vector(AS);
    free_double_vector(PS);
    
    for (i=0; i< Nt; i++)
    {
        t = TF[i];
        
        Amp[i] = gsl_spline_eval(Aspline, t, Aacc);
        // instaneous phase, frequency, fdot of the signal
        Phase[i] = gsl_spline_eval(Pspline, t, Pacc);
        fa[i] = gsl_spline_eval_deriv(Pspline, t, Pacc)/TPI;
        fda[i] = gsl_spline_eval_deriv2(Pspline, t, Pacc)/TPI;

    }
    
    gsl_spline_free(Aspline);
    gsl_interp_accel_free(Aacc);
    gsl_spline_free(Pspline);
    gsl_interp_accel_free(Pacc);
    
}





double *double_vector(int N)
{
    return malloc( (N+1) * sizeof(double) );
}

void free_double_vector(double *v)
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

double ***double_tensor(int N, int M, int L)
{
    int i,j;
    
    double ***t = malloc( (N+1) * sizeof(double **));
    for(i=0; i<N+1; i++)
    {
        t[i] = malloc( (M+1) * sizeof(double *));
        for(j=0; j<M+1; j++)
        {
            t[i][j] = malloc( (L+1) * sizeof(double));
        }
    }
    
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    int i;
    
    for(i=0; i<N+1; i++) free_double_matrix(t[i],M);
    
    free(t);
}


int *int_vector(int N)
{
    return malloc( (N+1) * sizeof(int) );
}

void free_int_vector(int *v)
{
    free(v);
}

void create_WDM(WDM *wdmpars)
{
    wdmpars->N = Nt*Nf;   // total points
    wdmpars->Tobs = dt*(double)(wdmpars->N);   // duration
    wdmpars->DT = dt*(double)(Nf);           // width of wavelet pixel in time
    wdmpars->DF = 1.0/(2.0*dt*(double)(Nf));   // width of wavelet pixel in frequency
    wdmpars->K = mult*2*Nf; // filter length
    wdmpars->Tw = dt*(double)(wdmpars->K); // filter duration
    wdmpars->L = 256;  // reduced filter length - must be a power of 2
    wdmpars->p = wdmpars->K/wdmpars->L; // downsample factor
    wdmpars->kx = 2*Nf/wdmpars->p;
    wdmpars->dom = TPI/wdmpars->Tw; // angular frequency spacing
    wdmpars->OM = PI/dt;  // Nyquist angular frequency
    wdmpars->DOM = wdmpars->OM/(double)(Nf); // 2 pi times DF
    wdmpars->insDOM = 1.0/sqrt(wdmpars->DOM);
    wdmpars->B = wdmpars->OM/(double)(2*Nf);
    wdmpars->A = (wdmpars->DOM-wdmpars->B)/2.0;
    wdmpars->BW = (wdmpars->A+wdmpars->B)/PI; // total width of wavelet in frequency
    // nonzero terms in phi transform (only need 0 and positive)
    wdmpars->np = (int)(wdmpars->BW*wdmpars->Tw/2);
    wdmpars->M = Nf*Nt/wdmpars->p;  // down-sampled total points
}


