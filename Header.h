typedef struct
{
  int N;
  int M;
  int K;
  int L;
  int p;
  int kx;
  int np;
  double Tobs;
  double DT;
  double DF;
  double Tw;
  double dom;
  double OM;
  double DOM;
  double insDOM;
  double A;
  double B;
  double BW;
    
}WDM;

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

int *int_vector(int N);
void free_int_vector(int *v);
double *double_vector(int N);
void free_double_vector(double *v);
double **double_matrix(int N, int M);
void free_double_matrix(double **m, int N);
double ***double_tensor(int N, int M, int L);
void free_double_tensor(double ***t, int N, int M);

double foft(double t, double *params);
double toff(double f, double *params);
void amp_phase_f(double f, double *params, double *phase, double *amp);
void amp_phase_t(double t, double *params, double *phase, double *amp);
double phitilde(double om, WDM *wdmpars);
void ChirpTime(WDM *wdmpars, double *params, double *h);
void wavelet_TaylorT(WDM *wdmpars, double *params, char *file);
void wavelet_TaylorF(WDM *wdmpars, double *params, char *file);
void wavelet_SparseT(WDM *wdmpars, double *params, char *file);
void wavelet_SparseF(WDM *wdmpars, double *params, char *file);
void TaylorTime(WDM *wdmpars, double *params, int NM, double df, int *Nfsam, double *fd, double ***lookup, int *Tlist, double *waveT, int *NU);
void TaylorFreq(WDM *wdmpars, double *params, int NMF, double delt, int *Nfsam, double *td, double ***lookup, int *Flist, double *waveF, int *NU);
void wavemaket(WDM *wdmpars, double df, int *Nfsam, double *fd, double *Phase, double *freq, double *freqd, double *Amp, double ***lookup, int *list, double *wave, int *NU);
void sparse_wavelet_freq(WDM *wdmpars, double *params, double *phif, int *list, double *wave, int *NU);
void sparse_wavelet_time(WDM *wdmpars, double *params, double *phi, int *list, double *wave, double **cH, double **sH, int *NU);
void wavemakef(WDM *wdmpars, double delt, int *Nfsam, double *td, double *Phase, double *ta, double *tpa, double *Amp, double ***lookup, int *list, double *wave, int *NU);
void ChirpWaveletT(WDM *wdmpars, double *params, double *TF, double *Amp, double *Phase, double *fa, double *fda);
void ChirpWaveletF(WDM *wdmpars, double *params, double *FT, double *Amp, double *Phase, double *ta, double *tpa);
void create_WDM(WDM *wdmpars);

