/*
	StencilProbe Heat Equation
	Implements 7pt stencil from Chombo's heattut example.
*/

#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SIZE
#define SIZE 500
#endif

volatile int log_ready = 0;
double log_data[SIZE*SIZE*SIZE];
pthread_t *log_id = 0;

void log_all_data(unsigned *size)
{
  static char logname[10];
  static int logindex = 0;

  for (;;) {
    while (!log_ready); 
    logindex = (logindex + 1)%2;
    sprintf(logname, "log_%d",logindex);
    FILE* logfile=fopen(logname, "w");
    fwrite(&log_ready, sizeof(log_ready), 1, logfile);
    fwrite(log_data, sizeof(double), *size, logfile);
    fclose(logfile);
    log_ready = 0;
  }
}

int recover_data(double* data, unsigned size)
{
  FILE* logfile0=fopen("log_0", "r");
  FILE* logfile1=fopen("log_1", "r");
  if (logfile0==0 || logfile1 == 0) return 0;
  int index0, index1, res;
  fread(&index0, sizeof(int), 1, logfile0);
  fread(&index1, sizeof(int), 1, logfile1);
  if (index0 < index1) {
     fread(data, sizeof(double), size, logfile0);
     res = index0;
  } 
  else {
     fread(data, sizeof(double), size, logfile1);
     res = index1;
  } 
  fclose(logfile0);
  fclose(logfile1);
  return res;
}

void jacobi7_3(const int nx,const int ny, int nz, const double alpha,double* A0,const int timesteps,const double* B,const int ldb, double* Anext,const int ldc) 
{
  double fac;
  double *temp_ptr;
  int i, j, k, t, tstart=0;	
  fac = 6.0/(A0[0]*A0[0]);
  double *l0, *lnext;

  unsigned size = nx*ny*nz;

  if (log_id == 0) {
    tstart=recover_data(A0, size);
    if (tstart > 0) 
       memcpy((void*)Anext,(void*)A0,size * sizeof(double));
    log_id = (pthread_t*) malloc(sizeof(pthread_t));
    pthread_create(log_id,NULL, log_all_data, &size); 
  }

  /*@;BEGIN(Nest1=Nest)@*/for (t = tstart; t < timesteps; t++) {
          int do_log = !(log_ready);
          if (t%2 == 0) { l0 = A0; lnext = Anext; }
          else {lnext = A0; l0 = Anext; }
    /*@;BEGIN(Nest2=Nest)@*/for (k = 1; k < nz - 1; k++) {
      /*@;BEGIN(Nest3=Nest)@*/for (j = 1; j < ny - 1; j++) {
	/*@;BEGIN(Nest4=Nest)@*/for (i = 1; i < nx - 1; i++) {
	  lnext[Index3D (nx, ny, i, j, k)] = 
	    l0[Index3D (nx, ny, i, j, k + 1)] +
	    l0[Index3D (nx, ny, i, j, k - 1)] +
	    l0[Index3D (nx, ny, i, j + 1, k)] +
	    l0[Index3D (nx, ny, i, j - 1, k)] +
	    l0[Index3D (nx, ny, i + 1, j, k)] +
	    l0[Index3D (nx, ny, i - 1, j, k)]
	    -  l0[Index3D (nx, ny, i, j, k)] *fac ;
	}
      }
    }
    if (do_log) {
       memcpy((void*)log_data,(void*)lnext,size * sizeof(double));
       log_ready=t+1;
    }
  }
}

