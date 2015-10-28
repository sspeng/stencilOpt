/*
	StencilProbe Heat Equation
	Implements 7pt stencil from Chombo's heattut example.
*/

#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#ifndef SIZE
#define SIZE 500
#endif

volatile int log_ready = 0;
int log_index1,log_index2;
int pt[4],logtime=0;
double log_data[SIZE*SIZE*SIZE],pmem[4][SIZE*SIZE*SIZE+1];
int st[4];
pthread_t *log_id = 0;

static inline void asm_clflush(volatile intptr_t *addr)
{
    __asm__ __volatile__ ("clflush %0 "::"m"(*addr));
}
static inline void asm_mfence(void)
{
    __asm__ __volatile__ ("mfence");
}
static inline uint64_t  rdtsc()
{
    unsigned long a,d;
    __asm__ __volatile__ ("cpuid; rdtsc":"=a"(a),"=d"(d)::"ebx","ecx");
    return a | ((uint64_t)d<<32);
}
double GetWallTime1(void)
{
    struct timeval tp1;
    static long start1=0, startu1;
    if (!start1)
        {
            gettimeofday(&tp1, NULL);
            start1 = tp1.tv_sec;
            startu1 = tp1.tv_usec;
            //printf("1-----%ld\n",start1);
            return(0.0);
        }
    gettimeofday(&tp1, NULL);
    //printf("2-----%ld\n",tp1.tv_sec);
    return( ((double) (tp1.tv_sec - start1)) + (tp1.tv_usec-startu1)/1000000.0 );
}

void flush(int index,unsigned size)
{
    int i,j;
    for (i=0;i<=size;i++)
    asm_clflush(&pmem[index][i]);
    asm_mfence();
    for (i=0;i<=3;i++)
    {
        asm_clflush(&st[i]);
        asm_clflush(&pt[i]);
    }
    asm_mfence();
}

void log_all_data(unsigned *size)
{
  static char logname[10];
  int t;
  //static int logindex = 0;

  for (;;) {
    while (!log_ready); 
    //logindex = (logindex + 1)%2;
    //sprintf(logname, "log_%d",logindex);
    //FILE* logfile=fopen(logname, "w");
    //fwrite(&log_ready, sizeof(log_ready), 1, logfile);
    //fwrite(log_data, sizeof(double), *size, logfile);
    //fclose(logfile);
    //memcpy((void*)pmem1,(void*)log_data,(*size) * sizeof(double));

    st[log_index1]=3;
    pt[log_index1]=log_ready;

    flush(log_index1,*size); //flush

    st[log_index1]=1;
    st[log_index2]=0;        // log_index1 become new point log_index2 invalid

    t=log_index1;
    log_index1=log_index2;
    log_index2=t;

    /*
    if (logindex)
    {
        memcpy((void*)pmem1,(void*)log_data,(*size) * sizeof(double));
        p1index=log_ready;
        //flush(pmem1,*size);
    }
    else
    {
        memcpy((void*)pmem2,(void*)log_data,(*size) * sizeof(double));
        p2index=log_ready;
        //flush(pmem2,*size);
    }
    */
    logtime++;
    log_ready = 0;
  }
}

int recover_data()
{
    int res,i;
    /*
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
    */
    for (i=0;i<=3;i++)
    if (st[i]==1)
    {
        res=i;
        //data=pmem[i];
        break;
    }
    if (i==4) return -1;
    /*
  if (p1index < p2index)
  {
      memcpy((void*)data,(void*)pmem2,size * sizeof(double));
      res = p1index;
  } 
  else
  {
      memcpy((void*)data,(void*)pmem1,size * sizeof(double));
      res = p2index;
  }
    */
  return res;
}

void jacobi7_3(const int nx,const int ny, int nz, const double alpha,double* A0,const int timesteps,const double* B,const int ldb, double* Anext,const int ldc) 
{
  double fac;
  double *temp_ptr;
  int i, j, k, t, tstart=0;	
  fac = 6.0/(A0[0]*A0[0]);
  double end=0,now;
  double *l0, *lnext,*change;
  int l0n,l1n,changen;
  static int start=0;

  unsigned size = nx*ny*nz;

  if (!start)
  {
      memset(pmem,0,sizeof(pmem));
      //memset(pmem2,0,sizeof(pmem2));
      start=1;
      //printf("$$$$\n");
  }
  if (!tstart)
  {
      memcpy((void *)pmem[3],(void *)A0,size*sizeof(double));
      memcpy((void *)pmem[1],(void *)A0,size*sizeof(double));
      memset(st,0,sizeof(st));
      memset(pt,0,sizeof(pt));
  }
  if (log_id == 0)
  {
      i=recover_data();
      if (i==-1) tstart=0;
      else tstart=pt[i];
      if (!tstart)
      {
          l0 = pmem[0];
          lnext = pmem[1];
          l0n = 0;
          l1n = 1;
          st[0]=st[1]=-1;
          st[2]=0;
          st[3]=1;
          log_index1=2;
          log_index2=3;
      }
      else
      {
          l0 = pmem[(i+1)%4];
          lnext = pmem[(i+2)%4];
          l0n = (i+1)%4;
          l1n = (i+2)%4;
          st[(i+1)%4]=st[(i+2)%4]=-1;
          st[(i+3)%4]=0;
          //st[3]=1;
          log_index1=(i+3)%4;
          log_index2=i;
          //printf("xxxxxxxxx\n");
      }
      log_id = (pthread_t*) malloc(sizeof(pthread_t));
      pthread_create(log_id,NULL, log_all_data, &size); 
  }
  else if (!tstart)
  {
      l0 = pmem[0];
      lnext = pmem[1];
      l0n = 0;
      l1n = 1;
      st[0]=st[1]=-1;
      st[2]=0;
      st[3]=1;
      log_index1=2;
      log_index2=3;
  }
  /*
  printf("************\n");
  for (i=0;i<=3;i++) printf("%d ",st[i]);
  printf("\n");
  for (i=0;i<=3;i++) printf("%d ",pt[i]);
  printf("\n");
  printf("++++++++++++++\n");
  */
  //now=GetWallTime1();
  /*@;BEGIN(Nest1=Nest)@*/for (t = tstart; t < timesteps; t++)
      {
          int do_log = !(log_ready);
          change=l0;
          l0=lnext;
          lnext=change;

          changen=l0n;
          l0n=l1n;
          l1n=changen;
          //else {lnext = A0; l0 = Anext; }
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
          /*
          printf("************\n");
          for (i=0;i<=3;i++) printf("%d ",st[i]);
          printf("\n");
          for (i=0;i<=3;i++) printf("%d ",pt[i]);
          printf("\n");
          printf("++++++++++++++\n");*/
    if (do_log)
    {
        st[log_index1]=-1;
        changen=l0n;
        l0n=log_index1;
        l0=pmem[log_index1];
        log_index1=changen;

        log_ready=t+1;
    }
    /*
    printf("************\n");
    for (i=0;i<=3;i++) printf("%d ",st[i]);
    printf("\n");
    for (i=0;i<=3;i++) printf("%d ",pt[i]);
    printf("\n");
    printf("++++++++++++++\n");
    */


          /*
    if (do_log) {
        now=GetWallTime1();
        //printf("%.15f\n",now);
        if ((end==0)|((now-end)>0.05))
        {
            memcpy((void*)log_data,(void*)lnext,size * sizeof(double));
            log_ready=t+1;
            //printf("%.15f %.15f\n",end,now);
            //printf("now %lf end %lf",now,end);
            //printf("time between logging:%lf \n",now-end);
            end=now;
            //printf("%.15f\n %.15f\n",end,now);
        }
        }*/
          
  }
  memcpy((void*)lnext,(void*)A0,size * sizeof(double));
  printf("logtime: %d\n",logtime);
}

