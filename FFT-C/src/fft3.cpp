#include <stdio.h>   
#include <stdlib.h>   
#include <math.h>   
#include "Header.h"
/*#define eps 1.0e-10*/   
unsigned int inverseq( unsigned int i1, unsigned int imax );   
void fft31( double A[2][N1], int ifft );   
void fft32( double A[2][N2], int ifft );   
void fft33( double A[2][N3], int ifft );   
void inital( void );   
static int ifirst=0,sign=-1;   
static unsigned M1,M2,M3,ijup[16],ikup[16],jjup[16],jkup[16],   
                kjup[16],kkup[16];   
static double pi,wri[N1],wii[N1],wrj[N2],wij[N2],wrk[N3],wik[N3];   
/*********** initial:   ********************************************/   
void initial ( void )   
{   
unsigned int i,j,ni;   
double unit,tt;   
ni=N1; M1=0;   
while( ni>1 ) { ni=ni/2; M1=M1+1; }   
ni=1;   
for(j=0;j<M1;j++) ni=ni*2;   
if( (ni-N1) != 0 )   
  { printf(" N1 isn't compatible!\n Strike any key to exit!\n");   
    getchar(); exit(1); }   
ni=N2; M2=0;   
while( ni>1 ) { ni=ni/2; M2=M2+1; }   
ni=1;   
for(j=0;j<M2;j++) ni=ni*2;   
if( (ni-N2) != 0 )   
  { printf(" N2 isn't compatible!\n Strike any key to exit!\n");   
    getchar(); exit(1); }   
ni=N3; M3=0;   
while( ni>1 ) { ni=ni/2; M3=M3+1; }   
ni=1;   
for(j=0;j<M3;j++) ni=ni*2;   
if( (ni-N3) != 0 )   
  { printf(" N3 isn't compatible!\n Strike any key to exit!\n");   
    getchar(); exit(1); }   
ikup[0]=N1/2;       ijup[0]=1;   /* kup[l]=2^(M1-l-1) */   
for(j=1;j<M1;j++)  { ikup[j]=ikup[j-1]/2; ijup[j]=ijup[j-1]*2; }   
pi=4.0*atan((double)1.0);  unit=2.0*pi/((double)N1);  tt=0.0;   
for(i=0;i<N1;i++)   
  { wri[i]=cos(tt); wii[i]=sin(tt); tt=tt+unit; }   
jkup[0]=N2/2;       jjup[0]=1;   /* jkup[l]=2^(M2-l-1) */   
for(j=1;j<M2;j++)  { jkup[j]=jkup[j-1]/2; jjup[j]=jjup[j-1]*2; }   
unit=2.0*pi/((double)N2);  tt=0.0;   
for(i=0;i<N2;i++)   
  { wrj[i]=cos(tt); wij[i]=sin(tt); tt=tt+unit; }   
kkup[0]=N3/2;       kjup[0]=1;   /* kkup[l]=2^(M3-l-1) */   
for(j=1;j<M3;j++)  { kkup[j]=kkup[j-1]/2; kjup[j]=kjup[j-1]*2; }   
unit=2.0*pi/((double)N3);  tt=0.0;   
for(i=0;i<N3;i++)   
  { wrk[i]=cos(tt); wik[i]=sin(tt); tt=tt+unit; }   
}   
/**   row: ifft=1 <--> FFT,ifft=-1 <--> IFFT ***************************/   
void fft31( double A[2][N1], int ifft )   
{   
unsigned int l,k,i,j,jl0,jl0k,jl1k,reseq,ij;   
double wr,wi,ar,ai,x[2][N1];   
if(ifirst==0) { initial(); ifirst=1; }   
if(ifft==-1) for(i=0;i<N1;i++)  A[1][i]=-A[1][i];   
for(l=0;l<M1;l++)   
  { for(j=0;j<ijup[l];j++)   
      { reseq=inverseq(j,l);   
        ij=ikup[l]*reseq;  jl0=2*ikup[l]*j;   
        wr=wri[ij];       wi=sign*wii[ij];   
        for(k=0;k<ikup[l];k++)   
          { jl0k=jl0+k;        jl1k=jl0+ikup[l]+k;   
            ar=A[0][jl1k]*wr-A[1][jl1k]*wi;   
            ai=A[0][jl1k]*wi+A[1][jl1k]*wr;   
            A[0][jl1k]=A[0][jl0k]-ar;   
            A[1][jl1k]=A[1][jl0k]-ai;   
            A[0][jl0k]=A[0][jl0k]+ar;   
            A[1][jl0k]=A[1][jl0k]+ai; }  }  }   
for(i=0;i<N1;i++)   
  { j=inverseq(i,M1);  x[0][i]=A[0][j];  x[1][i]=A[1][j]; }   
for(i=0;i<N1;i++)  { A[0][i]=x[0][i]; A[1][i]=x[1][i]; }   
if( ifft==-1 )   
  for(i=0;i<N1;i++)   
    { A[0][i]=A[0][i]/((double)N1); A[1][i]=-A[1][i]/((double)N1); }   
}   
/*column:  ***** ifft=1 <--> FFT, ifft=-1 <--> IFFT ******************/   
void fft32( double A[2][N2], int ifft )   
{   
unsigned int l,k,i,j,jl0,jl0k,jl1k,reseq,ij;   
double wr,wi,ar,ai,x[2][N2];   
if(ifirst==0) { initial(); ifirst=1; }   
if(ifft==-1) for(i=0;i<N2;i++)  A[1][i]=-A[1][i];   
for(l=0;l<M2;l++)   
  { for(j=0;j<jjup[l];j++)   
      { reseq=inverseq(j,l);   
        ij=jkup[l]*reseq;  jl0=2*jkup[l]*j;   
        wr=wrj[ij];       wi=sign*wij[ij];   
        for(k=0;k<jkup[l];k++)   
          { jl0k=jl0+k;        jl1k=jl0+jkup[l]+k;   
            ar=A[0][jl1k]*wr-A[1][jl1k]*wi;   
            ai=A[0][jl1k]*wi+A[1][jl1k]*wr;   
            A[0][jl1k]=A[0][jl0k]-ar;   
            A[1][jl1k]=A[1][jl0k]-ai;   
            A[0][jl0k]=A[0][jl0k]+ar;   
            A[1][jl0k]=A[1][jl0k]+ai; }  }  }   
for(i=0;i<N2;i++)   
  { j=inverseq(i,M2);  x[0][i]=A[0][j];  x[1][i]=A[1][j]; }   
for(i=0;i<N2;i++)  { A[0][i]=x[0][i]; A[1][i]=x[1][i]; }   
if( ifft==-1 )   
  for(i=0;i<N2;i++)   
    { A[0][i]=A[0][i]/((double)N2); A[1][i]=-A[1][i]/((double)N2); }   
}   
/************************************************************************/   
void fft33( double A[2][N3], int ifft )   
{   
unsigned int l,k,i,j,jl0,jl0k,jl1k,reseq,ij;   
double wr,wi,ar,ai,x[2][N3];   
if(ifirst==0) { initial(); ifirst=1; }   
if(ifft==-1) for(i=0;i<N3;i++)  A[1][i]=-A[1][i];   
for(l=0;l<M3;l++)   
  { for(j=0;j<kjup[l];j++)   
      { reseq=inverseq(j,l);   
        ij=kkup[l]*reseq;  jl0=2*kkup[l]*j;   
        wr=wrk[ij];       wi=sign*wik[ij];   
        for(k=0;k<kkup[l];k++)   
          { jl0k=jl0+k;        jl1k=jl0+kkup[l]+k;   
            ar=A[0][jl1k]*wr-A[1][jl1k]*wi;   
            ai=A[0][jl1k]*wi+A[1][jl1k]*wr;   
            A[0][jl1k]=A[0][jl0k]-ar;   
            A[1][jl1k]=A[1][jl0k]-ai;   
            A[0][jl0k]=A[0][jl0k]+ar;   
            A[1][jl0k]=A[1][jl0k]+ai; }  }  }   
for(i=0;i<N3;i++)   
  { j=inverseq(i,M3);  x[0][i]=A[0][j];  x[1][i]=A[1][j]; }   
for(i=0;i<N3;i++)  { A[0][i]=x[0][i]; A[1][i]=x[1][i]; }   
if( ifft==-1 )   
  for(i=0;i<N3;i++)   
    { A[0][i]=A[0][i]/((double)N3); A[1][i]=-A[1][i]/((double)N3); }   
}   
/********************************************************************/   
unsigned int inverseq( unsigned int i1, unsigned int imax )   
{   
unsigned int j,ii,i0,inv;   
i0=i1;  inv=0;   
for(j=0;j<imax;j++)  { ii=i0%2; inv=2*inv+ii; i0=(i0-ii)/2; }   
return(inv);   
}   