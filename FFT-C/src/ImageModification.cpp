#include <stdio.h>   
#include <stdlib.h>   
#include <math.h>  
#include "Header.h"
//#include "fft3r.c" 
#include "time.h" 
extern void fft3r(double A[2][N1][N2][N3], int ifft);
void functf(double A[2][N1][N2][N3]);
void main(void)
{
	clock_t totaltime;
	clock_t totaltime1;
	clock_t begin;
	clock_t end;
	unsigned int i, j, k;
	static double A[2][N1][N2][N3];
	FILE *fp;
	fp = fopen("fft3rm.d", "w");
	system("cls");
	functf(A);
	printf("The original data, Ak:\n");
	fprintf(fp, "The original data, Ak:\n");
	for (k = 0; k<N3; k++)
		for (j = 0; j<N2; j++)
			for (i = 0; i<N1; i++)
				fprintf(fp, "%4u,%4u,%4u: %15.9f\n", i, j, k, A[0][i][j][k]);
	begin = clock();
	fft3r(A, 1);
	end = clock();
	totaltime = end - begin;
	printf("fft: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("FFT, direct transform ( A -> x ), x:\n");
	fprintf(fp, "FFT, direct transform ( A -> x ), x:\n");
	for (k = 0; k<N3; k++)
		for (j = 0; j<N2; j++)
			for (i = 0; i<N1; i++)
				fprintf(fp, "%4u,%4u,%4u: %15.9f,%16.8e\n",
				i, j, k, A[0][i][j][k], A[1][i][j][k]);
	begin = clock();
	fft3r(A, -1);
	end = clock();
	totaltime1 = end - begin + totaltime;
	printf("inverse fft: %f seconds\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf("IFFT, inverse transform ( x -> A ), A:\n");
	fprintf(fp, "IFFT, inverse transform ( x -> A ), A:\n");
	for (k = 0; k<N3; k++)
		for (j = 0; j<N2; j++)
			for (i = 0; i<N1; i++)
				fprintf(fp, "%4u,%4u,%4u: %15.9f\n", i, j, k, A[0][i][j][k]);
	fclose(fp);
	printf("OK!\n");
	printf("total time: %f seconds\n", (double)(totaltime1) / CLOCKS_PER_SEC);
}
/*******************************************************************/
void functf(double A[2][N1][N2][N3])
{
	unsigned int i, j, k;
	double tti, ttj, ttk, dti, dtj, dtk, A0i, A0j, A0k, A1i, A1j, A1k;
	dti = 0.1;  dtj = 2.0*dti; dtk = dtj;
	A0i = exp(-N1*dti);  A0j = exp(-N2*dtj);     A0k = exp(-N3*dtk);
	for (i = 0; i<N1; i++)
	{
		tti = (double)i*dti; A1i = exp(-tti);
		for (j = 0; j<N2; j++)
		{
			ttj = (double)j*dtj; A1j = exp(-ttj);
			for (k = 0; k<N3; k++)
			{
				ttk = (double)k*dtk; A1k = exp(-ttk); A[1][i][j][k] = 0.0;
				A[0][i][j][k] = (A1i + A0i / A1i)*dti
					*(A1j + A0j / A1j)*dtj*(A1k + A0k / A1k)*dtk;
			}
		}
	}
}
