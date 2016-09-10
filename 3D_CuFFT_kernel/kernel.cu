
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>   
#include <stdlib.h>  
#include <math.h>  
#include <cuda.h> 
#include <cufft.h> 
#include <memory.h> 
#include<time.h>
#pragma comment( lib, "cufft.lib" )

#define N1 128 
#define N2 128
#define N3 16

#define CN3 ((int)N3/2+1) // half N3

void functf(double A[2][N1][N2][N3]);

void FFT3d_GPU(double A[2][N1][N2][N3]);

void main(void)
{

	static double A[2][N1][N2][N3];
	functf(A);

	FFT3d_GPU(A);
	printf("OK! \n");
}

void FFT3d_GPU(double A[2][N1][N2][N3])
{
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	const int    FFT_SIZE_R = N1 * N2 * N3  * sizeof(cufftDoubleReal);

	const int    FFT_SIZE_C = N1 * N2 * CN3 * sizeof(cufftDoubleComplex);

	cufftDoubleReal *h_oRealData = (cufftDoubleReal*)malloc(FFT_SIZE_R); 
	memset(h_oRealData, 0x00, FFT_SIZE_R);
	cufftDoubleComplex *h_otestComplexData = (cufftDoubleComplex*)malloc(FFT_SIZE_C);
	memset(h_otestComplexData, 0x00, FFT_SIZE_C);

	cufftDoubleReal *d_iRealData;
	cudaMalloc((void**)&d_iRealData, FFT_SIZE_R);
	cudaMemcpy(d_iRealData, A[0], FFT_SIZE_R, cudaMemcpyHostToDevice);

	cufftDoubleComplex *d_oComplexData;
	cudaMalloc((void**)&d_oComplexData, FFT_SIZE_C);

	cufftDoubleReal *d_oRealData;
	cudaMalloc((void**)&d_oRealData, FFT_SIZE_R);
	// /* Create a 3D FFT plan for D2Z */
	printf("cufft FFT,direct transform(x-->A), A:\n");
	cudaEventRecord(start, 0);
	cufftHandle planD2Z3D;
	cufftPlan3d(&planD2Z3D, N1, N2, N3, CUFFT_D2Z);
	cufftExecD2Z(planD2Z3D, (cufftDoubleReal*)d_iRealData, d_oComplexData);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("Kernel time: %.2f ms\n", elapsedTime);

	cudaMemcpy(h_otestComplexData, d_oComplexData, FFT_SIZE_C, cudaMemcpyDeviceToHost);
	cudaEventRecord(start, 0);
	//Create a 3D FFT plan for Z2D. 
	cufftHandle  planZ2D3D;
	cufftPlan3d(&planZ2D3D, N1, N2, N3, CUFFT_Z2D);
	//D2Z out of  place
	printf("\n\nAfter 3D C2R out of place : \n");
	printf("cuFFT IFFT,inverse transform(x-->A), A:\n");
	// Use the CUFFT plan to transform the signal out of place.  
	cufftExecZ2D(planZ2D3D, d_oComplexData, (cufftDoubleReal*)d_oRealData);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime1;
	cudaEventElapsedTime(&elapsedTime1, start, stop);
	printf("iCufft Kernel time: %.2f ms\n", elapsedTime1);
	printf("total Kernel time: %.2f ms\n", elapsedTime1 + elapsedTime);
	cudaMemcpy(h_oRealData, d_oRealData, FFT_SIZE_R, cudaMemcpyDeviceToHost);

//	printf("Total running time: %f seconds\n", (double)(totaltime1) / CLOCKS_PER_SEC);
	unsigned int i, j, k;
	FILE *fp;
	fp = fopen("fft3_cuda.d", "w");
	//system("cls");

	// origianl data
	printf("file copy\n");
	fprintf(fp, "the original data, AK:\n");
	for (k = 0; k<N3; k++)
		for (j = 0; j<N2; j++)
			for (i = 0; i<N1; i++)
				fprintf(fp, "%15.9f\n", A[0][i][j][k]);
	fprintf(fp, "\n");

	// direct transform
	fprintf(fp, "cufft FFT,direct transform(x-->A), A:\n");
	for (k = 0; k<CN3; k++)
		for (j = 0; j<N2; j++)
			for (i = 0; i<N1; i++)
				fprintf(fp, "%4u,%4u,%4u:%15.9f,%16.8e\n", i, j, k, h_otestComplexData[i*N2*CN3 + j*CN3 + k].x, h_otestComplexData[i*N2*CN3 + j*CN3 + k].y);
	fprintf(fp, "\n");

	//inverse transform

	fprintf(fp, "cu FFT IFFT,inverse transform(x-->A), A:\n");
	for (k = 0; k<N3; k++)
		for (j = 0; j<N2; j++)
			for (i = 0; i<N1; i++)
				fprintf(fp, "%4u,%4u,%4u:%15.9f\n", i, j, k, h_oRealData[i*N2*N3 + j*N3 + k] / (N1*N2*N3));
	fprintf(fp, "\n");

}
/*****************************************************/
void functf(double A[2][N1][N2][N3])
{
	unsigned int i, j, k;
	double tti, ttj, ttk, dti, dtj, dtk, A0i, A0j, A0k, A1i, A1j, A1k;
	dti = 0.1;  dtj = 2.0*dti;  dtk = dtj;
	A0i = exp(-N1*dti);
	A0j = exp(-N2*dtj);
	A0k = exp(-N3*dtk);

	for (i = 0; i<N1; i++)
	{
		tti = (double)i*dti;
		A1i = exp(-tti);
		for (j = 0; j<N2; j++)
		{
			ttj = (double)j*dtj;
			A1j = exp(-ttj);
			for (k = 0; k<N3; k++)
			{
				ttk = (double)k*dtk;
				A1k = exp(-ttk);

				A[0][i][j][k] = (A1i + A0i / A1i)*dti*(A1j + A0j / A1j)*dtj*(A1k + A0k / A1k)*dtk;
				A[1][i][j][k] = 0.0;
			}
		}
	}
}
