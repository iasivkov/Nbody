#include "Nbodygpu.cuh"
#include "Global_constants.h"

#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//ядро для подсчета взаимодействия
__device__ float4
bodyBodyInteraction(float4 bi, float4 bj, float4 ai, int NforDisc, int NforHalo, int NforBulge, int tile, int p)
{
	float3 r;
	int gtid = blockIdx.x * blockDim.x + threadIdx.x;
	float eps = 0.0f;
	
	/*
	if( gtid <  NforHalo + NforDisc  && gtid >= NforDisc )
	{
		ai.x += 0.0f;
		ai.y += 0.0f;
		ai.z += 0.0f;
		//ai.w += bj.w*bi.w/sqrt(distSqr);
		return ai;
	}
	*/

	if(gtid < NforDisc)
	{
		eps = 0.08f;
	}
	else if(gtid < NforHalo + NforDisc)
	{
		eps = 0.4f;
	}
	else if(gtid < NforHalo + NforDisc + NforBulge)
	{ 
		eps = 0.06f;
	}
	else 
	{
		eps = 0.08f;
	}
	
	r.x = bj.x - bi.x;
	r.y = bj.y - bi.y;
	r.z = bj.z - bi.z;
	// distSqr = dot(r_ij, r_ij) + EPS^2  
	float distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
	float R = sqrtf(distSqr);
	float u =R/eps; 
	float invDistCube = 1.0f;
	float fr=1.0f;
	float DistCube=R*R*R;
	if(u <= 1.0)
	{ 
		invDistCube = 1.0f/eps/eps/eps*(4.0f/3.0f - 6.0f/5.0f*u*u + 1.0f/20.0f*u*u*u*u*u);
	//	fr = -2.0f/eps*(1.0f/3.0f*u*u - 3.0f/20.0f*u*u*u*u);
	}
	else if(u <= 2.0)
	{
		invDistCube =1.0f/DistCube *(-1.0f/15.0f + 8.0f/3.0f*u*u*u - 3.0f*u*u*u*u + 6.0f/5.0f*u*u*u*u*u - 1.0f/6.0f*u*u*u*u*u*u);
	//	fr = -1.0f/15.0f*R - 1.0f/eps*(4.0f/3.0f*u*u - u*u*u + 3.0f/10.0f*u*u*u*u - 1.0f/30.0f*u*u*u*u*u) + 8.0f/5.0f*eps;
	}
	else 
	{
		invDistCube = 1.0f/DistCube;
		//	fr = 1/R;
	}


	float s = bj.w * invDistCube;
	// a_i =  a_i + s * r_ij [6 FLOPS]
	ai.x += r.x * s;
	ai.y += r.y * s;
	ai.z += r.z * s;
	//потенциал
	ai.w += bj.w*bi.w*fr;

	return ai;
}

//подсчет тайла
__device__ float4
tile_calculation(float4 myPosition, float4 accel,float4 *shPos, int NforDisc, int NforHalo, int NforBulge, int tile, int p)
{

	int i;
	
	for (i = 0; i < blockDim.x; i++) {
	accel = bodyBodyInteraction(myPosition,  shPos[i], accel, NforDisc,  NforHalo,  NforBulge, tile, p); 
}
return accel;
}


//-----------------------------------------------------------------------------

//подсчет ускорения

 void __global__ 
calculate_forces(float4 *devX, float4 *devA,float4 *devV, float4 *globalX, float4 *globalA, float4 *globalV, int N, int p, float timeStep,int NforDisc,int NforHalo,int NforBulge)
{

__shared__ float4 shPosition[512];
//globalX = devX;
//globalA = devA;
float4 myPosition;
int i, tile;
float4 acc = {0.0f, 0.0f, 0.0f, 0.0f};

int gtid = blockIdx.x * blockDim.x + threadIdx.x;

myPosition = devX[gtid];

for (i = 0, tile = 0; i < N; i += p, tile++) 
{
	int idx = tile * blockDim.x + threadIdx.x;
	shPosition[threadIdx.x] = devX[idx];

	__syncthreads();
	
	acc = tile_calculation(myPosition, acc, shPosition, NforDisc, NforHalo, NforBulge,tile,p);
	__syncthreads();
}

// Save the result in global memory for the integration step.

float4 acc4 = {acc.x, acc.y, acc.z, acc.w};
globalA[gtid] = acc4;

float b = devX[gtid].w * ((devV[gtid].x + acc.x * timeStep/2.0)*(devV[gtid].x + acc.x * timeStep/2.0) + (devV[gtid].y + acc.y * timeStep/2.0)*(devV[gtid].y + acc.y * timeStep/2.0) + (devV[gtid].z + acc.z * timeStep/2.0)*(devV[gtid].z + acc.z * timeStep/2.0))/2.0;
float velx = devV[gtid].x + globalA[gtid].x * timeStep ;
float vely = devV[gtid].y + globalA[gtid].y * timeStep ;
float velz = devV[gtid].z + globalA[gtid].z * timeStep ;
float4 vel4 = {velx, vely, velz, b };
globalV[gtid] = vel4;

float posx = devX[gtid].x + globalV[gtid].x * timeStep;
float posy = devX[gtid].y + globalV[gtid].y * timeStep;
float posz = devX[gtid].z + globalV[gtid].z * timeStep;
float4 pos = {posx, posy, posz, devX[gtid].w };
globalX[gtid] = pos;
}

//-----------------------------------------------------------------------------
// вычисление новых координат частиц
void Calculate(float timeStep, float4* X, float4* V, float4* A)
{		
		float4 *devX,*devA,*devV,*newX,*newA,*newV;

		cudaMalloc ( &devX, NParticles * sizeof(float4)); 
		cudaMalloc ( &devA, NParticles * sizeof(float4));
		cudaMalloc ( &devV, NParticles * sizeof(float4));
		cudaMalloc ( &newX, NParticles * sizeof(float4)); 
		cudaMalloc ( &newA, NParticles * sizeof(float4));
		cudaMalloc ( &newV, NParticles * sizeof(float4)); 
		
		checkCudaErrors(cudaMemcpy ( devX, X, NParticles * sizeof(float4), cudaMemcpyHostToDevice )); 
		checkCudaErrors(cudaMemcpy ( devA, A, NParticles * sizeof(float4), cudaMemcpyHostToDevice ));
		checkCudaErrors(cudaMemcpy( devV, V, NParticles * sizeof(float4), cudaMemcpyHostToDevice ));
		//printf("Iteration %f,\t%f\t\n", A[1].z,A[1].w);
		//printf("Velocity %f,\t%f\t%f\n", V[1000].x,V[1000].y,V[1000].z);
		// вычисляем новую скорость и позицию
	
		//cudaPrintfInit(25600000);
		
		//system("pause");
		printf("Velocity %f,\t%f\t%f\n", V[1000].x,V[1000].y,V[1000].z);
		printf("Position %f,\t%f\t%f\n", X[1000].x,X[1000].y,X[1000].z);
		
		calculate_forces<<<dim3((int)NParticles/BLOCK_SIZE),dim3(BLOCK_SIZE)>>> ( devX,  devA, devV, newX, newA, newV, NParticles, BLOCK_SIZE ,timeStep, NforDisc, NforHalo, NforBulge);
		printf("Tam-pam:", cudaGetErrorString(cudaPeekAtLastError()));
		cudaError_t f = cudaDeviceSynchronize();
		printf("Velocity %s\n",  cudaGetErrorString( f ));
		
		cudaMemcpy ( X, newX, NParticles * sizeof(float4), cudaMemcpyDeviceToHost ); 
		cudaMemcpy ( A, newA, NParticles * sizeof(float4), cudaMemcpyDeviceToHost ); 
		cudaMemcpy ( V, newV, NParticles * sizeof(float4), cudaMemcpyDeviceToHost );
		
		
		float E=0;
		for(int i = 0; i<NParticles; i++)
		{

			E+= A[i].w + V[i].w;



		/*	V[i].x += A[i].x * TimeStep ;
			V[i].y += A[i].y * TimeStep ;
			V[i].z += A[i].z * TimeStep ;
			
			X[i].x += V[i].x * TimeStep;
			X[i].y += V[i].y * TimeStep;
			X[i].z += V[i].z * TimeStep;
		*/
		}
		//printf("Position %f\r", E);
	//	printf("Velocity %f,\t%f\t%f\n", V[1000].x,V[1000].y,V[1000].z);
		// system("pause");
		//printf("\r");
		cudaFree (devX);
		cudaFree (devA);
		cudaFree (devV);
		cudaFree (newX);
		cudaFree (newA);
		cudaFree (newV);
	
    
}

void    WriteData(int i, int k, float4* X, float4* V, float4* A)
{
    char FileName[128];
	char FileName2[128];
	char FileName3[128];
	char FileName4[128];
	//char FileName5[128];

	if(k>0)
	{
		sprintf(FileName, "OutFiles_%i\\disc_%i.dat",k, i);
		sprintf(FileName2, "OutFiles_%i\\HALO_%i.dat", k,i);
		sprintf(FileName3, "OutFiles_%i\\bulge_%i.dat",k, i);
		sprintf(FileName4, "OutFiles_%i\\satellite_%i.dat",k, i);
		//sprintf(FileName5, "VR_%i\\Vr_%i.dat",j, i);
	}
	else
	{
		sprintf(FileName, "OutFiles\\disc_%i.dat", i);
		sprintf(FileName2, "OutFiles\\HALO_%i.dat", i);
		sprintf(FileName3, "OutFiles\\bulge_%i.dat", i);
		sprintf(FileName4, "OutFiles\\satellite_%i.dat", i);
	}


    FILE * out = fopen(FileName, "w");
	FILE * out2 = fopen(FileName2, "w");
	FILE * out3 = fopen(FileName3, "w");
	FILE * out4 = fopen(FileName4, "w");
	//FILE * out5 = fopen(FileName5, "w");
	/*float QdevSig[501];
	for(int j=0;j<501;j++)
	{
		float r = (j+1)*R/500.0;
		float K = sqrt(3.0*(dQ(r,0) + dQtot(r))/r + (ddQ(r,0) + ddQtot(r)));
		QdevSig[j] = K/(3.36 * G * rho0 * 2.0 * z0 * exp(-r/h));
	}
*/
    for(int j=0; j<NParticles; j++)
	{
		if(j<NforDisc)
		{
			fprintf(out, "%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",j, X[j].w, X[j].x, X[j].y, X[j].z, V[j].x, V[j].y, V[j].z);
			/*float r = sqrt(X[j].x * X[j].x + X[j].y * X[j].y);
			int l=r/R*500;
			//printf("halo: %f\t%i",r/R*1000.0,l);
			if(j<NforDisc/4)fprintf(out5, "%f\t%f\t%f\t%f\n", r, -V[j].x*X[j].x/r - V[j].y*X[j].y/r,QdevSig[l], -V[j].x*X[j].y/r + V[j].y*X[j].x/r );
			*/
		}
		
		if(j<(NforDisc + NforHalo) && j>=NforDisc)
		{
				fprintf(out2,"%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",j, X[j].w, X[j].x, X[j].y, X[j].z, V[j].x, V[j].y, V[j].z);
		}
			
		if(j <(NforDisc + NforHalo + NforBulge)&& j>=(NforDisc + NforHalo)) 
		{
						fprintf(out3,"%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",j, X[j].w, X[j].x, X[j].y, X[j].z, V[j].x, V[j].y, V[j].z);
		}
		if(j< (NforSatellite + NforDisc + NforHalo + NforBulge) && j>=(NforDisc + NforHalo + NforBulge))
		{
						fprintf(out4,"%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",j, X[j].w, X[j].x, X[j].y, X[j].z, V[j].x, V[j].y, V[j].z);
		}
						
	}
    fclose(out);
	fclose(out2);
	fclose(out3);
	fclose(out4);
	//fclose(out5);
	

}


namespace NbodyCu{
    void Run(float TimeStep, float4* X, float4* V, float4* A)
    {
        float GpuTime;
		printf("Start");
		int k = 0;

        WriteData(0, k, X, V, A);
        cudaEvent_t start,stop;

        cudaEventCreate(&start);
        cudaEventCreate(&stop);	

        cudaEventRecord(start,0);
        
        for (int i=0; i<500; i++)
        {
            if(i==0) Calculate(TimeStep/2.0, X, V, A);
            else Calculate(TimeStep, X, V, A);
            
            printf("Iteration %i,\n", i);
            WriteData(i + 1, k, X, V, A);
        }
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&GpuTime,start,stop);

            printf("Gpu Time is: %f\n", GpuTime/1000.0);

        cudaEventDestroy(start);
        cudaEventDestroy(stop);
    }
}