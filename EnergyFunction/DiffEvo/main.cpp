// Differential Evolution Test Program
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "DESolver.h"
#include <math.h>
#include <iterator>
#include <vector>


using namespace std;

double frequency;
double speed_of_sound;
double *CSM_real, *CSM_imag, *mic_positions;
size_t N_MIC;
size_t N_DIM;
size_t N_SRC;


// Polynomial fitting problem
class AcousticSourceSolver : public DESolver
{
    public:
        AcousticSourceSolver(int dim,int pop) : DESolver(dim,pop), count(0) {;}
        double EnergyFunction(double trial[],bool &bAtSolution);

    private:
        int count;
};

double AcousticSourceSolver::EnergyFunction(double *trial,bool &bAtSolution)
{
    double Cost = 0;
    double *r = new double[N_SRC*N_MIC];

    /*double** CSMmodelr = new double*[N_MIC];
    double** CSMmodeli = new double*[N_MIC];
    for(int i = 0; i < N_MIC; ++i) {
        CSMmodelr[i] = new double[N_MIC];
        CSMmodeli[i] = new double[N_MIC];
    }*/

    double CSMmodelr[32][32] = {{0}}, CSMmodeli[32][32] = {{0}};

    for (int i = 0; i < N_SRC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            r[N_MIC*i+j] = sqrt( pow(trial[4*i]-mic_positions[3*j],2) + pow(trial[4*i+1]-mic_positions[3*j+1],2) + pow(trial[4*i+2]-mic_positions[3*j+2],2));
        }
    }

    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            for (int k = 0; k < N_SRC; k++)
            {
                CSMmodelr[i][j] += 0.5*(trial[4*k+3]/(4*4*M_PI*M_PI))*(cos(-2.0*M_PI*frequency*(r[k*N_MIC+i]-r[k*N_MIC+j])/speed_of_sound)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
                CSMmodeli[i][j] += 0.5*(trial[4*k+3]/(4*4*M_PI*M_PI))*(sin(-2.0*M_PI*frequency*(r[k*N_MIC+i]-r[k*N_MIC+j])/speed_of_sound)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
            }
        }
    }

    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            Cost += pow(CSMmodelr[i][j]-CSM_real[N_MIC*j+i],2) +
                    pow(CSMmodeli[i][j]-CSM_imag[N_MIC*j+i],2);
        }
    }

    if (count++ % nPop == 0){
        printf("%.5f\t%.2f\n", M_PI, speed_of_sound);
		printf("Gen: %d\tEnergy: %.10f\n", count / nPop + 1, Energy());
		printf("GEN%d\t%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", count / nPop + 1, trial[0], trial[1], trial[2], trial[3], trial[4], trial[5], trial[6], trial[7]);
    }
    delete [] r;
    return(Cost);
}


double* loadData(const char* vijlnaam)
{
    size_t size;
    streampos begin,end;

    ifstream A(vijlnaam, ios::in|ios::binary);

    if (A) {
        begin = A.tellg();
        A.seekg (0, ios::end);
        end = A.tellg();
	}
	size = end - begin;

    char *memblock = new char [size];
    A.seekg(0, ios::beg);
    A.read(memblock, end-begin);
    A.close();

    double* out = (double*)memblock;
    return out;
}

#define N_POP 128
#define MAX_GENERATIONS	600
int main(void)
{
	int i;

    frequency = 2000.0;
    speed_of_sound = 343.0;
    N_MIC = 32;
    N_DIM = 8;
    N_SRC = 2;
    CSM_real = loadData("csmreal.bin");
    CSM_imag = loadData("csmimag.bin");
    mic_positions = loadData("micpos.bin");

	double min[] = {-1,-1,0,0,-1,-1,0,0};
	double max[] = {1,1,2,0.5,1,1,2,0.5};

	AcousticSourceSolver solver(N_DIM,N_POP);

	solver.Setup(min,max,stMyCase,0.4,0.75);

	printf("Calculating...\n\n");
	solver.Solve(MAX_GENERATIONS, min, max);

	double *solution = solver.Solution();

	printf("\n\nBest Coefficients:\n");
	for (i=0;i<N_DIM;i++)
		printf("[%d]: %lf\n",i,solution[i]);

}
