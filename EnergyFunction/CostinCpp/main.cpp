#include <iostream>
#include <math.h>
//#include <complex>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <memory.h>
#define CopyVector(a,b) memcpy((a),(b),8*sizeof(double))

#define N_MIC 3
#define N_SRC 2

using namespace std;

double testWhatever(void)
{
    /*complex<double> z(2, -5);
    const complex<double> I(0,1);
    cout << z << endl;
    cout << z.imag() << endl;
    cout << z.real() << endl;
    cout << conj(z) << endl;
    cout << exp(I*0.725) << endl;*/
    return 1.1;
}

void randomSqArray(double A[][N_MIC])
{
    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            A[i][j] = (double)rand() / RAND_MAX;
        }
    }
}

template< typename T, size_t N, size_t M >
void printArray( T(&theArray)[N][M]  ) {
    cout.setf(ios::fixed); cout.setf(ios::showpoint); cout.precision(6);
    for ( size_t x = 0; x < N; x ++ ) {
        for ( size_t y = 0; y < M; y++ ) {
            cout << theArray[x][y] << " ";
        }
        cout << endl;
    }
}

void doCalcs(double* Xm, double* Xs, double CSMr[][N_MIC], double CSMi[][N_MIC])
{
    //const complex<double> I(0,1);
    //complex<double> CSMmodel[N_MIC][N_MIC];
    double CSMmodelr[N_MIC][N_MIC] = {{0}}, CSMmodeli[N_MIC][N_MIC] = {{0}};

    double Cost = 0;
    double r[N_SRC*N_MIC];

    for (int i = 0; i < N_SRC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            r[N_MIC*i+j] = sqrt( pow(Xs[4*i]-Xm[3*j],2) + pow(Xs[4*i+1]-Xm[3*j+1],2) + pow(Xs[4*i+2]-Xm[3*j+2],2));
            cout << i + 1 << "," << j+1 << ": " << r[N_MIC*i+j] << endl;
        }
    }

    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            for (int k = 0; k < N_SRC; k++)
            {
                //CSMmodel[i][j] += 0.5*(Xs[4*k+3]/(4*4*M_PI*M_PI))*(exp(-I*2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
                CSMmodelr[i][j] += 0.5*(Xs[4*k+3]/(4*4*M_PI*M_PI))*(cos(-2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
                CSMmodeli[i][j] += 0.5*(Xs[4*k+3]/(4*4*M_PI*M_PI))*(sin(-2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
            }
        }
    }

    printArray(CSMmodelr);
    printArray(CSMmodeli);

    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            //Cost += pow(CSMmodel[i][j].real()-CSMr[i][j],2) +
            //        pow(CSMmodel[i][j].imag()-CSMi[i][j],2);
            Cost += pow(CSMmodelr[i][j]-CSMr[i][j],2) +
                    pow(CSMmodeli[i][j]-CSMi[i][j],2);
        }
    }
    cout.precision(15);
    cout << "DA ENERGY ISSSS....: " << Cost << endl;
}

void doCalcs2(double* mic_positions, double* trial, double CSMrr[], double CSMii[])
{
    double Cost = 0;
    double *r = new double[N_SRC*N_MIC];

    /*double** CSMmodelr = new double*[N_MIC];
    double** CSMmodeli = new double*[N_MIC];
    for(int i = 0; i < N_MIC; ++i) {
        CSMmodelr[i] = new double[N_MIC];
        CSMmodeli[i] = new double[N_MIC];
    }*/

    double CSMmodelr[N_MIC][N_MIC] = {{0}}, CSMmodeli[N_MIC][N_MIC] = {{0}};

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
                CSMmodelr[i][j] += 0.5*(trial[4*k+3]/(4*4*M_PI*M_PI))*(cos(-2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
                CSMmodeli[i][j] += 0.5*(trial[4*k+3]/(4*4*M_PI*M_PI))*(sin(-2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
            }
        }
    }

    //printArray(CSMmodelr);
    //printArray(CSMmodeli);

    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            Cost += pow(CSMmodelr[i][j]-CSMrr[N_MIC*j+i],2) +
                    pow(CSMmodeli[i][j]-CSMii[N_MIC*j+i],2);
        }
    }

    cout.precision(15);
    cout << "DA ENERGY ISSSS....: " << Cost << endl;
}

int main()
{
    double params[8] = {1.1, 2.1, 1.0, 1.0,
                        0.5, 1.1, 1.0, 2.0};
    double micpos[9] = {0.5, 0.5, 0, 1.0, -0.4, 0, -0.2, -0.3, 0};
    double destiny[8];

    //double CSMr[N_MIC][N_MIC];
    //double CSMi[N_MIC][N_MIC];
    //randomSqArray(CSMr); randomSqArray(CSMi);

    double CSMr[N_MIC][N_MIC] = {{0.08,0.05,0.03},{0.0408,0.03,0},{0.03,0.02,0.03}};
    double CSMi[N_MIC][N_MIC] = {{0,0.8,-0.004},{-0.8,0,-0.01},{0.004,0.01,0}};
    double CSMrr[N_MIC*N_MIC] = {0.08,0.0408,0.03,0.05,0.03,0.02,0.03,0,0.03};
    double CSMii[N_MIC*N_MIC] = {0,-0.8,0.004,0.8,0,0.01,-0.004,-0.01,0};

    //printArray(CSMi);
    doCalcs(micpos, params, CSMr, CSMi);

    CopyVector(destiny,params);
    printf("\t%.2f\t%.2f\t%.2f\n\n", destiny[0] , destiny[3] , destiny[7] );

    doCalcs2(micpos, params, CSMrr, CSMii);

    return 0;
}
