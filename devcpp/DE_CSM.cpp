// compile in matlab as: mex -I"O:\Array Benchmark\devcpp" DE_CSM.cpp DESolver.cpp

#define MAX_GENERATIONS 10
#define N_POP 5
#define p_c 0.75
#define F 0.4

#define  M_PI  3.1415926535897932

#include <mex.h>
#include <DESolver.h>
#include <math.h>

double frequency;
double speed_of_sound;
double *CSM_real;        // N_MIC by N_MIC
double *CSM_imag;        // N_MIC by N_MIC
double *mic_positions;   // N_MIC by 3
size_t N_DIM;
size_t N_MIC;

// Acoustic source problem using CSM
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
    //double CSMmodelr[N_MIC][N_MIC] = {{0}};
    //double CSMmodeli[N_MIC][N_MIC] = {{0}};
    //double r[N_DIM*N_MIC];
    double *r = new double[N_DIM*N_MIC];
    
    double** CSMmodelr = new double*[N_MIC];
    double** CSMmodeli = new double*[N_MIC];
    for(int i = 0; i < N_MIC; ++i) {
        CSMmodelr[i] = new double[N_MIC];
        CSMmodeli[i] = new double[N_MIC];
    }
    
    for (int i = 0; i < N_DIM; i++)
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
            for (int k = 0; k < N_DIM; k++)
            {
                CSMmodelr[i][j] += 0.5*(trial[4*k+3]/(4*4*M_PI*M_PI))*(cos(-2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
                CSMmodeli[i][j] += 0.5*(trial[4*k+3]/(4*4*M_PI*M_PI))*(sin(-2.0*M_PI*1000.0*(r[k*N_MIC+i]-r[k*N_MIC+j])/343.0)/(r[k*N_MIC+i]*r[k*N_MIC+j]));
            }
        }
    }
    
    for (int i = 0; i < N_MIC; i++)
    {
        for (int j = 0; j < N_MIC; j++)
        {
            Cost += pow(CSMmodelr[i][j]-CSM_real[N_MIC*i+j],2) +
                    pow(CSMmodeli[i][j]-CSM_imag[N_MIC*i+j],2);
        }
    }
    
    if (count++ % nPop == 0){
		//mexPrintf("Gen: %d\tEnergy: %lf\n", count / nPop + 1, Energy());
		mexPrintf("GEN%d\t%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", count / nPop + 1, trial[0], trial[1], trial[2], trial[3], trial[4], trial[5], trial[6], trial[7]);
    }
    delete [] r;
    return(Cost);
}

/* The computational routine */
void DE_CSM(double f, double c, double *CSM_real, double *CSM_imag, 
        double *lower_bounds, double *upper_bounds, 
        double *mic_positions, double *outBestCosts, double *outBestParams)
{
  int i;
  AcousticSourceSolver solver(N_DIM,N_POP);
  
  solver.Setup(lower_bounds,upper_bounds,stRand1Exp,F,p_c);
  solver.Solve(MAX_GENERATIONS);
  double *solution = solver.Solution();
  
  for (i=0; i<8; i++) {
    mexPrintf("\t%.6f\n", solution[i]);
  }
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Variable declarations
    // input scalar
    size_t M_CSM, N_CSM;
    
    // input arrays
    double *upper_bounds;    // vector with upper bounds
    double *lower_bounds;    // vector with lower bounds
    double *outBestCosts, *outBestParams;
    
    if(nrhs != 7) {
    mexErrMsgIdAndTxt("MyToolbox:DE_CSM:nrhs",
                      "Eight inputs required.");
    }
    if(nlhs != 2) {
    mexErrMsgIdAndTxt("MyToolbox:DE_CSM:nlhs",
                      "Two outputs required.");
    }
    
    frequency = mxGetScalar(prhs[0]);
    speed_of_sound = mxGetScalar(prhs[1]);
    mexPrintf("\tFreq: %.2f\n", frequency);
    mexPrintf("\tSpeed of sound: %.2f\n", speed_of_sound);
    
    CSM_real = mxGetPr(prhs[2]);
    CSM_imag = mxGetPr(prhs[3]);
    M_CSM = mxGetM(prhs[2]);
    N_CSM = mxGetN(prhs[2]);
    mexPrintf("\tCSMr: %d\t%d\t%.4f\n", M_CSM, N_CSM, CSM_real[1]);
    mexPrintf("\tCSMi: %d\t%d\t%.4f\n", M_CSM, N_CSM, CSM_imag[13]);
    
    lower_bounds = mxGetPr(prhs[4]);
    upper_bounds = mxGetPr(prhs[5]);
    N_DIM = mxGetM(prhs[4]);
    mexPrintf("\tBounds: %d\t%.2f\t%.2f\n", N_DIM, lower_bounds[1], upper_bounds[2]);
    
    mic_positions = mxGetPr(prhs[6]);
    N_MIC =  mxGetM(prhs[6]);
    mexPrintf("\tMics: %d\t%.4f\n", N_MIC, mic_positions[6]);

    plhs[0] = mxCreateDoubleMatrix(1,MAX_GENERATIONS,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,N_DIM,mxREAL);
    outBestCosts = mxGetPr(plhs[0]);
    outBestParams = mxGetPr(plhs[1]);
    
    DE_CSM(frequency,speed_of_sound,CSM_real,CSM_imag,
            lower_bounds,upper_bounds, mic_positions,
            outBestCosts, outBestParams);
    
//     plhs[0][0] = 1.0;
//     plhs[0][1] = -1.0;
//     plhs[0][2] = 0.0;
//     plhs[0][3] = 1.1;
//     plhs[0][4] = 1.5;
    
/* code here */
    
//     /* get the value of the scalar input  */
//     multiplier = mxGetScalar(prhs[0]);
//     /* create a pointer to the real data in the input matrix  */
//     inMatrix = mxGetPr(prhs[1]);
//     /* get dimensions of the input matrix */
//     ncols = mxGetN(prhs[1]);
//     
//     /* create the output matrix */
//     plhs[0] = mxCreateDoubleMatrix(1,ncols,mxREAL);
//     /* get a pointer to the real data in the output matrix */
//     outMatrix = mxGetPr(plhs[0]);
//     /* call the computational routine */
//     DE_CSM(multiplier,inMatrix,outMatrix,ncols);
}