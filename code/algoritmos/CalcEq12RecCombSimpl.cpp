#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793
#define Nmax 1

using namespace std;

void calcEva12RecIMP(int depth);
inline double sum(double * V, int length);

double * mFM, * W, * Fcoeff, * phiIM;
double Dhat, P0, alfa;
int Nc, offset, * nk;
complex <double> * Idet;
//                                  0  1   2      3        4       5       6    7       8 
// [Idet, I] = CalcEq12RecCombSimpl(t, f, mFMk, phiIMk, 2*pi*Fk, Fk/dff, alfa, beta2*L, P0);
//                                  t, f, mFM,  phiIM,  W,       Fcoeff, alfa, bL,      P0
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t M, N;
  int Nf, NN;
  double * IdetRe, * IdetIm, * Ir, *t, *f;
  double cft, sft, w;
  double bL;
  
  /* Check for proper number of arguments */
  if (nrhs != 9) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "São necessários 9 argumentos de entrada");
  } else if (nlhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "Esta função retorna apenas um argumento");
  }
  
  M = mxGetM(prhs[2]);
  N = mxGetN(prhs[2]);
  
  Nc = (int)N;
  alfa = mxGetScalar(prhs[6]);
  bL = mxGetScalar(prhs[7]);
  P0 = mxGetScalar(prhs[8]); 
    
  if(M != 1 || M != mxGetM(prhs[3]) || M != mxGetM(prhs[4]) || M != mxGetM(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de linhas diferente de 1 ou vetores de entrada com tamanhos diferentes!");
  else if(N < 2 || N != mxGetN(prhs[3]) || N != mxGetN(prhs[4]) || N != mxGetN(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de colunas diferente de Nc ou Vetores de entrada com tamanhos diferentes!");
    
  if(P0 < 0)
      mexErrMsgTxt("Erro!!! Argumentos de entrada inválidos");
  
  t = (double *) mxGetPr(prhs[0]);
  f = (double *) mxGetPr(prhs[1]);
  mFM = (double *) mxGetPr(prhs[2]);
  phiIM = (double *)mxGetPr(prhs[3]);
  W = (double *) mxGetPr(prhs[4]);
  Fcoeff = (double *) mxGetPr(prhs[5]);
  
  Nf = (int) mxGetN(prhs[1]);
  NN = (int) mxGetN(prhs[0]);
  
  // Realiza alguns cálculos antes de chamar a função recursiva
  offset = Nmax*((int)sum(Fcoeff, Nc));
 
  Dhat = -0.5*bL;
  
  Idet = new complex <double> [Nf];
  nk = new int[Nc];
  
  memset((void *)Idet, 0, Nf*sizeof(complex<double>));
  memset((void *)nk, 0, sizeof(int));
  
  // Chama função recursiva
  calcEva12RecIMP(0);
  
  // Passa resultado para o MatLab
  plhs[0] = mxCreateDoubleMatrix(1, Nf, mxCOMPLEX);
  plhs[1] = mxCreateDoubleMatrix(1, NN, mxREAL);
  
  IdetRe = mxGetPr(plhs[0]);
  IdetIm = mxGetPi(plhs[0]);
  
  Ir = mxGetPr(plhs[1]);
  // Ii = mxGetPi(plhs[1]);
    
  for(int i = 0; i < Nf; i++)
  {
      IdetRe[i] = P0*Idet[i].real();
      IdetIm[i] = P0*Idet[i].imag();
      w = 2*pi*f[i];
      for(int j = 0; j < NN; j++)
      {
          cft = cos(w*t[j]);
          sft = sin(w*t[j]);
          Ir[j] += IdetRe[i]*cft - IdetIm[i]*sft;
          // Ii[j] += IdetRe[i]*sft + IdetIm[i]*cft;
      }
  }
  
  // Limpeza
  delete Idet;
  delete nk;
  
  return;
}

void calcEva12RecIMP(int depth)
{
    if(depth == Nc-1)
	{
			int nn = offset;	
			double aux, phi;
            double S, P, u, theta, ProdJn;

            nk[depth] = -Nmax;
            depth++;
            
            for(int i = 0; i < depth; i++)
				nn += ((int)Fcoeff[i])*nk[i];
            	
			for(int k = -Nmax; k <= Nmax; k++)
			{
                if(k == 0 || nn == offset)
                {
                    nn += (int)Fcoeff[depth-1];
                    continue;
                }
                
				nk[depth-1] = k;

				// theta e exp(1j*nk.*dPSI)
                aux = 0;
                for(int i = depth-1; i >= 0; i--)
                    aux += W[i]*nk[i];
                   
                aux *= Dhat;
                
                S = 0.0;
                P = 1.0;
                phi = 0.0;
				for(int i = depth-1; i >= 0; i--)
				{
					theta = aux*W[i];
					u = 2*mFM[i]*sin(theta);
					theta = tan(theta);
                    
                    S += nk[i]/theta;
                    P *= (1 - abs(nk[i]) - u*nk[i]/2.0);
                    
                    phi += nk[i]*phiIM[i];
				}
				
                S /= alfa;
                
				Idet[nn] += P*(1 - S)*(std::complex <double> (cos(phi), sin(phi)));
                //mexPrintf("%d\n", nn);
                
				nn += (int)Fcoeff[depth-1];
			}
	}
	else
	{
		for(int k = -Nmax; k <= Nmax; k++)
		{
            if(k == 0)
                continue;
            else
            {                  
                nk[depth] = k;
                calcEva12RecIMP(depth+1);
            }
		}
	}
	return;
}

inline double sum(double * V, int length)
{
	double S = 0;
	for(int i = --length; i >= 0; i--)
		S += V[i];
	return S;
}