#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793

using namespace std;

void calcEva12RecIMP(int depth);
inline int sum(int * V, int length);
inline double sum(double * V, int length);

double * mIM, * mFM, * dPSI, * W, * Fcoeff, * phiIM;
double Dhat, P0;
int Nc, Nmax, offset, * nk;
complex <double> * epdPSI, * endPSI, * Idet;
//                                0  1   2      3    4       5       6       7       8       9    10
// [Idet, I] = CalcEq12RecCombCpp(t, f, mIMk, mFMk, dPSIk, phiIMk, 2*pi*Fk, Fk/dff, beta2*L, P0, Nmax);
//                                t, f, mIM,  mFM,  dPSI,  phiIM,  W,       Fcoeff, bL,      P0, Nmax
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t M, N;
  int Nf, NN;
  double * IdetRe, * IdetIm, * Ir, *t, *f;
  double cft, sft, w;
  double bL;
  
  /* Check for proper number of arguments */
  if (nrhs != 11) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "São necessários 11 argumentos de entrada");
  } else if (nlhs != 2) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "Esta função retorna apenas um argumento");
  }
  
  M = mxGetM(prhs[6]);
  N = mxGetN(prhs[6]);
  
  Nc = (int)N;
  bL = mxGetScalar(prhs[8]);
  P0 = mxGetScalar(prhs[9]);
  Nmax = (int)mxGetScalar(prhs[10]);
    
  if(M != 1 || M != mxGetM(prhs[2]) || M != mxGetM(prhs[3]) || M != mxGetM(prhs[4]) ||
          M != mxGetM(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de linhas diferente de 1 ou vetores de entrada com tamanhos diferentes!");
  else if(N < 2 || N != mxGetN(prhs[2]) || N != mxGetN(prhs[3]) || N != mxGetN(prhs[4]) ||
          N != mxGetN(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de colunas diferente de Nc ou Vetores de entrada com tamanhos diferentes!");
    
  if(P0 < 0 || Nmax < 0)
      mexErrMsgTxt("Erro!!! Argumentos de entrada inválidos");
  
  t = (double *) mxGetPr(prhs[0]);
  f = (double *) mxGetPr(prhs[1]);
  mIM = (double *) mxGetPr(prhs[2]);
  mFM = (double *) mxGetPr(prhs[3]);
  dPSI = (double *) mxGetPr(prhs[4]);
  phiIM = (double *)mxGetPr(prhs[5]);
  W = (double *) mxGetPr(prhs[6]);
  Fcoeff = (double *) mxGetPr(prhs[7]);
  
  Nf = (int) mxGetN(prhs[1]);
  NN = (int) mxGetN(prhs[0]);
  
  // Realiza alguns cálculos antes de chamar a função recursiva
  offset = Nmax*((int)sum(Fcoeff, Nc));

  epdPSI = new complex <double> [Nc];
  endPSI = new complex <double> [Nc];
  for(int i = 0; i < Nc; i++)
  {
    epdPSI[i] = complex<double>(cos(dPSI[i]), sin(dPSI[i]));
    endPSI[i] = complex<double>(cos(dPSI[i]), -sin(dPSI[i]));
  }
  
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
  delete epdPSI;
  delete endPSI;
  
  return;
}

void calcEva12RecIMP(int depth)
{
	if(depth == Nc-1)
	{
			int nn = offset;	
			double aux, nkdPSI, phi;
            double Jn, Jnp1, Jnn1, u, theta, ProdJn;

    		complex <double> DeltaJ(0,0);
            complex <double> SDeltaJ(0,0);
            const complex<double> jj (0,1);
            const complex <double> jj05(0,0.5);
			complex <double> jjn(1,0);
            const complex <double> umc (1,0);
            
            nk[depth] = Nmax-1;
      
			jjn = pow(jj, sum(nk, ++depth)); // Atention to ++detph
			
			for(int i = 0; i < depth; i++)
				nn += ((int)Fcoeff[i])*nk[i];
            	
			for(int k = 2*Nmax-1; k >= 0; k--)
			{
                if(k == Nmax || nn == offset)
                {
                    jjn *= -jj;
                    nn -= (int)Fcoeff[depth-1];
                    continue;
                }

                nk[depth-1] = k-Nmax;

				// theta e exp(1j*nk.*dPSI)
                aux = 0;
                nkdPSI = 0;
                for(int i = depth-1; i >= 0; i--)
                {
                    aux += W[i]*nk[i];
                    nkdPSI += dPSI[i]*nk[i];
                }
                   
                aux *= Dhat;
                
                ProdJn = 1.0;
                SDeltaJ = complex<double>(0,0);
                phi = 0;
                
				for(int i = depth-1; i >= 0; i--)
				{
					theta = aux*W[i];
					u = 2*mFM[i]*sin(theta);
					theta = cos(theta);
                    
                    // J(nk, u), J(nk+1,u) e J(nk-1,u)
                    Jn = _jn(nk[i],  u);
                    Jnn1 = _jn(nk[i]-1, u); 
                    Jnp1 = _jn(nk[i]+1, u);
                    
                    DeltaJ = endPSI[i]*Jnn1;
                    DeltaJ -= epdPSI[i]*Jnp1;
                    
                    DeltaJ *= mIM[i]*theta/Jn;
                    
                    SDeltaJ += DeltaJ;
                    ProdJn *= Jn;
                    
                    phi += nk[i]*phiIM[i];
				}
				
				Idet[nn] += ProdJn*jjn*(umc - jj05*SDeltaJ)*(complex <double> (cos(phi + nkdPSI), sin(phi + nkdPSI)));

				jjn *= -jj;
				nn -= (int)Fcoeff[depth-1];
			}
	}
	else
	{
		for(int k = 2*Nmax-1; k >= 0; k--)
		{
            if(k == Nmax)
                continue;
            else
            {
                nk[depth] = k-Nmax;
                calcEva12RecIMP(depth+1);
            }
		}
	}
	return;
}

inline int sum(int * V, int length)
{
	int S = 0;
	for(int i = --length; i >= 0; i--)
		S += V[i];
	return S;
}

inline double sum(double * V, int length)
{
	double S = 0;
	for(int i = --length; i >= 0; i--)
		S += V[i];
	return S;
}
