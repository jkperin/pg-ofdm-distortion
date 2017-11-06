#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793

using namespace std;

void calcEva17Rec(int depth);
inline int sum(int * V, int length);
inline double sum(double * V, int length);

double * mIM, * mFM, * dPSI, * W, * Fcoeff;
double Dhat, P0;
int Nc, Nmax, * nk, nfmax;
complex <double> * epdPSI, * endPSI, * Idet;


//eva12.IdetREC = calcEva17Rec(Nf, mIM, mFM, dPSI, W, Fcoeff, beta2, L, P0, Nmax, fmax/dll); 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t M, N;
  int Nf, offset;
  double * IdetRe, * IdetIm;
  double beta2, L;
  
  /* Check for proper number of arguments */
  if (nrhs != 11) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "São necessários 10 argumentos de entrada");
  } else if (nlhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "Esta função retorna apenas um argumento");
  }
  
  M = mxGetM(prhs[1]);
  N = mxGetN(prhs[1]);
  
  Nc = (int)N;
  Nf = (int)mxGetScalar(prhs[0]);
  beta2 = mxGetScalar(prhs[6]);
  L = mxGetScalar(prhs[7]);
  P0 = mxGetScalar(prhs[8]);
  Nmax = (int)mxGetScalar(prhs[9]);
  nfmax = (int)mxGetScalar(prhs[10]) + 1;
  
  if(M != 1 || M != mxGetM(prhs[2]) || M != mxGetM(prhs[3]) || M != mxGetM(prhs[4]) ||
          M != mxGetM(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de linhas diferente de 1 ou vetores de entrada com tamanhos diferentes!");
  else if(N < 2 || N != mxGetN(prhs[2]) || N != mxGetN(prhs[3]) || N != mxGetN(prhs[4]) ||
          N != mxGetN(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de colunas diferente de Nc ou Vetores de entrada com tamanhos diferentes!");
    
  if(L < 0 || P0 < 0 || Nmax < 0)
      mexErrMsgTxt("Erro!!! Argumentos de entrada inválidos");
  
  mIM = (double *) mxGetPr(prhs[1]);
  mFM = (double *) mxGetPr(prhs[2]);
  dPSI = (double *) mxGetPr(prhs[3]);
  W = (double *) mxGetPr(prhs[4]);
  Fcoeff = (double *) mxGetPr(prhs[5]);
    
  // Realiza alguns cálculos antes de chamar a função recursiva
  offset = Nmax*((int)sum(Fcoeff, Nc));

  epdPSI = new complex <double> [Nc];
  endPSI = new complex <double> [Nc];
  for(int i = 0; i < Nc; i++)
  {
    epdPSI[i] = complex<double>(cos(dPSI[i]), sin(dPSI[i]));
    endPSI[i] = complex<double>(cos(dPSI[i]), -sin(dPSI[i]));
  }
  
  Dhat = -0.5*beta2*L;
  
  Idet = new complex <double> [nfmax+1];
  nk = new int[Nc];
  
  memset((void *)Idet, 0, (nfmax+1)*sizeof(complex<double>));
  memset((void *)nk, 0, sizeof(int));
  
  // Chama função recursiva
  calcEva17Rec(0);
  
  // Passa resultado para o MatLab
  plhs[0] = mxCreateDoubleMatrix(1, 2*nfmax, mxCOMPLEX);
  IdetRe = mxGetPr(plhs[0]);
  IdetIm = mxGetPi(plhs[0]);
  
  for(int i = nfmax, ii = 0; i < 2*nfmax; i++, ii++)
  {
      IdetRe[i] = Idet[ii].real();
      IdetIm[i] = Idet[ii].imag();
  }
  for(int i = nfmax-1, ii = 1; i >= 0; i--, ii++)
  {
      IdetRe[i] = Idet[ii].real();
      IdetIm[i] = -Idet[ii].imag();
  }
  
  // Limpeza
  delete Idet;
  delete nk;
  delete epdPSI;
  delete endPSI;
  
  return;
}

void calcEva17Rec(int depth)
{
	if(depth == Nc-1)
	{
			int nn = 0;	
			double aux, aux2;
            double Jn, Jnp1, Jnn1, u, theta, ProdJn;

			complex <double> nkdPSI(0,0);
			complex <double> DeltaJ(0,0);
            complex <double> SDeltaJ(0,0);
            const complex<double> jj (0,1);
            const complex <double> jj05(0,0.5);
			complex <double> jjn(1,0);
            const complex <double> umc (1,0);
            
            nk[depth] = Nmax-1;
      
			jjn = pow(jj, sum(nk, ++depth));
			
			for(int i = 0; i < depth; i++)
				nn += ((int)Fcoeff[i])*nk[i];
            	
			for(int k = 2*Nmax-1; k >= 0; k--)
			{
                if(nn <= 0 || nn > nfmax)
                {
                    jjn *= -jj;
                    nn -= (int)Fcoeff[depth-1];
                    continue;
                }
                
				nk[depth-1] = k-Nmax;

				// theta e exp(1j*nk.*dPSI)
                aux = 0;
                aux2 = 0;
                for(int i = depth-1; i >= 0; i--)
                {
                    aux += W[i]*nk[i];
                    aux2 += dPSI[i]*nk[i];
                }
                   
                aux *= Dhat;
                nkdPSI = complex<double>(cos(aux2), sin(aux2));
                
                ProdJn = 1.0;
                SDeltaJ = complex<double>(0,0);
                
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
				}
				
				Idet[nn] += (P0*ProdJn)*jjn*nkdPSI*(umc - jj05*SDeltaJ);

				jjn *= -jj;
				nn -= (int)Fcoeff[depth-1];
			}
	}
	else
	{
		for(int k = 2*Nmax-1; k >= 0; k--)
		{
			nk[depth] = k-Nmax;
			calcEva17Rec(depth+1);
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
