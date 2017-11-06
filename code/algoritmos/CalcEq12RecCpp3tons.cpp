#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793
#define Nc 3

using namespace std;

void calcEva12(complex <double> *Idet);
inline double sum(double * V, int length);

double * mIM, * mFM, * dPSI, * W, * Fcoeff, * phiIM;
double Dhat, P0;
int Nmax, offset;
complex <double> * epdPSI, * endPSI;


//eva12.IdetREC = calcEva17Rec(Nf, mIM, mFM, dPSI, W, Fcoeff, beta2, L, P0, Nmax); 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t M, N;
  int Nf;
  double * IdetRe, * IdetIm;
  complex <double> * Idet;
  double beta2, L;
  
  /* Check for proper number of arguments */
  if (nrhs != 11) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "São necessários 11 argumentos de entrada");
  } else if (nlhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "Esta função retorna apenas um argumento");
  }
  
  M = mxGetM(prhs[1]);
  N = mxGetN(prhs[1]);
  
  Nf = (int)mxGetScalar(prhs[0]);
  beta2 = mxGetScalar(prhs[6]);
  L = mxGetScalar(prhs[7]);
  P0 = mxGetScalar(prhs[8]);
  Nmax = (int)mxGetScalar(prhs[9]);
  
  if(M != 1 || M != mxGetM(prhs[2]) || M != mxGetM(prhs[3]) || M != mxGetM(prhs[4]) ||
          M != mxGetM(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de linhas diferente de 1 ou vetores de entrada com tamanhos diferentes!");
  else if(N != Nc || N != mxGetN(prhs[2]) || N != mxGetN(prhs[3]) || N != mxGetN(prhs[4]) ||
          N != mxGetN(prhs[5]))
    mexErrMsgTxt("Erro!!! Número de colunas diferente de Nc ou Vetores de entrada com tamanhos diferentes!");
    
  if(L < 0 || P0 < 0 || Nmax < 0)
      mexErrMsgTxt("Erro!!! Argumentos de entrada inválidos");
  
  mIM = (double *) mxGetPr(prhs[1]);
  mFM = (double *) mxGetPr(prhs[2]);
  dPSI = (double *) mxGetPr(prhs[3]);
  W = (double *) mxGetPr(prhs[4]);
  Fcoeff = (double *) mxGetPr(prhs[5]);
  phiIM = (double *) mxGetPr(prhs[10]);
    
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
  
  Idet = new complex <double> [Nf];
  
  memset((void *)Idet, 0, Nf*sizeof(complex<double>));
  
  // Chama função recursiva
  calcEva12(Idet);
  
  // Passa resultado para o MatLab
  plhs[0] = mxCreateDoubleMatrix(1, Nf, mxCOMPLEX);
  IdetRe = mxGetPr(plhs[0]);
  IdetIm = mxGetPi(plhs[0]);
  
  for(int i = 0; i < Nf; i++)
  {
      IdetRe[i] = Idet[i].real();
      IdetIm[i] = Idet[i].imag();
  }
  
  // Limpeza
  delete Idet;
  delete epdPSI;
  delete endPSI;
  
  return;
}

void calcEva12(complex <double> *Idet)
{
    double Jn, Jnp1, Jnn1, u, theta, ProdJn, Sw, SdPSI, phi;
    int nn = 0;
    
    complex <double> nkdPSI(0,0);
    complex <double> DeltaJ(0,0);
    complex <double> SDeltaJ(0,0);
    const complex <double> jj05(0,0.5);
    complex <double> jjn(1,0);
    const complex <double> umc (1,0);
    
    int nk[] = {-Nmax,-Nmax,-Nmax};
    
    for(int * n1 = nk; *n1 < Nmax; *n1++)
        for(int * n2 = nk+1; *n2 < Nmax; *n2++)
            for(int * n3 = nk+2; *n3 < Nmax; *n3++)
            {
                jjn = complex <double> (cos((*n1+*n2+*n3)*pi/2), sin((*n1+*n2+*n3)*pi/2));
                nn = *n1*((int)Fcoeff[0]) + *n2*((int)Fcoeff[1]) + *n3*((int)Fcoeff[2]) + offset;
                
                Sw = Dhat*(*n1*W[0] + *n2*W[1] + *n3*W[2]);
                
                SdPSI = *n1*dPSI[0] + *n2*dPSI[1] + *n3*dPSI[2];
                nkdPSI = complex<double>(cos(SdPSI), sin(SdPSI));
                phi = 0;              
                for(int i = 2; i >= 0; i--)
				{
					theta = Sw*W[i];
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
                
                Idet[nn] += (P0*ProdJn)*jjn*nkdPSI*(umc - jj05*SDeltaJ)*(complex <double> (cos(phi), sin(phi)));
            }
    return;
}

//     teoricoc3b.Idet = zeros(size(fff));
//     mIM = mIM(1);
//     for n1 = -Nmax:Nmax-1
//         for n2 = -Nmax:Nmax-1
//             for n3 = -Nmax:Nmax-1
//                 nk = [n1 n2 n3];
//                 nn = sum(nk.*F/min(F)) + sum(N2/2*F/min(F)) + 1;
// //                 nn = sum(nk.*Find/dff) + sum(Nmax*Find/dff) + 1;
//                 theta = -(1/2)*beta2*L*Wind*sum(nk.*Wind);
//                 u = 2*mFMind.*sin(theta);
//                 teoricoc3b.Idet(nn) = teoricoc3b.Idet(nn) + P0*1j^(sum(nk))*exp(1j*sum(nk.*dPSIind))*prod(J(nk,u))*...
//                     (1 - sum(1j*mIM/2.*cos(theta).*(J(nk-1, u).*exp(-1j*dPSIind) - J(nk+1, u).*exp(1j*dPSIind))./J(nk,u)));
//             end
//         end
//     end

inline double sum(double * V, int length)
{
	double S = 0;
	for(int i = --length; i >= 0; i--)
		S += V[i];
	return S;
}
                
