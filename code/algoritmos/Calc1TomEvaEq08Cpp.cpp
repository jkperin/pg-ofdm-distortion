/** Calcula corrente Idet usando o modelo do artigo da Eva Peral equações 7 e 8 **/

#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793

using namespace std;

inline void calcEvaEq08(complex <double> * Idet, double P0, double mIM, double mFM, double dPSI, double theta, int N2);

//Idet = Calc1TomEvaEq08Cpp(mIM, mFM, dPSI, theta, P0, N2); 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int N2;
  double * IdetRe, * IdetIm;
  double mIM, mFM, dPSI, theta, P0;
  
  /* Check for proper number of arguments */
  if (nrhs != 6) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "São necessários 6 argumentos de entrada");
  } else if (nlhs != 1) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "Esta função retorna apenas um argumento");
  }
   
  mIM = mxGetScalar(prhs[0]);
  mFM = mxGetScalar(prhs[1]);
  dPSI = mxGetScalar(prhs[2]);
  theta = mxGetScalar(prhs[3]);
  P0 = mxGetScalar(prhs[4]);
  N2 = (int)mxGetScalar(prhs[5]);
    
  if(mIM < 0 || N2 < 0 || P0 < 0)
      mexErrMsgTxt("Erro!!! Argumentos de entrada inválidos");
    
  
  complex <double> * Idet = new complex <double> [2*N2];
  
  memset((void *)Idet, 0, 2*N2*sizeof(complex<double>));
  
  // Chama função recursiva
  calcEvaEq08(Idet, P0, mIM, mFM, dPSI, theta, N2);
  
  // Passa resultado para o MatLab
  plhs[0] = mxCreateDoubleMatrix(1, 2*N2, mxCOMPLEX);
  IdetRe = mxGetPr(plhs[0]);
  IdetIm = mxGetPi(plhs[0]);
  
  for(int i = 2*N2-1; i >= 0; i--)
  {
      IdetRe[i] = Idet[i].real();
      IdetIm[i] = Idet[i].imag();
  }
  
  // Limpeza
  delete Idet;
 
  return;
}

inline void calcEvaEq08(complex <double> * Idet, double P0, double mIM, double mFM, double dPSI, double theta, int N2)
{
    double u = 0.0;
    const complex <double> jj(0,1);
    complex <double> jjn(cos(-N2*pi/2),sin(-N2*pi/2));
      
    const complex <double> edPSIp (cos(dPSI),sin(dPSI));
    const complex <double> edPSIn (cos(dPSI),-sin(dPSI));
    complex <double> nedPSI (cos(-N2*dPSI),sin(-N2*dPSI));
    const complex <double> umc(1,0);
    
    for(int n = -N2, nn = 0; n < N2; n++)
    {
        u = 2*mFM*sin(n*theta);
        Idet[nn] = P0*jjn*nedPSI*(_jn(n,u)*umc - jj*(mIM/2*cos(n*theta))*(_jn(n-1,u)*edPSIn - _jn(n+1,u)*edPSIp));
                
        nn++;
        jjn *= jj;
        nedPSI *= edPSIp;
    }
    return;
}

// for n = -N2/2:N2/2-1
//     u = 2*mFM*sin(n*theta);
//     eva08.Idet(nn) = P0*(1j)^n*exp(1j*n*dPSI)*(J(n, u) - 1j*mIM/2*cos(n*theta)*(J(n-1,u)*exp(-1j*dPSI) - J(n+1,u)*exp(1j*dPSI)));
//     eva08.I = eva08.I + eva08.Idet(nn)*exp(1j*n*(w*t + phiIM));
//     nn = nn + 1;
// end
        
