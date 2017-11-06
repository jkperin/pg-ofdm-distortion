/** Calcula corrente Idet usando o modelo do artigo da Eva Peral equações 7 e 8 **/

#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793
#define MAXIT 10

using namespace std;

inline void calcEvaEq17(complex <double> * Idet, double P0, double mIM, double mFM, double dPSI, double theta, int N2);
double M_l(int l, double mIM);
double fatorial(int N);
double dmCm_mCk(int m, int k);

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
  calcEvaEq17(Idet, P0, mIM, mFM, dPSI, theta, N2);
  
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

inline void calcEvaEq17(complex <double> * Idet, double P0, double mIM, double mFM, double dPSI, double theta, int N2)
{
    double u = 0.0;
    
    for(int n = -N2, nn = 0; n < N2; n++, nn++)
    {
        u = 2*mFM*sin(n*theta);
        for(int l1 = -MAXIT; l1 < MAXIT; l1++)
            for(int l2 = -MAXIT; l2 < MAXIT; l2++)
                Idet[nn] += (P0*M_l(l1,mIM)*M_l(l2,mIM)*_jn(n+l1-l2,u))*complex<double>(cos((n + l1 - l2)*(pi/2 + dPSI) + n*theta*(l1+l2)), 
                        sin((n + l1 - l2)*(pi/2 + dPSI) + n*theta*(l1+l2)));
    }
    
    return;
}

// for n = -N2/2:N2/2-1
//     u = 2*mFM*sin(n*theta);
//     for l1 = -MAXIT:MAXIT
//         for l2 = -MAXIT:MAXIT
//             ml1 = 0;
//             ml2 = 0;
//             for k = 0:MAXIT
//                 m = 2*k + abs(l1);
//                 ml1 = ml1 + 1/(1-2*m)*binomialc(2*m,m)*binomialc(m,k)*(-mIM/8)^m;
//                 m = 2*k + abs(l2);
//                 ml2 = ml2 + 1/(1-2*m)*binomialc(2*m,m)*binomialc(m,k)*(-mIM/8)^m;
//             end
//             eva17.Idet(nn) = eva17.Idet(nn) + P0*(1j)^(n + l1 - l2)*ml1*ml2*exp(1j*(n+l1-l2)*dPSI)*exp(1j*n*theta*(l1+l2))*J((n+l1-l2),u);
//         end
//     end
//     eva17.I = eva17.I + eva17.Idet(nn)*exp(1j*n*(w*t + phiIM));
//     nn = nn+1;
// end


double M_l(int l, double mIM)
{
    double ml = 0.0;
    l = abs(l);
    for(int k = 0, m = 0; k < 2*MAXIT; k++)
    {
        m = 2*k + l;
        ml += (1.0/(1.0-2.0*(double)m))*(dmCm_mCk(m,k))*pow(-mIM/8.0,m);
    }
    return ml;
}
                
double dmCm_mCk(int m, int k)
{
    double R = 1.0;
    for(int i = 2*m; i > m; i--)
        R *= i;
              
    return (R/(fatorial(k)*fatorial(m-k)));
}

double fatorial(int N)
{
    if (N < 0)
    {
        mexPrintf("Argumento inválido para fatorial");
        return 0.0;
    }
    else if (N == 0 || N == 1)
        return 1.0;
    else
        return ((double)N)*fatorial(N-1);
}
           
