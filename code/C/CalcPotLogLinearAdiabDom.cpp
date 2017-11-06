/* 	5. CalcPotLogLinearAdiabDom:
		a. Função para calcular componentes de intermodulação que serão utilizados  no cálculo da potência do sinal detectado.
		b. As IMPs calculadas são somadas em potência. Portanto, o resultado final é a potência média para todas as realizações possíveis.
		c. Todos os IMPs são computados, quer caiam na banda do sinal ou não.
		d. Função calcula até componente de 5ª ordem, sendo o mínimo de 3ª ordem. 
		e. Esta função não pode ser utilizada para o cálculo do SIR, pois não há compatibilidade com o modelo numérico.
        f. O cálculo é feito considerando-se chirp adiabático como dominante (chirp transiente é nulo).
*/

#include <math.h>
#include <complex>
#include "mex.h"

#define pi 3.141592653589793
#define Nmax 1
#define N3 5
#define N4 9
#define N5 21

#define ORD4
#define ORD5

using namespace std;

void calcEva12Var(double * Idet2, double * Idet3, double * Idet4, double * Idet5);
double calcvar(const int * nk, int nn, double w);
void J(double &Jn, double &Jnn1, double &Jnp1, int n, double x);
double J(int n, double x);
double dJ(int n, double x);

double * F, * W, * mIM, * mFM, * dPSI, * Idet2, * Idet3, * Idet4, * Idet5;
double Dhat;
int Nc, offset, * Fi;
complex <double> * epdPSI, * endPSI;
const complex <double> umc (1,0);
const complex <double> jj (0,1);

const int cc3[5][3] = {{-1, -1, +1},
                        {-1, +1, +1},
                        {+1, -1, +1},
                        {+1, +1, -1},
                        {+1, +1, +1}};
                    
const int cc4[9][4] = {{-1, -1, -1, +1},
                     {-1, -1, +1, +1},
                     {-1, +1, -1, +1},
                     {-1, +1, +1, +1},
                     {+1, -1, -1, +1},
                     {+1, -1, +1, +1},
                     {+1, +1, -1, +1},
                     {+1, +1, +1, -1},
                     {+1, +1, +1, +1}};

const int cc5[21][5] = {{-1, -1, -1, -1, +1},{-1, -1, -1, +1, +1},{-1, -1, +1, -1, +1},
                       {-1, -1, +1, +1, -1},{-1, -1, +1, +1, +1},{-1, +1, -1, -1, +1},
                       {-1, +1, -1, +1, +1},{-1, +1, +1, -1, +1},{-1, +1, +1, +1, -1},
                       {-1, +1, +1, +1, +1},{+1, -1, -1, -1, +1},{+1, -1, -1, +1, +1},
                       {+1, -1, +1, -1, +1},{+1, -1, +1, +1, -1},{+1, -1, +1, +1, +1},
                       {+1, +1, -1, -1, +1},{+1, +1, -1, +1, +1},{+1, +1, +1, -1, -1},
                       {+1, +1, +1, -1, +1},{+1, +1, +1, +1, -1},{+1, +1, +1, +1, +1}};                     

// F, Fi, mIM, mFM, dPSI, beta2L, P0, N 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t M, N;
  int Nf;
  double * Fid, * Idet2Re, * Idet3Re, * Idet4Re, * Idet5Re;
  double beta2L, P0;
  
  /* Check for proper number of arguments */
  if (nrhs != 8) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "São necessários 11 argumentos de entrada");
  } else if (nlhs != 4) {
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
            "Esta função retorna apenas um argumento");
  }
  
  M = mxGetM(prhs[0]);
  N = mxGetN(prhs[0]);
  Nc = (int)N;
  
  beta2L = mxGetScalar(prhs[5]);
  P0 = mxGetScalar(prhs[6]);
  Nf = (int)mxGetScalar(prhs[7]);
  
  if(M != 1 || M != mxGetM(prhs[1]) || M != mxGetM(prhs[2]) || M != mxGetM(prhs[3]) || M != mxGetM(prhs[4]))
    mexErrMsgTxt("Erro!!! Número de linhas diferente de 1 ou vetores de entrada com tamanhos diferentes!");
  else if(N < 2 || N != mxGetN(prhs[1]) || N != mxGetN(prhs[2]) || N != mxGetN(prhs[3]) || N != mxGetN(prhs[4]))
    mexErrMsgTxt("Erro!!! Número de colunas diferente de Nc ou Vetores de entrada com tamanhos diferentes!");
  else if(P0 < 0 || Nf < 0)
      mexErrMsgTxt("Erro!!! Argumentos de entrada inválidos");
  
  F = (double *) mxGetPr(prhs[0]);
  Fid = (double *) mxGetPr(prhs[1]);
  mIM = (double *) mxGetPr(prhs[2]);
  mFM = (double *) mxGetPr(prhs[3]);
  dPSI = (double *) mxGetPr(prhs[4]);
  
  // Realiza alguns cálculos antes de chamar a função recursiva

  epdPSI = new complex <double> [Nc];
  endPSI = new complex <double> [Nc];
  W = new double [Nc];
  Fi = new int [Nc];
  for(int i = 0; i < Nc; i++)
  {
    epdPSI[i] = complex<double>(cos(dPSI[i]), sin(dPSI[i]));
    endPSI[i] = complex<double>(cos(dPSI[i]), -sin(dPSI[i]));
    W[i] = 2*pi*F[i];
    Fi[i] = (int)Fid[i];
  }
  
  Dhat = -0.5*beta2L;
  
  Idet2 = new double [Nf/2];
  Idet3 = new double [Nf/2];
  Idet4 = new double [Nf/2];
  Idet5 = new double [Nf/2];
 
  memset((void *)Idet2, 0, Nf/2*sizeof(double));
  memset((void *)Idet3, 0, Nf/2*sizeof(double));
  memset((void *)Idet4, 0, Nf/2*sizeof(double));
  memset((void *)Idet5, 0, Nf/2*sizeof(double));
  
  #ifdef ORD5
    mexPrintf("Máxima ordem de IMPs = 5\n");
  #elif ORD4
    mexPrintf("Máxima ordem de IMPs = 4\n");
  #else
    mexPrintf("Máxima ordem de IMPs = 3\n");
  #endif
  
  // Chama função recursiva
  calcEva12Var(Idet2, Idet3, Idet4, Idet5);
  // Passa resultado para o MatLab
  plhs[0] = mxCreateDoubleMatrix(1, Nf, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, Nf, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, Nf, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1, Nf, mxREAL);
  
  Idet2Re = mxGetPr(plhs[0]);
  Idet3Re = mxGetPr(plhs[1]);
  Idet4Re = mxGetPr(plhs[2]);
  Idet5Re = mxGetPr(plhs[3]);
  
  P0 *= P0;
  offset = Nf/2;
  for(int i = 1; i < Nf/2; i++)
  {
      Idet2Re[i] = P0*Idet2[Nf/2-i];
      Idet3Re[i] = P0*Idet3[Nf/2-i];
      Idet4Re[i] = P0*Idet4[Nf/2-i];
      Idet5Re[i] = P0*Idet5[Nf/2-i];
      
      Idet2Re[i+offset] = P0*Idet2[i];
      Idet3Re[i+offset] = P0*Idet3[i];
      Idet4Re[i+offset] = P0*Idet4[i];
      Idet5Re[i+offset] = P0*Idet5[i];    
  }
  
  // Limpeza
  delete Idet2;
  delete Idet3;
  delete Idet4;
  delete Idet5;
  delete epdPSI;
  delete endPSI;
  delete W;
  delete Fi;
  
  return;
}

void calcEva12Var(double * Idet2, double * Idet3, double * Idet4, double * Idet5)
{
    int nn;
    double fk, w, theta, ctheta;
    int * nk = new int [Nc];
    memset(nk, 0, Nc*sizeof(int));
       
    for(int k1 = 0; k1 < Nc; k1++)
    {
        for(int k2 = k1 + 1; k2 < Nc; k2++)
        {   
            fk = -F[k1] + F[k2];
            if(fk > 0)
            {
                w = 2*pi*fk;
                nn = -Fi[k1] + Fi[k2];
                nk[k1] = -1; nk[k2] = 1;
                Idet2[nn] += calcvar(nk, nn, w);
                nk[k1] = 0; nk[k2] = 0;
            }
            fk = F[k1] + F[k2];
            if(fk > 0)
            {
                w = 2*pi*fk;
                nn = Fi[k1] + Fi[k2];
                nk[k1] = 1; nk[k2] = 1;
                Idet2[nn] += calcvar(nk, nn, w);
                nk[k1] = 0; nk[k2] = 0;
            }
            
            for(int k3 = k2 + 1; k3 < Nc; k3++)
            {
                for(int i = 0; i < N3; i++)
                {
                    fk = cc3[i][0]*F[k1] + cc3[i][1]*F[k2] + cc3[i][2]*F[k3];
                    w = 2*pi*fk;
                    if(fk > 0)
                    {
                        nn = cc3[i][0]*Fi[k1] + cc3[i][1]*Fi[k2] + cc3[i][2]*Fi[k3];
                        nk[k1] = cc3[i][0]; nk[k2] = cc3[i][1]; nk[k3] = cc3[i][2];
                        Idet3[nn] += calcvar(nk, nn, w);
                        nk[k1] = 0; nk[k2] = 0; nk[k3] = 0;
                    }
                }
                #ifdef ORD4
                for(int k4 = k3 + 1; k4 < Nc; k4++)
                {
                    for(int i = 0; i < N4; i++)
                    {
                        fk = cc4[i][0]*F[k1] + cc4[i][1]*F[k2] + cc4[i][2]*F[k3] + cc4[i][3]*F[k4];
                        w = 2*pi*fk;
                        if(fk > 0)
                        {
                            nn = cc4[i][0]*Fi[k1] + cc4[i][1]*Fi[k2] + cc4[i][2]*Fi[k3] + cc4[i][3]*Fi[k4];
                            nk[k1] = cc4[i][0]; nk[k2] = cc4[i][1]; nk[k3] = cc4[i][2]; nk[k4] = cc4[i][3];
                            Idet4[nn] += calcvar(nk, nn, w);
                            nk[k1] = 0; nk[k2] = 0; nk[k3] = 0; nk[k4] = 0;
                        }
                    }
                    #ifdef ORD5
                    for(int k5 = k4 + 1; k5 < Nc; k5++)
                    {
                        for(int i = 0; i < N5; i++)
                        {
                            fk = cc5[i][0]*F[k1] + cc5[i][1]*F[k2] + cc5[i][2]*F[k3] + cc5[i][3]*F[k4] + cc5[i][4]*F[k5];
                            w = 2*pi*fk;
                            if(fk > 0)
                            {
                                nn = cc5[i][0]*Fi[k1] + cc5[i][1]*Fi[k2] + cc5[i][2]*Fi[k3] + cc5[i][3]*Fi[k4] + cc5[i][4]*Fi[k5];
                                nk[k1] = cc5[i][0]; nk[k2] = cc5[i][1]; nk[k3] = cc5[i][2]; nk[k4] = cc5[i][3]; nk[k5] = cc5[i][4];
                                Idet5[nn] += calcvar(nk, nn, w);
                                nk[k1] = 0; nk[k2] = 0; nk[k3] = 0; nk[k4] = 0; nk[k5] = 0;
                            }
                        }
                    }
                    #endif
                }
                #endif
            }
        }
    }
    delete nk;
    return;
}

double calcvar(const int *nk, int nn, double w)
{
    double Jn, dJn, u, theta, ctheta, ProdJn = 1.0;
    complex <double> SDeltaJ = complex<double>(0,0);
    for(int k = 0; k < Nc; k++)
    {
        theta = Dhat*W[k]*w;
        ctheta = cos(theta);
        u = 2*mFM[k]*sin(theta);

        Jn = J(nk[k], u);
        dJn = dJ(nk[k],u);

        SDeltaJ += (mIM[k]*ctheta*dJn/Jn);
        ProdJn *= Jn;
    }
    SDeltaJ = ProdJn*(umc - jj*SDeltaJ);    
    return (abs(SDeltaJ)*abs(SDeltaJ));  
}

void J(double &Jn, double &Jnn1, double &Jnp1, int n, double x)
{
    if(n == 0)
    {
        Jn = 1.0 - (x*x)/4.0;
        Jnp1 = x/2.0;
        Jnn1 = -Jnp1;
    }
    else if(n == 1)
    {
        Jn = x/2.0;
        Jnp1 = 0;
        Jnn1 = 1.0 - (x*x)/4.0;
    }
    else
    {
        Jn = -x/2.0;
        Jnp1 = 1.0 - (x*x)/4.0;
        Jnn1 = 0;
    }
    return;
}

double J(int n, double x)
{
    double Jn;
    if(n == 0)
        Jn = 1- (x*x)/4.0;
    else
        Jn = n*x/2.0;
    return Jn;
}

double dJ(int n, double x)
{
    double Jn;
    if(n == 0)
        Jn = -x/2.0;
    else
        Jn = n/2.0;
    return Jn;
}