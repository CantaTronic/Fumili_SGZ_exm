//-*-mode: C++; c-basic-offset: 1; -*-
 #include <TFile.h>
 #include <TTree.h>
 #include <GEvent.hh>
 #include <GTrack.hh>
 #include <TChain.h>
 #include "run_conditions.h"
 #include "declarations.h"
 #include <RKFitter.h>
 #include "TF_fwd.h"
 #include "GAnkeGeom.hh"
 #include <TLorentzVector.h>
 #include "GClusterListDC.hh"
 #include "TH1.h"
 #include "TH1.h"
 #include "TH2.h"
 #include "WireChamber.h"
 #include "SinglePlane.h"
 #include <string>
 #include <vector>
 #include <DriftWireCalib.hh>
 #include <utility>
 #include"ExpSet.hh"
 #include"FdSetup.hh"
 #include <stdio.h>
 #include <float.h>
  double data[100000][2];
 using namespace std;
#include "./Referencies"
int idebug;
int SGZ_Vika(int m, double &S, double A[], double PL[],
      double G[], double Z[])
{
 S = 0.;
 for (int i = 0; i < m; i++) {
  for (int j = 0; j <= i; j++) Z[(i+1)*i/2 + j ] = .0;
  G[i] = .0;
 }
 if(idebug)cout << " SGZ_Vika, m                   = " << m << endl;
 if(idebug)cout << " SGZ_Vika, Parameters at entry = " << endl;
 for(int i = 0; i < m;i++)
  {
   if(idebug)cout << " SGZ_Vika, i,A[i],PL[i]  = " << i << "  " << A[i] << "  " << PL[i] << endl;
  }
 for (int  iev = 0; iev < 100000; iev++)
  {
   double CH = (1.+ A[0]*data[iev][0] +A[1]*data[iev][1]);
   double ZN = (1. + 0.5*A[0] +0.5*A[1]);
   double pdf = CH/ZN;//pdf
   // Derivatives
   double df[2];
   for(int i = 0; i < 2;i++)
    {
     double dch = data[iev][i];
     double dzn = 0.5;
     double ys = (dch - pdf*dzn)/ZN;
     df[i] = ys/pdf;
    }
   S -= log(pdf);
   int il = 0;
   {
    for (int  ip = 0; ip < m; ip++)
     {
      // Calculation of gradient
      if (PL[ip] > .0)
       {
	G[ip] = G[ip] - df[ip];
	//if(idebug)cout << " G[ip],Z[il]: " << G[ip] << "  " << endl;
	{
	 for (int jp = 0; jp <= ip; jp++)
	  {
	   if (PL[jp] > .0)
	    {
	     Z[il] = Z[il] + df[ip]*df[jp];
	     //if(idebug)cout << Z[il]  << "  " << endl;
	     il = il + 1;
	    }
	  }
	}
       }
     }
   }
  }
 return 2;
}
void VikaC(int M,double A[],
				     double *psis,double **dpsis)
{
   /* Conservation laws for pp->pp. First proton is in Fd, second in STT */
   /* First - cleaning derivatives */
   //  if(idebug)cout << " get_psis_and_derivatives_for_All_pp_to_pp,M  : " << M << endl;
   //cout << " MissingMassSquared1Track,M  : " << M << endl;
   for(int j = 0;j < M;j++)
      {
	 if(idebug)cout << " j,A[j] " <<j << "  " << A[j] << endl; ; 
      }
   for(int i = 0;i < 1;i++)
      {
	 for(int j = 0;j < M;j++)
	    {
	       dpsis[i][j] = .0; 
	    }
      }
   psis[0] = A[0] + A[1] - 1.3;
   dpsis[0][0] = 1.;
   dpsis[0][1] = 1.;
   for(int i = 0; i < 1; i++)
    {
     if(idebug)cout << " i,psis[i],derivatives.. = " << i  << "  " << psis[i] << "  " ; 
     for(int j = 0; j < M; j++)
      {
       if(idebug)cout << dpsis[i][j] << "  " ; 
      }
     if(idebug)cout << endl;
     if(idebug)cout << "*******************************" << endl;
    }
}

   //idebug = 0;
/*

  Comparison of two fits - constrained and non constrained

 UnConstrained Fit: Fitted Values and their Errors 
0  0.514879 +/-  0.0229163
1  0.815076 +/-  0.0259613

 Constrained Fit: Fitted Values and their Errors 
0  0.501023 +/-  0.0126582
1  0.798977 +/-  0.0126582
 */
int main(int arc, char **argv)
 {
  psis_Out = new double[1];
  double *psis = psis_Out;
  ifstream ifs("unif.dat");
  int count = 0;
  for(int i= 0;i < 100000;i++)
   {
    ifs >> data[i][0] >> data[i][1] ;
    /*
    if(count < 5)
     {
      cout << data[i][0] << "  " << data[i][1] << endl;
     }*/
    count++;
   }
  double AMN[2],AMX[2],PL0[2],A[2],R[2],SIGMA[2],VL[2*2],akappa,
   S;
  for(int i = 0; i < 2;i++)
   {
    AMN[i] = -10.;
    AMX[i] = 10.;
    PL0[i] = .1;
    A[i] = .0;
   }
  
  int N1 = 1, N2 =  1 , N3 = 30, IT = -30;
  double EPS = .1;
  //   if(idebug)cout << "FitDataVDY.fnPar = " << FitDataVDY.fnPar << endl;
  //cout << "FitDataVDY.fnPar = " << FitDataVDY.fnPar << endl;
   idebug = 0;
   for(int i = 0; i < 2;i++)
      {
	 if(idebug)cout << " i,A,PL0,AMX,AMN "
	      << i << "  "
	      << A[i] << "  "
	      << PL0[i] << "  "
	      << AMX[i] << "  "
	      << AMN[i] << endl;
      }
   idebug = 0;
   int ii;
   ii = fumiliSK(2,  S, N1,  N2,  N3,
		 EPS, IT, A, PL0,
		 AMX, AMN, R, SIGMA,
		 SGZ_Vika, akappa, VL,
		 1,VikaC);
   if(akappa < EPS)
    {
     cout << " Constrained Fit: Fitted Values and their Errors " << endl;
     for(int i = 0;i < 2;i++)
      {
       cout << i << "  " << A[i] << " +/-  " << SIGMA[i] << endl;
      }
    }
  
return 1;
}

       
