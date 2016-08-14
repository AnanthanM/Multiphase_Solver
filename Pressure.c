//Every function for solving Pressure Poisson 
//Equation is in this file

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

void Get_RHS_of_Pressure_Poisson_Eqn(Domain domain,Field * RHS,Constant constant)
{
  Field * u   = domain.u;
  Field * v   = domain.v;
  Field * p   = domain.p;
   
  double dx = constant.dx;
  double dy = constant.dy;
  double dt = constant.dt;
  
  int Nx_v      = v->N_x;
  int Nx_u      = u->N_x; 
  int N_cells_C = p->N;
  int Nx_C      = p->N_x;
  int Ny_C      = p->N_y;

/********************************/
//           --- vN-----
//          |     |    |
//          uW---pij---uE
//          |     |    |
//          -----vS-----
/********************************/
  int l,i,j;
  int row;
  double  * val_uE, * val_uW, * val_vN, * val_vS;
  double uE,uW,vN,vS;

  for(l = 0;l<N_cells_C;l++)
  {
    RHS->val[l] = 0.0;
  }

  for(j=1;j<(Ny_C-1);j++)
  {
    row = j*Nx_C;              // interms of colocated grid

    val_uE = &u->val[row-j];
    val_uW = &u->val[row-j-1];
    val_vN = &v->val[row];       
    val_vS = &v->val[row-Nx_v];

    for(i=1;i<(Nx_C-1);i++)
    {
      l = j*Nx_C + i;

      uE = val_uE[i];
      uW = val_uW[i];
      vN = val_vN[i];
      vS = val_vS[i];

      RHS->val[l] = ( 1/(2*dt) ) * ( (uE - uW)/dx + (vN - vS)/dy ) ;

    }
  }
  
  return;
}

void Pressure_Poisson(Field * p, Constant constant,double * Ap_vector,Domain domain)
{
  Field * rho = domain.rho;

  double dx = constant.dx;
  double dy = constant.dy;

  int N_cells_C    = rho->N;
  int Nx_C         = rho->N_x;
  int Ny_C         = rho->N_y;

  int l,i,j;
  
  double   pP,
           pN,pS,pE,pW;                  
  
  double   rhoP,
           rhoN,rhoS,rhoE,rhoW;                  

  double   aP,
           aN,aS,aE,aW;    

  int      row;

  double   * val_pP, * val_pE, * val_pW, * val_pN, * val_pS;
  double   * val_rhoP, * val_rhoE, * val_rhoW, * val_rhoN, * val_rhoS;


  for(l=0;l<N_cells_C;l++)
    Ap_vector[l] = 0.0;

  //BC Setup
  
  //West Side Neumann 
  
  i = 0;
  for(j=0;j<Ny_C;j++)
  {
    l = j*Nx_C + i;
    p->val[l] = p->val[l+1];
  }

  //East Side Neumann
  
  i = Nx_C - 1;
  for(j=0;j<Ny_C;j++)
  {
    l = j*Nx_C + i;
    p->val[l] = p->val[l-1];
  }
  
  //South side Neumann

  j = 0;
  for(i=0;i<Nx_C;i++)
  {
    l = j*Nx_C + i;
    p->val[l] = p->val[l+Nx_C];
  }
  
  //North side Neumann

  j = Ny_C - 1;
  for(i=0;i<Nx_C;i++)
  {
    l = j*Nx_C + i;
    p->val[l] = p->val[l-Nx_C];
  }
  
  // Loops for all the inner cells 
  
  for(j=1;j<(Ny_C-1);j++)
  {
    row = j*Nx_C;

    val_pP = &p->val[row];
    val_pE = &p->val[row+1];
    val_pW = &p->val[row-1];
    val_pN = &p->val[row+Nx_C];
    val_pS = &p->val[row-Nx_C];

    val_rhoP = &rho->val[row];
    val_rhoE = &rho->val[row+1];
    val_rhoW = &rho->val[row-1];
    val_rhoN = &rho->val[row+Nx_C];
    val_rhoS = &rho->val[row-Nx_C];
    
    for(i=1;i<(Nx_C-1);i++)
    {
      l = j*Nx_C + i;

      pP = val_pP[i];
      pE = val_pE[i];
      pW = val_pW[i];
      pN = val_pN[i];
      pS = val_pS[i];

      rhoP = val_rhoP[i];
      rhoE = val_rhoE[i];
      rhoW = val_rhoW[i];
      rhoN = val_rhoN[i];
      rhoS = val_rhoS[i]; 
      
      aE = ( (1/(dx*dx))*(1/(rhoE+rhoP)) );  
      aW = ( (1/(dx*dx))*(1/(rhoP+rhoW)) );  
      aN = ( (1/(dy*dy))*(1/(rhoN+rhoP)) );  
      aS = ( (1/(dy*dy))*(1/(rhoP+rhoS)) );
      aP = aE + aW + aN + aS ;

      Ap_vector[l] = -pP*aP + pE*aE + pW*aW + pN*aN + pS*aS; 
    }
  }
  return;

}  
