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
  int N,S,W,E;
  double uE,uW,vN,vS;

  for(l = 0;l<N_cells_C;l++)
  {
    if(p->bc_type[l] == NONE)
    {
        i = l%Nx_C;
        j = (int) l/Nx_C;

        N = l;
        S = l - Nx_v;
        E = l - j;
        W = E - 1;

        uE = u->val[E];
        uW = u->val[W];
        vN = v->val[N];
        vS = v->val[S];

        RHS->val[l] = ( 1/(2*dt) ) * ( (uE - uW)/dx + (vN - vS)/dy ) ;
    }
    else
    {
        RHS->val[l] = 0.0;
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

  int      N,S,E,W;                       //for rho values
  
  double   pP,
           pN,pS,pE,pW;                  
  
  double   rhoP,
           rhoN,rhoS,rhoE,rhoW;                  

  double   aP,
           aN,aS,aE,aW;                  
  
  for(l = 0;l<N_cells_C;l++)
  {
    if(p->bc_type[l] == NONE)
    {
        i = l%Nx_C;
        j = (int) l/Nx_C;
  
        S = (j-1)*Nx_C + i;
        N = (j+1)*Nx_C + i;
        W = j*Nx_C + (i-1);      // or l-1
        E = j*Nx_C + (i+1);      // or l+1
        
        rhoP =rho->val[l];
        rhoE =rho->val[E];
        rhoW =rho->val[W];
        rhoN =rho->val[N];
        rhoS =rho->val[S];
        
        aE = ( (1/(dx*dx))*(1/(rhoE+rhoP)) );  
        aW = ( (1/(dx*dx))*(1/(rhoP+rhoW)) );  
        aN = ( (1/(dy*dy))*(1/(rhoN+rhoP)) );  
        aS = ( (1/(dy*dy))*(1/(rhoP+rhoS)) );
        aP = aE + aW + aN + aS ;
 
        pP =p->val[l];
        pE =p->val[E];
        pW =p->val[W];
        pN =p->val[N];
        pS =p->val[S];

        if(p->bc_type[N] != NONE)
          pN = pP;

        if(p->bc_type[S] != NONE)
          pS = pP;

        if(p->bc_type[E] != NONE)
          pE = pP;

        if(p->bc_type[W] != NONE)
          pW = pP;
        
        Ap_vector[l] = -pP*aP + pE*aE + pW*aW + pN*aN + pS*aS; 

    }
    else
    {
        Ap_vector[l] = 0.0;
    }
    
  }

  return;

}  
