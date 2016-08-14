/*************** Diffusion *****************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

//
//Diffusion in u velocity
//the result is updated to advection part in u_T
//

void Diffusion_u(Domain domain, Constant constant, Field * u_T)
{
  Field * u   = domain.u;
  Field * v   = domain.v;
  Field * rho = domain.rho;
  Field * mu  = domain.mu;

  double dx = constant.dx;
  double dy = constant.dy;

  int N_cells_u    = u->N;
  int Nx_u         = u->N_x;
  int Ny_u         = u->N_y;

  int N_cells_v    = v->N;
  int Nx_v         = v->N_x;
  int Ny_v         = v->N_y;

  int N_cells_C    = mu->N;
  int Nx_C         = mu->N_x;
  int Ny_C         = mu->N_y;

/********************************/  
//                uN
//                |
//         vNW-- muN--vNE
//          |     |    |
//      uW--muW--uij--muE--uE
//          |     |    |
//         vSW-- muS--vSE
//                |
//                uS
//
//  we need to interpolate and find 
//  muN and muS , everything else is readily available
//
//  we also need rhoE and rhoW in the same positions
//  as muE and muW
/********************************/  
  
  int l,i,j;

  int row;

  double   uP,
           uN,uS,uE,uW,                    
           vNW,vNE,vSE,vSW,                
           mu_NE,mu_NW,mu_E,mu_W,mu_SE,mu_SW, //actual mu values
           mu_N,mu_S,                         //interpolated mu values 
           rho_E,rho_W;
  
  double   * val_uP,
           * val_uN,* val_uS,* val_uE,* val_uW,                    
           * val_vNW,* val_vNE,* val_vSE,* val_vSW,                
           * val_mu_NE,* val_mu_NW,* val_mu_E,* val_mu_W,* val_mu_SE,* val_mu_SW, //actual mu values
           * val_rho_E,* val_rho_W;


  //BC Setup
  
  //West Side 
  
  i = 0;
  for(j=0;j<Ny_u;j++)
  {
    l = j*Nx_u + i;
    u->val[l] = u->BC_Value[XMIN];
  }

  //East Side 
  
  i = Nx_u - 1;
  for(j=0;j<Ny_u;j++)
  {
    l = j*Nx_u + i;
    u->val[l] = u->BC_Value[XMAX];
  }
  

  //South side 

  j = 0;
  for(i=0;i<Nx_u;i++)
  {
    l = j*Nx_u + i;
    u->val[l]   = (2*u->BC_Value[YMIN]) - u->val[l+Nx_u];
    v->val[l+1] = v->BC_Value[YMIN];
  }
  
  //North side 

  j = Ny_u - 1;
  for(i=0;i<Nx_u;i++)
  {
    l = j*Nx_u + i;
    u->val[l]   = (2*u->BC_Value[YMAX]) - u->val[l-Nx_u];
    v->val[l-1] = v->BC_Value[YMAX];
  }
  
  // Loops for all the inner cells 
  
  for(j=1;j<(Ny_u-1);j++)
  {
    row = j*Nx_u;

    val_uP = &u->val[row];
    val_uE = &u->val[row+1];
    val_uW = &u->val[row-1];
    val_uN = &u->val[row+Nx_u];
    val_uS = &u->val[row-Nx_u];

    val_vNW = &v->val[row+j];
    val_vNE = &v->val[row+j+1];
    val_vSW = &v->val[row+j-Nx_v];
    val_vSE = &v->val[row+j-Nx_v+1];

    val_mu_W  = &mu->val[row+j];
    val_mu_E  = &mu->val[row+j+1];
    val_mu_NW = &mu->val[row+j+Nx_C];
    val_mu_NE = &mu->val[row+j+Nx_C+1];
    val_mu_SW = &mu->val[row+j-Nx_C];
    val_mu_SE = &mu->val[row+j-Nx_C+1];
    
    val_rho_W  = &rho->val[row+j];
    val_rho_E  = &rho->val[row+j+1];
    
    for(i=1;i<(Nx_u-1);i++)
    {
      l = j*Nx_u + i;
     
      uP = val_uP[i];
      uE = val_uE[i];
      uW = val_uW[i];
      uN = val_uN[i];
      uS = val_uS[i];

      vNW = val_vNW[i];
      vNE = val_vNE[i];
      vSW = val_vSW[i];
      vSE = val_vSE[i];

      mu_W  = val_mu_W[i];  
      mu_E  = val_mu_E[i];
      mu_NW = val_mu_NW[i];
      mu_NE = val_mu_NE[i];
      mu_SW = val_mu_SW[i];
      mu_SE = val_mu_SE[i];
              
      rho_W = val_rho_W[i];
      rho_E = val_rho_E[i];

      mu_N  = 0.25*( mu_NE + mu_NW + mu_W  + mu_E  );
      mu_S  = 0.25*( mu_E  + mu_W  + mu_SW + mu_SE );                        //interpolated mu values 
      

      u_T->val[l] = u_T->val[l] + 
                      (
                      (2.0/(rho_E+rho_W))*
                      ( (1/dx)*( 2*mu_E*( (uE - uP)/dx ) - 
                                 2*mu_W*( (uP - uW)/dx ) )  +
                        (1/dy)*( mu_N*( ( (uN-uP)/dy ) + ( (vNE-vNW)/dx ) ) - 
                                 mu_S*( ( (uP-uS)/dy ) + ( (vSE-vSW)/dx ) ) ) )
                      );


    }
  }
  return; 

}


//
//Diffusion in v velocity
//the result is updated to advection part in v_T
//

void Diffusion_v(Domain domain, Constant constant, Field * v_T)
{
  Field * u   = domain.u;
  Field * v   = domain.v;
  Field * rho = domain.rho;
  Field * mu  = domain.mu;

  double dx = constant.dx;
  double dy = constant.dy;

  int N_cells_u    = u->N;
  int Nx_u         = u->N_x;
  int Ny_u         = u->N_y;

  int N_cells_v    = v->N;
  int Nx_v         = v->N_x;
  int Ny_v         = v->N_y;

  int N_cells_C    = mu->N;
  int Nx_C         = mu->N_x;
  int Ny_C         = mu->N_y;

/********************************/  
//                vN
//                |
//         uNW-- muN--uNE
//          |     |    |
//      vW--muW--vij--muE--vE
//          |     |    |
//         uSW-- muS--uSE
//                |
//                vS
//
//  we need to interpolate and find 
//  muE and muW , everything else is readily available
//
//  we also need rhoN and rhoS in the same positions
//  as muN and muS
/********************************/  
  
  int l,i,j;

  int row;

  double   vP,
           vN,vS,vE,vW,                    
           uNW,uNE,uSE,uSW,                
           mu_N,mu_NW,mu_NE,mu_S,mu_SE,mu_SW, //actual mu values
           mu_W,mu_E,                         //interpolated mu values 
           rho_N,rho_S;
 
  double   * val_vP,
           * val_vN,* val_vS,* val_vE,* val_vW,                    
           * val_uNW,* val_uNE,* val_uSE,* val_uSW,                
           * val_mu_N,* val_mu_NW,* val_mu_NE,* val_mu_S,* val_mu_SE,* val_mu_SW, //actual mu values
           * val_rho_N,* val_rho_S;

  //BC Setup
  
  //West Side 
  
  i = 0;
  for(j=0;j<Ny_v;j++)
  {
    l = j*Nx_v + i;
    v->val[l]   = (2*v->BC_Value[XMIN]) - v->val[l+1] ;
    u->val[l-j] = u->BC_Value[XMIN];
  }

  //East Side 
  
  i = Nx_v - 1;
  for(j=0;j<Ny_v;j++)
  {
    l = j*Nx_v + i;
    v->val[l]     = (2*v->BC_Value[XMAX]) - v->val[l-1] ;
    u->val[l-j-1] = u->BC_Value[XMAX];
  }
  

  //South side 

  j = 0;
  for(i=0;i<Nx_v;i++)
  {
    l = j*Nx_v + i;
    v->val[l] = v->BC_Value[YMIN];
  }
  
  //North side 

  j = Ny_v - 1;
  for(i=0;i<Nx_v;i++)
  {
    l = j*Nx_v + i;
    v->val[l] = v->BC_Value[YMAX];
  }
  
  // Loops for all the inner cells 
  
  for(j=1;j<(Ny_v-1);j++)
  {
    row = j*Nx_v;

    val_vP = &v->val[row];
    val_vE = &v->val[row+1];
    val_vW = &v->val[row-1];
    val_vN = &v->val[row+Nx_v];
    val_vS = &v->val[row-Nx_v];

    val_uNE = &u->val[row+Nx_u-j];
    val_uNW = &u->val[row+Nx_u-j-1];
    val_uSE = &u->val[row-j];
    val_uSW = &u->val[row-j-1];

    val_mu_S  = &mu->val[row];
    val_mu_SE = &mu->val[row+1];
    val_mu_SW = &mu->val[row-1];
    val_mu_N  = &mu->val[row+Nx_C];
    val_mu_NE = &mu->val[row+Nx_C+1];
    val_mu_NW = &mu->val[row+Nx_C-1];

    val_rho_S  = &rho->val[row];
    val_rho_N  = &rho->val[row+Nx_C];

    
    for(i=1;i<(Nx_v-1);i++)
    {
      l = j*Nx_v + i;
     
      vP = val_vP[i];
      vE = val_vE[i];
      vW = val_vW[i];
      vN = val_vN[i];
      vS = val_vS[i];

      uNE = val_uNE[i];
      uNW = val_uNW[i];
      uSE = val_uSE[i];
      uSW = val_uSW[i];

      mu_S  = val_mu_S[i]; 
      mu_SE = val_mu_SE[i];
      mu_SW = val_mu_SW[i];
      mu_N  = val_mu_N[i]; 
      mu_NE = val_mu_NE[i];
      mu_NW = val_mu_NW[i];

      rho_S  = val_rho_S[i]; 
      rho_N  = val_rho_N[i]; 

      mu_W  = 0.25*( mu_N  + mu_NW + mu_S  + mu_SW  );
      mu_E  = 0.25*( mu_N  + mu_NE + mu_SE + mu_S   );                        //interpolated mu values 


      v_T->val[l] = v_T->val[l] + 
                      (
                      (2.0/(rho_N+rho_S))*
                      ( (1/dx)*( mu_E*( ( (uNE-uSE)/dy ) + ( (vE-vP)/dx ) ) -
                                 mu_W*( ( (uNW-uSW)/dy ) + ( (vP-vW)/dx ) ) ) + 
                        (1/dy)*( 2*mu_N*( (vN - vP)/dy ) - 
                                 2*mu_S*( (vP - vS)/dy ) ) )
                      );


    }
  }

  return; 

}
