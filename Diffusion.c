#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

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

  int      N,S,E,W,                       //for u velocities
           NW,NE,SE,SW,                   //for v velocities
           C_NE,C_NW,C_E,C_W,C_SE,C_SW;   //for colocated variables(mu & rho) 
  double   uP,
           uN,uS,uE,uW,                    
           vNW,vNE,vSE,vSW,                
           mu_NE,mu_NW,mu_E,mu_W,mu_SE,mu_SW, //actual mu values
           mu_N,mu_S,                         //interpolated mu values 
           rho_E,rho_W;
  

  for(l = 0;l<N_cells_u;l++)
  {
    if(u->bc_type[l] == NONE)
    {
        i = l%Nx_u;
        j = (int) l/Nx_u;
  
        S = (j-1)*Nx_u + i;
        N = (j+1)*Nx_u + i;
        W = j*Nx_u + (i-1);      // or l-1
        E = j*Nx_u + (i+1);      // or l+1
  
        NW =  l + j;
        NE = NW + 1;
        SW = NW - Nx_v;
        SE = SW + 1;

        C_W  = l + j; 
        C_NW = C_W + Nx_C;
        C_NE = C_NW + 1;
        C_E  = C_W + 1;
        C_SE = C_E - Nx_C;
        C_SW = C_SE - 1;     //for colocated variables(mu & rho) 
  
        uP = u->val[l];
        
        uE = u->val[E];
        uW = u->val[W];
        uN = u->val[N];
        uS = u->val[S];
  
        vNE = v->val[NE];
        vNW = v->val[NW];
        vSE = v->val[SE];
        vSW = v->val[SW];
        
        mu_NE = mu->val[C_NE];    
        mu_NW = mu->val[C_NW];    
        mu_E  = mu->val[C_E];  
        mu_W  = mu->val[C_W]; 
        mu_SE = mu->val[C_SE];  
        mu_SW = mu->val[C_SW];        //actual mu values
        
        mu_N  = 0.25*( mu_NE + mu_NW + mu_W  + mu_E  );
        mu_S  = 0.25*( mu_E  + mu_W  + mu_SW + mu_SE );                        //interpolated mu values 
        
        rho_E = rho->val[C_E];
        rho_W = rho->val[C_W];

        if( i == 1 || i == (Nx_u-2) || j == 1 || j == (Ny_u-2) ) // Next to boundary 
        {
            if(u->bc_type[N] != NONE)
            {
              uN  = (2* u->BC_Value[YMAX]) - uP;
              vNW = v->BC_Value[YMAX];
              vNE = v->BC_Value[YMAX];
            }
            if(u->bc_type[S] != NONE)
            {
              uS  = (2* u->BC_Value[YMIN]) - uP;
              vSW = v->BC_Value[YMIN];
              vSE = v->BC_Value[YMIN];
            }
            if(u->bc_type[W] != NONE)
            { 
  
              uW  = u->BC_Value[XMIN];
  
            }
            if(u->bc_type[E] != NONE)
            {
  
              uE  = u->BC_Value[XMAX];
  
            }
        }
        
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

  int      N,S,E,W,                       //for v velocities
           NW,NE,SE,SW,                   //for u velocities
           C_NE,C_NW,C_N,C_S,C_SE,C_SW;   //for colocated variables(mu & rho) 
  double   vP,
           vN,vS,vE,vW,                    
           uNW,uNE,uSE,uSW,                
           mu_N,mu_NW,mu_NE,mu_S,mu_SE,mu_SW, //actual mu values
           mu_W,mu_E,                         //interpolated mu values 
           rho_N,rho_S;
  

  for(l = 0;l<N_cells_v;l++)
  {
    if(v->bc_type[l] == NONE)
    {
        i = l%Nx_v;
        j = (int) l/Nx_v;
  
        S = (j-1)*Nx_v + i;
        N = (j+1)*Nx_v + i;
        W = j*Nx_v + (i-1);      // or l-1
        E = j*Nx_v + (i+1);      // or l+1
  
        SE =  l - j;
        SW = SE - 1;
        NE = SE + Nx_u;
        NW = NE - 1;

        C_S  = l ; 
        C_SW = C_S - 1;
        C_SE = C_S + 1;
        C_N  = C_S + Nx_C;
        C_NE = C_N + 1;
        C_NW = C_N - 1;     //for colocated variables(mu & rho) 
  
        vP = v->val[l];
        
        vE = v->val[E];
        vW = v->val[W];
        vN = v->val[N];
        vS = v->val[S];
  
        uNE = u->val[NE];
        uNW = u->val[NW];
        uSE = u->val[SE];
        uSW = u->val[SW];
        
        mu_NE = mu->val[C_NE];    
        mu_NW = mu->val[C_NW];    
        mu_N  = mu->val[C_N];  
        mu_S  = mu->val[C_S]; 
        mu_SE = mu->val[C_SE];  
        mu_SW = mu->val[C_SW];        //actual mu values
        
        mu_W  = 0.25*( mu_N  + mu_NW + mu_S  + mu_SW  );
        mu_E  = 0.25*( mu_N  + mu_NE + mu_SE + mu_S   );                        //interpolated mu values 
        
        rho_N = rho->val[C_N];
        rho_S = rho->val[C_S];

        
        if( i == 1 || i == (Nx_v-2) || j == 1 || j == (Ny_v-2) ) // Next to boundary 
        {
            if(v->bc_type[W] != NONE)
            {
              vW  = (2* v->BC_Value[XMIN]) - vP;
              uNW = u->BC_Value[XMIN];
              uSW = u->BC_Value[XMIN];
            }
            if(v->bc_type[E] != NONE)
            {
              vE  = (2* v->BC_Value[XMAX]) - vP;
              uSE = u->BC_Value[XMAX];
              uNE = u->BC_Value[XMAX];
            }
            if(v->bc_type[N] != NONE)
            { 
  
              vN  = v->BC_Value[YMAX];
  
            }
            if(v->bc_type[S] != NONE)
            {
  
              vS  = v->BC_Value[YMIN];
  
            }
        }

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
