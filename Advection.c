/*****************Advection*****************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

/********* Advection of velocity u   ********************
 *  
 *  the advected u velocity is stored in u_T, a temporary array
 *
 * *********************************************************/

void Advection_u(Domain domain, Constant constant, Field * u_T)
{
  Field * u = domain.u;
  Field * v = domain.v;

  double dx = constant.dx;
  double dy = constant.dy;

  int N_cells_u    = u->N;
  int Nx_u         = u->N_x;
  int Ny_u         = u->N_y;
  int Nx_v         = v->N_x; 
  int Ny_v         = v->N_y; 

/********************************/
           
//            un
// vij = vnw_ |_ vne                 
//        |   |   |
//    uw - - uij - - ue 
//        |   |   |               
//       vsw  |  vse               
//            us               

/*******************************/  
  int i,j,l;

  int      E,  W,  N,  S,      // for u velocities   
          NE, NW, SE, SW;      // for v velocities
  double u_p,
         u_e, u_w, u_n, u_s,
         v_ne,v_nw,v_se,v_sw;

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
  
        u_p = u->val[l];
        
        u_e = u->val[E];
        u_w = u->val[W];
        u_n = u->val[N];
        u_s = u->val[S];
  
        v_ne = v->val[NE];
        v_nw = v->val[NW];
        v_se = v->val[SE];
        v_sw = v->val[SW];
        
        if( i == 1 || i == (Nx_u-2) || j == 1 || j == (Ny_u-2) ) // Next to boundary 
        {
            if(u->bc_type[N] != NONE)
            {
              u_n  = (2* u->BC_Value[YMAX]) - u_p;
              v_nw = v->BC_Value[YMAX];
              v_ne = v->BC_Value[YMAX];
            }
            if(u->bc_type[S] != NONE)
            {
              u_s  = (2* u->BC_Value[YMIN]) - u_p;
              v_sw = v->BC_Value[YMIN];
              v_se = v->BC_Value[YMIN];
            }
            if(u->bc_type[W] != NONE)
            { 
  
              u_w  = u->BC_Value[XMIN];
  
            }
            if(u->bc_type[E] != NONE)
            {
  
              u_e  = u->BC_Value[XMAX];
  
            }
        }
  
        u_T->val[l] = (1/dx)*( ((u_e+u_p)/2)*((u_e+u_p)/2) - ((u_p+u_w)/2)*((u_p+u_w)/2) ) + 
                      (1/dy)*( ((u_n+u_p)/2)*((v_ne+v_nw)/2) - ((u_p+u_s)/2)*((v_se+v_sw)/2) );
      
    }
  }

  return;

}


/********* Advection of velocity v   ********************
 *  
 *  the advected v velocity is stored in v_T, a temporary array
 *
 * *********************************************************/

void Advection_v(Domain domain, Constant constant, Field * v_T)
{
  Field * v = domain.v;
  Field * u = domain.u;

  double dx = constant.dx;
  double dy = constant.dy;

  int N_cells_v    = v->N;
  int Nx_v         = v->N_x;
  int Ny_v         = v->N_y;
  int Nx_u         = u->N_x; 
  int Ny_u         = u->N_y; 

/********************************/
           
//            vn
//       unw_ |_ une                 
//        |   |   |
//    vw - - vij - - ve 
//        |   |   |               
//       usw  |  use = uij              
//            vs               

/*******************************/  
  int i,j,l;
  int      E,  W,  N,  S,      // for v velocities   
          NE, NW, SE, SW;      // for u velocities
  double v_p,
         v_e, v_w, v_n, v_s,
         u_ne,u_nw,u_se,u_sw;

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

        v_p = v->val[l];
        
        v_e = v->val[E];
        v_w = v->val[W];
        v_n = v->val[N];
        v_s = v->val[S];
  
        u_ne = u->val[NE];
        u_nw = u->val[NW];
        u_se = u->val[SE];
        u_sw = u->val[SW];
        
        if( i == 1 || i == (Nx_v-2) || j == 1 || j == (Ny_v-2) ) // Next to boundary 
        {
            if(v->bc_type[W] != NONE)
            {
              v_w  = (2* v->BC_Value[XMIN]) - v_p;
              u_nw = u->BC_Value[XMIN];
              u_sw = u->BC_Value[XMIN];
            }
            if(v->bc_type[E] != NONE)
            {
              v_e  = (2* v->BC_Value[XMAX]) - v_p;
              u_se = u->BC_Value[XMAX];
              u_ne = u->BC_Value[XMAX];
            }
            if(v->bc_type[N] != NONE)
            { 
  
              v_n  = v->BC_Value[YMAX];
  
            }
            if(v->bc_type[S] != NONE)
            {
  
              v_s  = v->BC_Value[YMIN];
  
            }
        }
  
        v_T->val[l] = (1/dx)*( ((u_se+u_ne)/2)*((v_e+v_p)/2) - ((u_nw+u_sw)/2)*((v_p+v_w)/2) ) + 
                      (1/dy)*( ((v_n+v_p)/2)*((v_n+v_p)/2) - ((v_p+v_s)/2)*((v_p+v_s)/2) );
      
    }
  }

  return;

}
