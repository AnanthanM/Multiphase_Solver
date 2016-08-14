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

  int row;

  double u_p,
         u_e, u_w, u_n, u_s,
         v_ne,v_nw,v_se,v_sw;

  double * val_u_p,
         * val_u_e, * val_u_w, * val_u_n, * val_u_s,
         * val_v_ne,* val_v_nw,* val_v_se,* val_v_sw;
 
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

    val_u_p = &u->val[row];
    val_u_e = &u->val[row+1];
    val_u_w = &u->val[row-1];
    val_u_n = &u->val[row+Nx_u];
    val_u_s = &u->val[row-Nx_u];

    val_v_nw = &v->val[row+j];
    val_v_ne = &v->val[row+j+1];
    val_v_sw = &v->val[row+j-Nx_v];
    val_v_se = &v->val[row+j-Nx_v+1];
    
    for(i=1;i<(Nx_u-1);i++)
    {
      l = j*Nx_u + i;
     
      u_p = val_u_p[i];
      u_e = val_u_e[i];
      u_w = val_u_w[i];
      u_n = val_u_n[i];
      u_s = val_u_s[i];

      v_nw = val_v_nw[i];
      v_ne = val_v_ne[i];
      v_sw = val_v_sw[i];
      v_se = val_v_se[i];

      u_T->val[l] = -1*( (1/dx)*( ((u_e+u_p)/2)*((u_e+u_p)/2) - ((u_p+u_w)/2)*((u_p+u_w)/2) ) + 
                           (1/dy)*( ((u_n+u_p)/2)*((v_ne+v_nw)/2) - ((u_p+u_s)/2)*((v_se+v_sw)/2) ) ) ;

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

  int row;

  double v_p,
         v_e, v_w, v_n, v_s,
         u_ne,u_nw,u_se,u_sw;
  
  double * val_v_p,
         * val_v_e, * val_v_w, * val_v_n, * val_v_s,
         * val_u_ne,* val_u_nw,* val_u_se,* val_u_sw;
  
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

    val_v_p = &v->val[row];
    val_v_e = &v->val[row+1];
    val_v_w = &v->val[row-1];
    val_v_n = &v->val[row+Nx_v];
    val_v_s = &v->val[row-Nx_v];

    val_u_ne = &u->val[row+Nx_u-j];
    val_u_nw = &u->val[row+Nx_u-j-1];
    val_u_se = &u->val[row-j];
    val_u_sw = &u->val[row-j-1];
    
    for(i=1;i<(Nx_v-1);i++)
    {
      l = j*Nx_v + i;
     
      v_p = val_v_p[i];
      v_e = val_v_e[i];
      v_w = val_v_w[i];
      v_n = val_v_n[i];
      v_s = val_v_s[i];

      u_nw = val_u_nw[i];
      u_ne = val_u_ne[i];
      u_sw = val_u_sw[i];
      u_se = val_u_se[i];

      v_T->val[l] = -1*(  (1/dx)*( ((u_se+u_ne)/2)*((v_e+v_p)/2) - ((u_nw+u_sw)/2)*((v_p+v_w)/2) ) + 
                          (1/dy)*( ((v_n+v_p)/2)*((v_n+v_p)/2)   - ((v_p+v_s)/2)*((v_p+v_s)/2) ) );

    }
  }

  return;

}
