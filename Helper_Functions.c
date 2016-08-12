/*************** Helper Functions *****************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

// Allocate Field takes 3 integers as input and initialise 
// a Field with those values as Nx, Ny and Nz points
// in x,y and z direction.It returs the adress of the field.

Field * Allocate_Field(int N_x,int N_y,int N_z)
{
  Field * phi;

  phi           = malloc(sizeof(Field));
  phi->N_x      = N_x;
  phi->N_y      = N_y;
  phi->N_z      = N_z;
  phi->N        = N_x*N_y*N_z;
  phi->bc_type  = malloc(phi->N * sizeof(BC_Type)); 
  phi->val      = malloc(phi->N * sizeof(double));

  return phi;
}

void Correction_Velocities(Domain domain,Constant constant)
{
  Field * u   = domain.u;
  Field * v   = domain.v;
  Field * p   = domain.p;
  Field * rho = domain.rho;

  double dx = constant.dx;
  double dy = constant.dy;
  double dt = constant.dt;

  int N_cells = u->N;      // Since total no of cells for both u and v are same 
  int Nx_u    = u->N_x;
  int Ny_u    = u->N_y;
  int Nx_v    = v->N_x; 
  int Ny_v    = v->N_y;
  int N_x_p   = p->N_x;

  int l,i,j;
  int E,W;       //p values for u correction
  int N,S;       //p values for v correction
 
  double pN,pS;  // for v
  double pE,pW;  // for u

  double rhoN,rhoS;  // for v
  double rhoE,rhoW;  // for u
  
  for(l = 0;l<N_cells;l++)
  {
    if(u->bc_type[l] == NONE)
    {
        i = l%Nx_u;
        j = (int) l/Nx_u;

        W = l + j;
        E = W + 1;

        pE = p->val[E];
        pW = p->val[W];
        
        rhoE = rho->val[E];
        rhoW = rho->val[W];

        u->val[l] = u->val[l] - ( dt/(0.5*(rhoE+rhoW)) )*( (pE-pW)/dx );

    }

    if(v->bc_type[l] == NONE)
    {
        i = l%Nx_v;
        j = (int) l/Nx_v;

        S = l ;
        N = S + N_x_p;

        pN = p->val[N];
        pS = p->val[S];
        
        rhoN = rho->val[N];
        rhoS = rho->val[S];

        v->val[l] = v->val[l] - ( dt/(0.5*(rhoN+rhoS)) )*( (pN-pS)/dy );

    }    
  }

  return;

}

// The grid we are using is staggered for u and v but for ploting 
// output we need cell centered values of all u and v

void Making_u_v_Collocated(Domain domain,Constant constant)
{
  Field * u   = domain.u;
  Field * u_C = domain.u_C;
  Field * v   = domain.v;
  Field * v_C = domain.v_C;
  Field * p   = domain.p;

  int N_x_v   = v->N_x;
  int N_x_C   = p->N_x;
  int N_y_C   = p->N_y;

  int N_cells_C  = N_x_C*N_y_C; 
  
  int l,i,j;
  for(l = 0;l < N_cells_C;l++)
  {
    i = l%N_x_C;
    j = (int) l/N_x_C;

    if( i==0 )
      u_C->val[l] = u->val[l-j];
    else if( i == (N_x_C-1) )
      u_C->val[l] = u->val[l-j-1];
    else
      u_C->val[l] = 0.5*( u->val[l-j] + u->val[l-j-1] );

    if( j==0 )
      v_C->val[l] = v->val[l];
    else if ( j == (N_y_C - 1) )
      v_C->val[l] = v->val[l-N_x_v];
    else
      v_C->val[l] = 0.5*( v->val[l] + v->val[l-N_x_v] );
  }
  
  return;
}


void Validation_Data(Domain domain, Constant constant)
{
  char filename1[30];
  char filename2[30];

  sprintf(filename1, "line_data_u.dat");
  sprintf(filename2, "line_data_v.dat");

  FILE *fp1 = fopen(filename1, "w");
  FILE *fp2 = fopen(filename2, "w");
  
  int N_y = domain.u_x->N_y;
  int N_x = domain.u_y->N_x;
  
  int i,j,mid;

  mid = (N_x-1)/2 + 2;
  
  for(j = 1; j <N_y-1; j ++)
  {
    i = N_x*j+mid;
    fprintf(fp1, "%2.8lf \t %2.8lf \n", constant.dy*((j-1) + 0.5), domain.u_C->val[i]);
  }
  
  mid = (N_y-1)/2 + 2;
  
  for(j = 1; j <N_x-1; j ++)
  {
    i = mid*N_x + j;
    fprintf(fp2, "%2.8lf \t %2.8lf \n", constant.dx*((j-1) + 0.5), domain.v_C->val[i]);
  }

}

















