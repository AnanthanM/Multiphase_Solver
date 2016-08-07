#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

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



















