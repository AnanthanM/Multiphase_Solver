#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

void set_ghost_cells_type(Domain domain)
{
  int i,j,l;

  Field * p   = domain.p;
  Field * u   = domain.u;
  Field * v   = domain.v;
  Field * rho = domain.rho;
  Field * mu  = domain.mu;
  Field * C   = domain.C;
  Field * nx  = domain.nx;
  Field * ny  = domain.ny;

  int N_Cells_x = p->N_x;
  int N_Cells_y = p->N_y;
  int N_Cells   = N_Cells_x*N_Cells_y;
  
  for(l=0;l<N_Cells;l++)
  {
    i = l%N_Cells_x;
    j = (int) l/N_Cells_x;

    if( i == 0 || i == (N_Cells_x-1) || j == 0 || j == (N_Cells_y-1) )
    {
      p->bc_type[l]   = NEUMANN;
      mu->bc_type[l]  = NEUMANN;
      rho->bc_type[l] = NEUMANN;
      u->bc_type[l]   = DIRICHLET;
      v->bc_type[l]   = DIRICHLET;     
      C->bc_type[l]   = DIRICHLET;     
      nx->bc_type[l]  = DIRICHLET;     
      ny->bc_type[l]  = DIRICHLET;     
    }
    else
    {
      p->bc_type[l]   = NONE; 
      mu->bc_type[l]  = NONE;
      rho->bc_type[l] = NONE;
      u->bc_type[l]   = NONE;
      v->bc_type[l]   = NONE;    
      C->bc_type[l]   = NONE;     
      nx->bc_type[l]  = NONE;     
      ny->bc_type[l]  = NONE;     
    }
  }

  return;

}

void set_ghost_cells_value(Field * phi)
{
  int i,j,l;

  int N_Cells_x = phi->N_x;
  int N_Cells_y = phi->N_y;
  int N_Cells   = N_Cells_x*N_Cells_y;
  
  for(l=0;l<N_Cells;l++)
  {
    i = l%N_Cells_x;
    j = (int) l/N_Cells_x;
    
    if(phi->bc_type[l] != NONE)
    {
      if(i==0)
        phi->val[l] = phi->BC_Value[XMIN];
      else if(i==(N_Cells_x-1))
        phi->val[l] = phi->BC_Value[XMAX];
      if(j==0)
        phi->val[l] = phi->BC_Value[YMIN];
      else if(j==(N_Cells_y-1))
        phi->val[l] = phi->BC_Value[YMAX];
    }
  }

  return;
}
