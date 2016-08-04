#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"


int main(int argc, char *argv[])
{
  Constant constant;
  Domain   domain; 

  double L_x = 1.0;
  double L_y = 1.0;
  
  int N_cells_x = 100+2;
  int N_cells_y = 100+2;
  int N_cells_z = 1;
  int N_cells   = N_cells_x * N_cells_y * N_cells_z ;

  constant.L_x = L_x;
  constant.L_y = L_y; 
  constant.dx  = L_x/(N_cells_x-2);
  constant.dy  = L_y/(N_cells_y-2);
  constant.dz  = 1.0;

  constant.dt = 0.02;

  Field * p    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * u    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * v    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * rho  = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * mu   = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * C    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * nx   = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * ny   = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  
  domain.p   = p;
  domain.u   = u;
  domain.v   = v;
  domain.rho = rho;
  domain.mu  = mu;
  domain.C   = C;
  domain.nx  = nx;
  domain.ny  = ny;

  set_ghost_cells_type(domain);

  int i;
// i = 0 -> XMIN -> left
// i = 1 -> XMAX -> right
// i = 2 -> YMIN -> bottom
// i = 3 -> YMAX -> top
  for(i=0;i<4;i++)
  {
    p->BC_Value[i]   = 1000;
    rho->BC_Value[i] = 100;
    mu->BC_Value[i]  = 5;
    u->BC_Value[i]   = 10;
    v->BC_Value[i]   = 15;
    C->BC_Value[i]   = 0.0;
    nx->BC_Value[i]  = 0;
    ny->BC_Value[i]  = 0;
  }

  p->BC_Value[XMAX] = 10000;

  set_ghost_cells_value(p);
  set_ghost_cells_value(u);
  set_ghost_cells_value(v);
  set_ghost_cells_value(rho);
  set_ghost_cells_value(mu);
  set_ghost_cells_value(C);
  set_ghost_cells_value(nx);
  set_ghost_cells_value(ny);

  Initiation_Void_Fraction(domain,constant);

  for(i =0;i < N_cells;i++)
  {
    if(p->bc_type[i] == NONE)
    {
      p->val[i]   = 100.0; 
      u->val[i]   = 1.0; 
      v->val[i]   = 1.5; 
//      C->val[i]   = Void_Fraction already initialised 
      rho->val[i] = 10.0*C->val[i] + 1.0*(1-C->val[i]);  
      mu->val[i]  = 3.0*C->val[i] + 0.3*(1-C->val[i]); 
      nx->val[i]  = 4.0; 
      ny->val[i]  = 3.0; 
    }
  }
  
  double final_time = 1.0;
  double time       = 0.0;

  int si_no =0;
  while(time < (final_time+constant.dt/2) ) 
  {
    Write_VTK(si_no,domain,constant);
    printf("At time %2.8lf VTK file is written \n",time);
    time += constant.dt;
    si_no ++;
  }

  return 0;
}



















