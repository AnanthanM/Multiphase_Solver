#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

#define RHO1 1      //For which C = 1
#define RHO2 1      //For which C = 0
#define MU1  0.0025 //For which C = 1
#define MU2  0.0025 //For which C = 0
#define g_x  0.0    //Body forces in x-direction
#define g_y  10.0   //Body forces in y direction


int main(int argc, char *argv[])
{
  Constant constant;
  Domain   domain; 

  double L_x = 1.0;
  double L_y = 1.0;
  
  int N_cells_x   = 5+2; //For the entire domain i.e for p,mu.rho.C,nx,ny
  int N_cells_x_u = N_cells_x -1;//For staggered u velocity
  int N_cells_y   = 5+2; //For the entire domain i.e for p,mu.rho.C,nx,ny
  int N_cells_y_v = N_cells_y -1;//For staggered v velocity
  int N_cells_z   = 1;
  int N_cells     = N_cells_x * N_cells_y * N_cells_z ;

  constant.L_x = L_x;
  constant.L_y = L_y; 
  constant.dx  = L_x/(N_cells_x-2);
  constant.dy  = L_y/(N_cells_y-2);
  constant.dz  = 1.0;

  constant.dt = 0.02;

  Field * p    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );

  Field * u    = Allocate_Field( N_cells_x_u, N_cells_y, N_cells_z );
  Field * u_T  = Allocate_Field( N_cells_x_u, N_cells_y, N_cells_z ); 
  
  Field * v    = Allocate_Field( N_cells_x, N_cells_y_v, N_cells_z );
  Field * v_T  = Allocate_Field( N_cells_x, N_cells_y_v, N_cells_z ); 
  
  Field * u_C  = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );    // Collocated u and v values 
  Field * v_C  = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );    // for plotting purposes
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
  domain.u_C = u_C;
  domain.v_C = v_C;

  set_ghost_cells_type(domain);

  int i;
// i = 0 -> XMIN -> left
// i = 1 -> XMAX -> right
// i = 2 -> YMIN -> bottom
// i = 3 -> YMAX -> top
  for(i=0;i<4;i++)
  {
//  For Staggered grid we dont need value of pressure at ghost cells
    rho->BC_Value[i] = RHO2;
    mu->BC_Value[i]  = MU2;
    u->BC_Value[i]   = 0.0;
    v->BC_Value[i]   = 0.0;
    C->BC_Value[i]   = 0.0;
    nx->BC_Value[i]  = 0.0;
    ny->BC_Value[i]  = 0.0;
  }
  
  v->BC_Value[YMIN] = 1;

// GCs values not needed for p
  set_ghost_cells_value(u);
  set_ghost_cells_value(v);
  set_ghost_cells_value(rho);
  set_ghost_cells_value(mu);
  set_ghost_cells_value(C);
  set_ghost_cells_value(nx);
  set_ghost_cells_value(ny);

//  Initiation_Void_Fraction(domain,constant);

// Normals_Using_Youngs_Method(domain,constant);


  for(i =0;i < N_cells;i++)
  {
    if(p->bc_type[i] == NONE)
    {
      p->val[i]   = 0.0; 
//      C->val[i]   = Void_Fraction already initialised 
      C->val[i]   = 0.0; 
      rho->val[i] = RHO1*C->val[i] + RHO2*(1-C->val[i]);  
      mu->val[i]  = MU1*C->val[i] + MU2*(1-C->val[i]); 
//      nx->val[i]  = Normal_x already found 
//      ny->val[i]  = Normal_y already found
      nx->val[i]  = 0.0; 
      ny->val[i]  = 0.0;
    }
  }

  int N_cells_uv = u->N;
  for(i=0;i < N_cells_uv;i++)
  {
    if(u->bc_type[i] == NONE)
      u->val[i] = 0.0;
    if(v->bc_type[i] == NONE)
      v->val[i] = 0.0;
  }

  
  double final_time = 1.0;
  double time       = 0.0;

/******STARTING*******/  
  Write_VTK(0,domain,constant);
  
  int si_no     = 1;
  int N_cells_u = u->N;
  int N_cells_v = v->N;
  double dt     = constant.dt;

  while(time < (final_time+dt/2) ) 
  { 

    /********ADVECTION*******************/
    for(i=0;i<N_cells_u;i++)
        u_T->val[i] = 0.0;
    Advection_u(domain,constant,u_T);

    for(i=0;i<N_cells_v;i++)
        v_T->val[i] = 0.0;
    Advection_v(domain,constant,v_T);
     
    for(i=0;i<N_cells_u;i++)
    {
       if(u->bc_type[i] == NONE)
         u->val[i] = u->val[i] + dt*(u_T->val[i]) ;    
       if(v->bc_type[i] == NONE)
         v->val[i] = v->val[i] + dt*(v_T->val[i]) ;
    }
    
    Write_VTK(si_no,domain,constant);
    
    printf("At time %2.8lf VTK file is written \n",time);
    
    time += dt;
    
    si_no ++;
  }

  return 0;
}



















