/**************************************************************************/
/**     A Finite Volume Based 1 Fluid Multiphase Flow Solver              */
/**     -> Structured Uniform Grid, Staggered 
 *      -> Projection Method 
 *      -> For Initialising Void Fraction an External Library :VOFI
 *      ->Volume of Fluid for Tracking Volume with PLIC
 *      ->Our own internal BiCGSTAB Solver for Solving Pressure Poisson Equation
 *      ->Explicit Advection and Diffusion 
 *
 *      Written By : Ananthan M
 *                 : Dr. Gaurav Tomar
 *                 Indian Institute of Science, Mechanical Engg. Dept
 *  
 ****************************************************************************/  

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include <omp.h>

#include "fvm.h"

#define RHO1     1      //For which C = 1
#define RHO2     1      //For which C = 0
#define MU1      0.001  //For which C = 1
#define MU2      0.001  //For which C = 0
#define g_x      0.0    //Body forces in x-direction
#define g_y      10.0   //Body forces in y direction


int main(int argc, char *argv[])
{
  double start_wall_time = omp_get_wtime();
  
  Constant constant;
  Domain   domain; 

  double L_x = 1.0;                                                    // Domain Dimensions
  double L_y = 1.0;

  double epsilon  = 1.0E-7;

  int N_cells_x   = 140+2;                                            //For the entire domain for colocated variables
  int N_cells_x_u = N_cells_x - 1;                                     //For staggered u velocity
  int N_cells_y   = 140+2;                                            //For the entire domain for colocated variablesy
  int N_cells_y_v = N_cells_y - 1;                                     //For staggered v velocity
  int N_cells_z   = 1;
  int N_cells     = N_cells_x * N_cells_y * N_cells_z ;               //Total Number of Colocated cells

  /******Initialising the Constants **************/
  
  constant.L_x = L_x;
  constant.L_y = L_y; 
  constant.dx  = L_x/(N_cells_x-2);
  constant.dy  = L_y/(N_cells_y-2);
  constant.dz  = 1.0;
  constant.dt  = 0.002;

  /*****Initialising the Field Structs**********/

  // p   -> Pressure
  // u   -> horizontal velocity
  // v   -> vertical velocity
  // rho -> density
  // mu  -> viscosity
  // C   -> Color function/ void fraction
  // nx  -> Normal component of interface in x
  // ny  -> Normal component of interface in y
  

  Field * p      = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * RHS    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );   //Acts as RHS of Pressure Poisson equation

  Field * u      = Allocate_Field( N_cells_x_u, N_cells_y, N_cells_z ); //Note no of cells in x direction is different 
  Field * ustar  = Allocate_Field( N_cells_x_u, N_cells_y, N_cells_z ); 
  
  Field * v      = Allocate_Field( N_cells_x, N_cells_y_v, N_cells_z ); //Note no of cells in y direction is different
  Field * vstar  = Allocate_Field( N_cells_x, N_cells_y_v, N_cells_z ); 
  
  Field * u_C    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );   // Collocated u and v values 
  Field * v_C    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );   // for plotting purposes

  Field * rho    = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );   // All other colocated variables 
  Field * mu     = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * C      = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * nx     = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  Field * ny     = Allocate_Field( N_cells_x, N_cells_y, N_cells_z );
  
  Bicgstab * PP = Allocate_Bicgstab( N_cells ); 
 
  /*****Initialising the Domain struct***************/

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

  domain.PP  = PP;
  
  /****Ghost cell types are set for each Field*********/

  set_ghost_cells_type(domain);

  /****Boundary Cell values are set for each Field******/

  int i;
// i = 0 -> XMIN -> left
// i = 1 -> XMAX -> right
// i = 2 -> YMIN -> bottom
// i = 3 -> YMAX -> top
  for(i=0;i<4;i++)
  {
    p->BC_Value[i]   = 0.0;
    rho->BC_Value[i] = RHO2;
    mu->BC_Value[i]  = MU2;
    u->BC_Value[i]   = 0.0;
    v->BC_Value[i]   = 0.0;
    C->BC_Value[i]   = 0.0;
    nx->BC_Value[i]  = 0.0;
    ny->BC_Value[i]  = 0.0;
  }
  
  u->BC_Value[YMAX] = 1;         // Means that u velocity would have  
                                 // magnitude 1 in top horizontal side

  set_ghost_cells_value(p);
  set_ghost_cells_value(u);
  set_ghost_cells_value(v);
  set_ghost_cells_value(rho);
  set_ghost_cells_value(mu);
  set_ghost_cells_value(C);
  set_ghost_cells_value(nx);
  set_ghost_cells_value(ny);

// Initiation_Void_Fraction(domain,constant);    // we use VOFI library to 
                                                 // initialise C Field

// Normals_Using_Youngs_Method(domain,constant); //We call a function to 
                                                 // initialise normal fields 


  for(i =0;i < N_cells;i++)
  {
    if(p->bc_type[i] == NONE)
    {
      p->val[i]   = 0.0; 
//      C->val[i]   = Void_Fraction already initialised 
      C->val[i]   = 0.0; 
//      rho->val[i] = RHO1*C->val[i] + RHO2*(1-C->val[i]);  
//      mu->val[i]  = MU1*C->val[i] + MU2*(1-C->val[i]); 
      rho->val[i] = RHO2;  
      mu->val[i]  = MU2; 
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
 
  double final_time =  66.0;
  double time       =  0.0;

/******STARTING*******/  
  Write_VTK(0,domain,constant);                    //Writing Output at the time t=0
  
  int si_no     = 1;
  int test;      
  int N_cells_u = u->N;
  int N_cells_v = v->N;
  double dt     = constant.dt;
  double norm;

  while(time < (final_time+dt/2) ) 
  { 

    for(i=0;i<N_cells_u;i++)
        ustar->val[i] = 0.0;

    for(i=0;i<N_cells_v;i++)
        vstar->val[i] = 0.0;

    /********ADVECTION*******************/
    Advection_u(domain,constant,ustar);
    Advection_v(domain,constant,vstar);
    
    /********DIFFUSION********************/
    Diffusion_u(domain,constant,ustar);
    Diffusion_v(domain,constant,vstar);

    for(i=0;i<N_cells_u;i++)
    {
       if(u->bc_type[i] == NONE)
         u->val[i] = u->val[i] + dt*(ustar->val[i]) ;    
       if(v->bc_type[i] == NONE)
         v->val[i] = v->val[i] + dt*(vstar->val[i]) ;
    }
    
    /******Solving Pressure Poisson ********/
    Get_RHS_of_Pressure_Poisson_Eqn(domain,RHS,constant);
    test = solve_Pressure_Poisson_BiCGSTAB(p,RHS,constant,domain);
    if(test)
    {
      printf("\n Problem Solving Pressure Poisson Equation\n");
      break;
    }
    
    /********Getting Final Velocities**********/
    Correction_Velocities(domain,constant);

    /********Continuity Test*********************/
    norm = Continuity_Test(domain,constant);
    if( norm >= epsilon )
    {
      printf("At time %2.8lf ERROR N-S eqn Crashed\n",time);
      break;
    }
    printf("At time %2.8lf RMS of DivU is %2.8lf \n",time,norm);
    
    /*********Writing Output*************************/
    printf("\n At %2.8lf N-S equations solved \n ",time);
    if( (si_no-1)%5000 == 0)
    { 
      set_ghost_cells_value(p);
      set_ghost_cells_value(u);
      set_ghost_cells_value(v);
      Write_VTK(si_no,domain,constant); 
      printf("At time %2.8lf VTK file is written \n",time);
    }

    time += dt;
    
    si_no ++;
  }

  /***********Write Validation data**************/
  Validation_Data(domain,constant); 

  printf("Time taken for the simulation: %lf", omp_get_wtime() - start_wall_time);
  return 0;
}



















