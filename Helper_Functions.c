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

// Allocates struct Bicgstab which contains arrays necessary
// for BiCGSTAB Solver, Takes input as an integer

Bicgstab * Allocate_Bicgstab(int N)
{
  Bicgstab * phi;

  phi             = malloc(sizeof(Bicgstab));

  phi->r0_cap     = malloc(N*sizeof(double));
  phi->sj         = malloc(N*sizeof(double));
  phi->rj         = malloc(N*sizeof(double));
  phi->pj         = malloc(N*sizeof(double));
  phi->pcap       = malloc(N*sizeof(double));
  phi->scap       = malloc(N*sizeof(double));
  phi->Ax_vector  = malloc(N*sizeof(double));
  phi->h          = malloc(N*sizeof(double));
  phi->vj         = malloc(N*sizeof(double));

  return phi;
}

// Intermediate velocities ustar and vstar are obtained

void Get_VelocityStar(Domain domain,Constant constant,Field * ustar,Field * vstar)
{
  
  Field * u   = domain.u;
  Field * v   = domain.v;

  double dt = constant.dt;

  int Nx_u    = u->N_x;
  int Ny_u    = u->N_y;
  int Nx_v    = v->N_x; 
  int Ny_v    = v->N_y;

  int i,j,l;


  // Loops for all the inner cells for u velocity 
  
  for(j=1;j<(Ny_u-1);j++)
  {
     
    for(i=1;i<(Nx_u-1);i++)
    {
      l = j*Nx_u + i;

      u->val[l] = u->val[l] + dt*(ustar->val[l]) ;    
    }
  }

  // Loops for all the inner cells for velocity v 
  
  for(j=1;j<(Ny_v-1);j++)
  {

    for(i=1;i<(Nx_v-1);i++)
    {
      l = j*Nx_v + i;
     
      v->val[l] = v->val[l] + dt*(vstar->val[l]) ;
    }
  }
  
  return;
}

//The Final velocity un+1 is obtained 
void Correction_Velocities(Domain domain,Constant constant)
{
  Field * u   = domain.u;
  Field * v   = domain.v;
  Field * p   = domain.p;
  Field * rho = domain.rho;

  double dx = constant.dx;
  double dy = constant.dy;
  double dt = constant.dt;

  int Nx_u    = u->N_x;
  int Ny_u    = u->N_y;
  int Nx_v    = v->N_x; 
  int Ny_v    = v->N_y;
  int N_x_p   = p->N_x;

  int i,j,l;

  int row;
 
  double pN,pS;  // for v
  double pE,pW;  // for u

  double rhoN,rhoS;  // for v
  double rhoE,rhoW;  // for u

  double * val_pN,* val_pS;  // for v
  double * val_pE,* val_pW;  // for u

  double * val_rhoN,* val_rhoS;  // for v
  double * val_rhoE,* val_rhoW;  // for u
  
  // Loops for all the inner cells for u velocity 
  
  for(j=1;j<(Ny_u-1);j++)
  {
    row = j*Nx_u;
    
    val_pW = &p->val[row+j];
    val_pE = &p->val[row+j+1];
    
    val_rhoW  = &rho->val[row+j];
    val_rhoE  = &rho->val[row+j+1];
    
    for(i=1;i<(Nx_u-1);i++)
    {
      l = j*Nx_u + i;

      pW = val_pW[i];
      pE = val_pE[i];
                   
      rhoW = val_rhoW[i];
      rhoE = val_rhoE[i];

      u->val[l] = u->val[l] - ( dt/(0.5*(rhoE+rhoW)) )*( (pE-pW)/dx );
    }
  }


  // Loops for all the inner cells for velocity v 
  
  for(j=1;j<(Ny_v-1);j++)
  {
    row = j*Nx_v;

    val_rhoS  = &rho->val[row];
    val_rhoN  = &rho->val[row+N_x_p];

    val_pS  = &p->val[row];
    val_pN  = &p->val[row+N_x_p];
    
    for(i=1;i<(Nx_v-1);i++)
    {
      l = j*Nx_v + i;
     
      pS  = val_pS[i]; 
      pN  = val_pN[i]; 

      rhoS  = val_rhoS[i]; 
      rhoN  = val_rhoN[i]; 

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


//Writing data at x and y crosssection of the domain
//so that results can be plotted along with Ghia etal paper

void Validation_Data(Domain domain, Constant constant)
{
  char filename1[30];
  char filename2[30];

  sprintf(filename1, "line_data_u.dat");
  sprintf(filename2, "line_data_v.dat");

  FILE *fp1 = fopen(filename1, "w");
  FILE *fp2 = fopen(filename2, "w");
  
  int N_y = domain.u_C->N_y;
  int N_x = domain.u_C->N_x;
  
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

//Continuity Test
//After each sucessfull iteration un+1 values are used to
//check the continuity eaquation

double Continuity_Test(double * DivU,Domain domain, Constant constant)
{
  Field * u = domain.u;
  Field * v = domain.v;
  Field * p = domain.p; 

  double dx = constant.dx;
  double dy = constant.dy;
  
  int Nx_v      = v->N_x;
  int N_cells_C = p->N;
  int Nx_C      = p->N_x;
  int Ny_C      = p->N_y;

/********************************/
//           --- vN-----
//          |     |    |
//          uW---pij---uE
//          |     |    |
//          -----vS-----
/********************************/
  int l,i,j;

  int row;

  double uE,uW,vN,vS;
 
  double * val_uE,* val_uW,* val_vN,* val_vS;
  
  double norm = 0.0;


  for(l = 0;l<N_cells_C;l++)
  {
    DivU[l] = 0.0;
  }

  for(j=1;j<(Ny_C-1);j++)
  {
    row = j*Nx_C;              // interms of colocated grid

    val_uE = &u->val[row-j];
    val_uW = &u->val[row-j-1];
    val_vN = &v->val[row];       
    val_vS = &v->val[row-Nx_v];

    for(i=1;i<(Nx_C-1);i++)
    {
      l = j*Nx_C + i;

      uE = val_uE[i];
      uW = val_uW[i];
      vN = val_vN[i];
      vS = val_vS[i];

      DivU[l] = ( (uE - uW)/dx + (vN - vS)/dy ) ;

    }
  }

  for(l=0;l<N_cells_C;l++)
    norm += DivU[l]*DivU[l];

  norm = sqrt(norm/((double)N_cells_C)) ;

  return norm;

}














