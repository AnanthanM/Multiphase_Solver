#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include<vofi.h>

#include "fvm.h"

/**********Function to find Void Fraction******************************/

typedef const double creal;
typedef const int cint;
typedef double real;

//Impicit Function
real Implicit_Function(creal xy[])
{
  real  x,y,r,xc,yc,f0;
  
  r  = 0.25;
  xc = 0.5;
  yc = 0.5;

  x = xy[0];  
  y = xy[1];

  f0 = (x-xc)*(x-xc) + (y-yc)*(y-yc) - r*r;

  return f0;
}
//Finding Void Fractions Initialy Using VOFI library

void Initiation_Void_Fraction(Domain domain,Constant constant)
{
  Field * C = domain.C;

  int N_cells_x = C->N_x;
  int N_cells_y = C->N_y;
  int N_cells   = N_cells_x * N_cells_y;
  
  double L_x = constant.L_x;
  double L_y = constant.L_y;
  double dx  = constant.dx;
  double dy  = constant.dy;

  vofi_real  cc,x0[3],xloc[3];

  double fh,h0;
  int    ndim0 = 2;
  int    itrue = 1;

  if(!(dx == dy))
  {
    printf("Error in Void Fraction Initialisation \n Grid is not Cubic \n");
    return;
  }

  h0 = dx;

  /* put starting point in the center of the square */
  x0[0] = L_x/2;
  x0[1] = L_y/2; 
  x0[2] = 0.; 

  /* get the characteristic value fh of the implicit function */
  fh = vofi_Get_fh(Implicit_Function,x0,h0,ndim0,itrue);

  /* put now starting point in (X0,Y0) to initialize the color function */
  x0[0] = 0.0; 
  x0[1] = 0.0;

  /* xloc: minor vertex of each cell of the grid */
  xloc[2] = 0.;    

  int l,i,j;

  for(l=0;l<N_cells;l++)
  {
    if(C->bc_type[l] == NONE)
    {
      i = l%N_cells_x;
      j = (int) l/N_cells_x;

      xloc[0] = x0[0] + (i-1)*dx;
      xloc[1] = x0[1] + (j-1)*dy;
      
      cc        = vofi_Get_cc(Implicit_Function,xloc,h0,fh,ndim0);
      C->val[l] = cc;
    }
  }    
  
  return;
} 

/************************END*******************************************/

/*******************FIND NORMAL nx and ny using Youngs FD method*********/

// Normal_x(i,j) = (1/dx)*{ C(i+1,j+1) + 2*C(i+1,j) + C(i+1,j-1) - C(i-1,j+1) - 2*C(i-1,j) - C(i-1,j-1) }
  
// Normal_y(i,j) = (1/dy)*{ C(i+1,j+1) + 2*C(i,j+1) + C(i-1,j+1) - C(i+1,j-1) - 2*C(i,j-1) - C(i-1,j-1) } 

void Normals_Using_Youngs_Method(Domain domain,Constant constant)
{
  Field * C  = domain.C;
  Field * nx = domain.nx;
  Field * ny = domain.ny;

  int N_cells_x = C->N_x;
  int N_cells_y = C->N_y;
  int N_cells   = N_cells_x * N_cells_y;
  
  double dx  = constant.dx;
  double dy  = constant.dy;

  int l;

  double Cip1jp1,Cip1j,Cip1jm1,
         Cijp1,Cijm1,
         Cim1jp1,Cim1j,Cim1jm1;

  for(l=0;l<N_cells;l++)
  {
    if(C->bc_type[l] == NONE)
    {
      Cip1jp1 = C->val[l+N_cells_x+1];
      Cip1j   = C->val[l+1];
      Cip1jm1 = C->val[l-N_cells_x+1];
      Cijp1   = C->val[l+N_cells_x];
      Cijm1   = C->val[l-N_cells_x];
      Cim1jp1 = C->val[l+N_cells_x-1];
      Cim1j   = C->val[l-1];
      Cim1jm1 = C->val[l-N_cells_x-1];

      nx->val[l] = (1/dx)*( Cip1jp1 + 2*Cip1j + Cip1jm1 - Cim1jp1 - 2*Cim1j - Cim1jm1 );
      ny->val[l] = (1/dy)*( Cip1jp1 + 2*Cijp1 + Cim1jp1 - Cip1jm1 - 2*Cijm1 - Cim1jm1 );
    }
  }    

  return;
}























/*************************END*******************************************/
