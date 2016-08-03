#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

Field *Allocate_Field(int N_x,int N_y,int N_z)
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
