#ifndef FVM_H_
#define FVM_H_

typedef struct
{
  double dx, dy, dz;
  double dt;
} Constant; 

typedef enum
{
  NONE,
  DIRICHLET,
  NEUMANN
} BC_Type;

typedef enum
{
  XMIN,
  XMAX,
  YMAX,
  YMIN
} BC_Position;

typedef struct
{
  int N_x,N_y,N_z;
  int N;
  BC_Type *bc_type;
  double BC_Value[4];
  double *val;
} Field;

typedef struct
{
  Field * p   ;
  Field * u   ;
  Field * v   ;
  Field * rho ;
  Field * mu  ;
} Domain;

void  Write_VTK(int,Domain,Constant);

Field * Allocate_Field(int,int,int);

void set_ghost_cells_type(Domain);

void set_ghost_cells_value(Field *);

#endif
