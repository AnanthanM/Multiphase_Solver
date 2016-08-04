#ifndef FVM_H_
#define FVM_H_

typedef struct
{
  double L_x,L_y;
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
  Field * C   ;
  Field * nx  ;
  Field * ny  ;
} Domain;

void  Write_VTK(int,Domain,Constant);

Field * Allocate_Field(int,int,int);

void set_ghost_cells_type(Domain);

void set_ghost_cells_value(Field *);

void Initiation_Void_Fraction(Domain,Constant);

void Normals_Using_Youngs_Method(Domain,Constant);

#endif
