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
  ON_BOUNDARY,
  INSIDE_GC,
  NONE
} BC_Type;

typedef enum
{
  XMIN,
  XMAX,
  YMIN,
  YMAX
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
  Field * u_C ;
  Field * v_C ;
} Domain;

void  Write_VTK(int,Domain,Constant);

Field * Allocate_Field(int,int,int);

void set_ghost_cells_type(Domain);

void set_ghost_cells_value(Field *);

void Initiation_Void_Fraction(Domain,Constant);

void Normals_Using_Youngs_Method(Domain,Constant);

void Advection_u(Domain,Constant,Field *);

void Advection_v(Domain,Constant,Field *);

void Diffusion_u(Domain,Constant,Field *);

void Diffusion_v(Domain,Constant,Field *);

int  solve_Pressure_Poisson_BiCGSTAB(Field *,Field *,Constant,Domain); 

void Get_RHS_of_Pressure_Poisson_Eqn(Domain,Field *,Constant);

void Pressure_Poisson(Field *,Constant,double *,Domain);

void Correction_Velocities(Domain,Constant);

void Making_u_v_Collocated(Domain,Constant);

#endif
