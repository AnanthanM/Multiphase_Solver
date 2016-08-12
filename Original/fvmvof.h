#ifndef FVMVOF_H_
#define FVMVOF_H_

#define DELTA 1e-10
#define PI 2.0*acos(0.0)

#define EAST i+1
#define EEAST i+2
#define WEST i-1
#define WWEST i-2
#define NORTH i+N_x
#define NNORTH i+2*N_x
#define SOUTH i-N_x
#define SSOUTH i-2*N_x
#define BOTTOM i-N_x*N_y
#define BBOTTOM i-2*N_x*N_y
#define TOP i+N_x*N_y
#define TTOP i+2*N_x*N_y
#define P i

typedef enum{
  INSIDE,
  XMIN,
  XMAX,
  YMIN,
  YMAX,
  ZMIN,
  ZMAX,
  SOLID
} PATCH_type;
PATCH_type patch;

typedef enum{
  CENTERED,
  X_STAGGERED,
  Y_STAGGERED
}data_location;

typedef enum{
  NONE,
  DIRICHLET,
  NEUMANN,
} BC_type;

typedef enum{
  CENTRAL,
  MIM,
  UPWIND,
  QUICK
} Discretization_type;

typedef enum{
  POISSON,
  HELMHOLTZ
} Equation_type;

typedef struct _Field Field;
struct _Field{
  data_location grid;  
  int N_x, N_y, N_z;
  int N;
  double bc_val[8];
  BC_type *bc;
  double *val;
};

typedef struct _Domain Domain;
struct _Domain{
  Field * p;
  Field * u_x;
  Field * u_y;
  Field * u_z;
  Field * phi;
};

typedef struct _Constant Constant;
struct _Constant{
  double dx, dy, dz;  
  double dt;
  double nu;
  double rho;
};

double mu, CFL, p_ref, rho_ref, u_ref;

/* Declare Functions */
void advection(Field *, Field *, Field *, Field *, Field *, Constant );
void diffusion(Field * ,double ,Field *, Constant );
void diffusion_implicit( Field * phi, Constant constant,double * tmp );
void set_ghosts(Domain);
void set_bc(Field * phi);
void divergence(Field * div, Field * u_x, Field * u_y, Constant constant);

void laplacian( Field * phi, Constant constant,
    double * tmp
    );
void gradient(Field * p, Field * temp_x, Field * temp_y, Constant constant);

#endif
