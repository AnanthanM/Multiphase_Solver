#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include "fvmvof.h"

//void Compute_AX(double * );
//int solve_BiCGSTAB(void);
int solve_BiCGSTAB(Field *p, Field *div, Constant constant,Equation_type) ;
static void write_vtk(int, Domain, Constant);
static void write_line_data(Domain domain, Constant constant);
static Field * allocate_field( int, int, int);

int main(int argc, char *argv[])
{
#ifdef _OPENMP
double start_wall_time = omp_get_wtime();
#endif
  Domain domain;
  Constant constant;
  double l_x = 1.0;
  double l_y = 1.0;
  int N_cells_x = 127+2, N_cells_y = 127+2, N_cells_z = 1;
  int N_cells = N_cells_x *  N_cells_y * N_cells_z;
  
  constant.dx = l_x / (N_cells_x-2); constant.dy = l_y / (N_cells_y-2);
  constant.dz = 1.0; constant.dt = 0.002; constant.nu = 0.01;
  
  Field * p   = allocate_field( N_cells_x, N_cells_y, N_cells_z);
  Field * u_x = allocate_field( N_cells_x, N_cells_y, N_cells_z);
  Field * u_y = allocate_field( N_cells_x, N_cells_y, N_cells_z);
  Field * u_z = allocate_field( N_cells_x, N_cells_y, N_cells_z);
  Field * div = allocate_field( N_cells_x, N_cells_y, N_cells_z);
  Field * temp_x = allocate_field(N_cells_x, N_cells_y, N_cells_z);
  Field * temp_y = allocate_field(N_cells_x, N_cells_y, N_cells_z);
  
  int i,l,m; 
  
  for(i=0;i<8;i++){
    p->bc_val[i] = 0.0;
    u_x->bc_val[i] = 0.0;
    u_y->bc_val[i] = 0.0;
    temp_x->bc_val[i] = 0.0;
    temp_y->bc_val[i] = 0.0;
  }
  u_x->bc_val[YMAX] = 1.0;
  domain.p = p; domain.u_x = u_x; domain.u_y = u_y; 
  domain.u_z = u_z; 
  set_ghosts(domain);

  //initialization 
  for(i=0;i<N_cells;i++){
      l= i%N_cells_x;
      m =(int) i/N_cells_x;
    p->val[i] = 0.0;
    u_x->val[i] = 0.0;
    u_y->val[i] = 0.0;
    u_z->val[i] = 0.0;
    temp_x->val[i] = 0.0;
    temp_y->val[i] = 0.0;
    div->val[i] = 0.0;
  }
  
  set_bc(u_x);
  set_bc(u_y);
  set_bc(u_z);
  set_bc(p); 
/*
  double peclet = 1.0*1.0 *constant.dx/constant.nu;
  if(peclet>=2.0){
    printf("peclet number not less than 2.0 \n");
    exit(1);
  }*/
  if(constant.dt>=0.3*constant.dx/1.0){
    printf("t not within CFL criterion \n");
    exit(1);
  }
    //advection
  int qq=0, test;
  double time = 0.0;
  double end_time = 10.0;

  while (time <= end_time){
    if(qq%500 == 0)  
      write_vtk(qq, domain, constant); 

    for(i=0;i<N_cells;i++){
      temp_x->val[i] = 0.0;
      temp_y->val[i] = 0.0;
    }
    advection(u_x, u_x, u_y, u_z, temp_x, constant );
    advection(u_y, u_x, u_y, u_z, temp_y, constant );

    for( i=0;i<N_cells;i++){
      temp_x->val[i] = u_x->val[i] + constant.dt*temp_x->val[i]/(constant.dx*constant.dy);
      temp_y->val[i] = u_y->val[i] + constant.dt*temp_y->val[i]/(constant.dx*constant.dy);
    }

    test  = solve_BiCGSTAB(u_x, temp_x, constant,HELMHOLTZ);
    test  = solve_BiCGSTAB(u_y, temp_y, constant,HELMHOLTZ);
    // for(i=0;i<N_cells;i++)
    // if(u_x->val[i] > DELTA && u_x->bc[i]==NONE) printf("hello %lf %d\n", u_x->val[i], i);
    divergence(div,u_x,u_y,constant);
    test  = solve_BiCGSTAB(p, div, constant,POISSON);
    gradient(p,temp_x,temp_y, constant);
    for( i=0;i<N_cells;i++){
      u_x->val[i] = u_x->val[i] - constant.dt*temp_x->val[i]/(constant.dx*constant.dy);
      u_y->val[i] = u_y->val[i] - constant.dt*temp_y->val[i]/(constant.dx*constant.dy);
    }
    qq++;
    time += constant.dt;
    //set_bc(phi);
    // write_vtk(qq, domain, constant); 
printf("Time step %d, time %2.8f", qq, time);
printf("\r");
fflush(stdout);
  }
#ifdef _OPENMP
  printf("Time taken for the simulation: %lf", omp_get_wtime() - start_wall_time);
#endif

  write_line_data(domain,constant);
  return 0;
}

int solve_BiCGSTAB(Field *p, Field *div, Constant constant,Equation_type eqn) 
{
  int i ;
  double *r0_star, *sj, *rj, *pj, *pstar, *sstar ; 
  double *Temp, *Uj, *Var  ;
  double rhoj_Minus=0.0, alphaj=0.0, omegaj =0.0, rhoj, betaj, H1, H2 ;
  double norm, BICGEPS = 1.0E-12; 
  int BICG_ITER ;
  bool STOP = false ;
  int N =  p->N;
 
  double* b = div->val;

  r0_star = malloc(N*sizeof(double));
  sj = malloc(N*sizeof(double));
  rj = malloc(N*sizeof(double));
  pj = malloc(N*sizeof(double));
  pstar = malloc(N*sizeof(double));
  sstar = malloc(N*sizeof(double));
  Temp = malloc(N*sizeof(double));
  Uj = malloc(N*sizeof(double));
  Var = malloc(N*sizeof(double));

  // Start BICGSTAB iterations
  // set initial solution vector x_0 = (Uj, Vj)
  
  for(i = 0 ; i < N ; i++) {
Temp[i] = 0.0;
  }
  if(eqn == POISSON)
    laplacian(p, constant,Temp) ;
  else if(eqn == HELMHOLTZ)
    diffusion_implicit(p,constant, Temp);
  //Compute_AX(p, Temp, constant) ;
  // Initial vector r_0 = b - Ax_0, and r0* = r_0
  for(i = 0 ; i < N ; i++) {
    Uj[i] = p->val[i] ;
    rj[i] = b[i] - Temp[i] ; 
    r0_star[i] = rj[i] ;
  }
  BICG_ITER = 0 ; norm = 0.0 ;
  do {
    // compute rhoj = (r0, r0*)
    rhoj = 0.0 ;
    for(i = 0 ; i < N ; i++) 
      rhoj += rj[i]*r0_star[i] ; 
    if( sqrt(rhoj/((double)N)) < BICGEPS ) STOP = true ;
    else {
      if( BICG_ITER == 0 ) {
        for(i = 0 ; i < N ; i++) 
          pj[i] = rj[i]; // p0 = r0 
      } else {
        betaj = (rhoj/rhoj_Minus)*(alphaj/omegaj) ;
        for(i = 0 ; i < N ; i++) 
          pj[i] = rj[i] + betaj*(pj[i] - omegaj*Var[i]);
      }
      // Solve for Upstar, Vpstar from Ku* = u...., where K is the preconditioning matrix
      // No preconditioning
      for(i = 0 ; i < N ; i++)
        pstar[i] = pj[i] ;
      // compute vj = A*pstar
      for(i = 0 ; i < N ; i++) 
        p->val[i] = pstar[i] ;
      //    laplacian(p,constant,Temp) ;
      if(eqn == POISSON)
        laplacian(p, constant,Temp) ;
      else if(eqn == HELMHOLTZ)
        diffusion_implicit(p,constant, Temp);
      for(i = 0 ; i < N ; i++) 
        Var[i] = Temp[i] ;
      H1 = 0.0 ;
      for(i = 0 ; i < N ; i++) 
        H1 += Var[i]*r0_star[i] ;
      alphaj = rhoj/H1 ;
      // find sj
      for(i = 0 ; i < N ; i++) 
        sj[i] = rj[i] - alphaj*Var[i] ;
      // Solve for Upstar, Vpstar from Ku* = u..., where K is the preconditioning matrix
      // No preconditioning
      for(i = 0 ; i < N ; i++) 
        sstar[i] = sj[i] ;
      norm = 0.0 ;
      for(i = 0 ; i < N ; i++) 
        norm += sstar[i]*sstar[i] ;
      norm = sqrt(norm/((double)N)) ;
      if( norm < BICGEPS) {
        STOP = true ; //if ||s||_2 is small x_i = x_{i-1}+alphai*p_i
        for(i = 0 ; i < N ; i++) 
          Uj[i] += alphaj*pstar[i] ;
      } else {
        // compute t = As
        for(i = 0 ; i < N ; i++) 
          p->val[i] = sstar[i] ; 
        //laplacian(p,constant,Temp) ;
        if(eqn == POISSON)
          laplacian(p, constant,Temp) ;
        else if(eqn == HELMHOLTZ)
          diffusion_implicit(p,constant, Temp);
        H1 = H2 = 0.0 ;
        for(i =0 ; i < N ; i++) {
          H1 += Temp[i]*sj[i] ;
          H2 += Temp[i]*Temp[i] ;
        }	
        omegaj = H1/H2;
        // find xj 
        norm = 0.0 ;
        for(i = 0 ; i < N ; i++) {
          H1 = (alphaj*pstar[i] + omegaj*sstar[i]) ;
          Uj[i] += H1 ;
          norm += H1*H1 ;
        }
        norm = sqrt(norm/((double)N)) ;
        if(norm < BICGEPS) STOP = true ;
        // find rjplusone
        for(i = 0 ; i < N ; i++) 
          rj[i] = sj[i] - omegaj*Temp[i];
        rhoj_Minus = rhoj ;
      }
    }
    BICG_ITER++;
//    if(BICG_ITER%100 == 0) printf("%d \t %lf \n", BICG_ITER, norm );
  }while( (BICG_ITER < 10000) && (!STOP) ) ;
  for(i = 0 ; i < N ; i++)
    p->val[i] = Uj[i] ;
  
  free(r0_star);
  free(sj);
  free(rj);
  free(pj);
  free(pstar);
  free(sstar);
  free(Temp);
  free(Uj);
  free(Var);

  return 0;
}
/*
void Compute_AX(double * Temp){
  int i,l,m;
  double phi_w, phi_e, phi_n, phi_s;
  for(i=0;i<N_cells;i++){
    if(p_bc[i] == DIRICHLET){
      Temp[i] = p[i];
    }else{
      l= i%N_cells_x;
      m =(int) i/N_cells_x;
      int south = (m-1)*N_cells_x + l, north =(m+1)*N_cells_x + l,
          west =m*N_cells_x + (l-1), east = m*N_cells_x + (l+1);
        if(p_bc[south] == DIRICHLET)
          phi_s = 2.0*p[south] - p[i];
        else if(p_bc[south] == NEUMANN)
          phi_s = p[i];
        else phi_s = p[south];

      if(p_bc[north] == DIRICHLET)
        phi_n = 2.0*p[north] - p[i];
      else if(p_bc[north] == NEUMANN)
        phi_n = p[i];
      else phi_n = p[north];

      if(p_bc[east] == DIRICHLET)
        phi_e = 2.0*p[east] - p[i];
      else if(p_bc[east] == NEUMANN)
        phi_e = p[i];
      else phi_e = p[east];

      if(p_bc[west] == DIRICHLET)
        phi_w = 2.0*p[west] - p[i];
      else if(p_bc[west] == NEUMANN)
        phi_w = p[i];
      else phi_w = p[west];

      Temp[i] = -2.0*(dy/dx + dx/dy)*p[i] + phi_w*dy/dx + phi_e * dy/dx + phi_s * dx/dy + phi_n * dx/dy ;    
    }
  }
  return;
}
*/
static void write_line_data(Domain domain, Constant constant){
  char filename1[30];
  char filename2[30];
  sprintf(filename1, "line_data_u.dat");
  sprintf(filename2, "line_data_v.dat");
  FILE *fp1 = fopen(filename1, "w");
  FILE *fp2 = fopen(filename2, "w");
  int N_y = domain.u_x->N_y;
  int N_x = domain.u_y->N_x;
  int i,j,mid ;
  double loc=0.0;
  mid = (N_x-1)/2 + 2;
  for(j = 1; j <N_y-1; j ++){
    i = N_x*j+mid;
    fprintf(fp1, "%2.8lf \t %2.8lf \n", constant.dy*((j-1) + 0.5), domain.u_x->val[i]);
  }
  mid = (N_y-1)/2 + 2;
  for(j = 1; j <N_x-1; j ++){
    i = mid*N_x + j;
    fprintf(fp2, "%2.8lf \t %2.8lf \n", constant.dx*((j-1) + 0.5), domain.u_y->val[i]);
  }
}
static void write_vtk(int q, Domain domain, Constant constant)
{
  char filename[30]; 
  sprintf(filename, "output_%05d.vtk",q);
  FILE *fp = fopen(filename, "w");

  int Nx = domain.p->N_x+1;
  int Ny = domain.p->N_y+1;
  int Nz = 2; //N_cells_z+1;
  int N_cells_y = domain.p->N_y;
  int N_cells_x = domain.p->N_x;
  int N_cells_z = domain.p->N_z;
  int N_cells = N_cells_x * N_cells_y * N_cells_z;
  fprintf(fp,"# vtk DataFile Version 3.0\n");     
  fprintf(fp,"particle point data\n");           
  fprintf(fp,"ASCII\n");                         
  fprintf(fp,"DATASET STRUCTURED_GRID\n");       
  fprintf(fp,"DIMENSIONS %d %d %d\n",Nx,Ny,Nz);  
  fprintf(fp,"POINTS %d double\n",Nx*Ny*Nz);
  int l,m,n;
  for(n = 0; n<Nz; n++){
    for(m = 0; m<Ny; m++){
      for( l = 0; l<Nx ; l ++){
        fprintf(fp,"%2.8lf %2.8lf %2.8lf\n",l*constant.dx , m*constant.dy, n*1.0);
      }
    }
  }
  fprintf(fp,"CELL_DATA %d\n SCALARS pressure double 1\n LOOKUP_TABLE default\n",N_cells);  
  for( l = 0; l<N_cells_y ; l ++){
    for(m = 0; m<N_cells_x; m++){
      fprintf(fp,"%2.8lf\n",domain.p->val[l*N_cells_x + m]);
    }
  }
/*  fprintf(fp,"SCALARS phi double 1\n LOOKUP_TABLE default\n");  
  for( l = 0; l<N_cells_y ; l ++){
    for(m = 0; m<N_cells_x; m++){
      fprintf(fp,"%2.8lf\n",domain.phi->val[l*N_cells_x + m]);
    }
  }*/
  fprintf(fp,"SCALARS boundary int 1\n LOOKUP_TABLE default\n");  
  for( l = 0; l<N_cells_y ; l ++){
    for(m = 0; m<N_cells_x; m++){
      fprintf(fp,"%d\n",domain.p->bc[l*N_cells_x + m]);
    }
  }
  fprintf(fp,"VECTORS velocity double \n");  
  for( l = 0; l<domain.u_x->N_y ; l ++){
    for(m = 0; m<domain.u_x->N_x; m++){
      fprintf(fp,"%2.8lf %2.8lf %2.8lf \n",domain.u_x->val[l*domain.u_x->N_x + m], domain.u_y->val[l*domain.u_y->N_x + m], 0.0 );
    }
  }
  //  fprintf(fp,"%2.8lf %2.8lf %2.8lf\n", x.x, x.y, x.z);
  return;
}

static Field * allocate_field( int N_x, int  N_y, int N_z)
{
  Field * phi;
  phi      = malloc(sizeof(Field));
  phi->grid = CENTERED;
  phi->N_x = N_x;
  phi->N_y = N_y;
  phi->N_z = N_z;
  phi->N   = N_x*N_y*N_z;
  phi->val = malloc(phi->N * sizeof(double));
  phi->bc  = malloc(phi->N * sizeof(BC_type));
  return phi;
}
