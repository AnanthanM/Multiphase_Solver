#include<stdio.h>                                                               
#include<stdlib.h>                                                              
#include<math.h>                                                                
#include<string.h>  
#include<stdbool.h>                                                             
#include "fvmvof.h"

void set_ghosts(Domain domain)
{ 
  int i,l,m;
  Field *p = domain.p;
  Field *u_x = domain.u_x;
  Field *u_y = domain.u_y;
  Field *u_z = domain.u_z;
 // Field *phi = domain.phi;

int N_cells = p->N_x * p->N_y;
int N_cells_x = p->N_x;
int N_cells_y = p->N_y;

  for(i=0;i<N_cells;i++){         
    l = i%N_cells_x;               
    m = (int) i/N_cells_x;         
    if(l==0 || l == N_cells_x-1 || m==0 || m== N_cells_y-1){ 
      //bc[i]=AMBIENT;             
      p->bc[i]   = NEUMANN;         
      u_x->bc[i] = DIRICHLET;       
      u_y->bc[i] = DIRICHLET;       
      u_z->bc[i] = DIRICHLET;       
//      phi->bc[i] = DIRICHLET;       
    }else{                         
      //bc[i] = NONE;              
      p->bc[i]   = NONE;            
     // phi->bc[i]   = NONE;            
      u_x->bc[i] = NONE;            
      u_y->bc[i] = NONE;            
      u_z->bc[i] = NONE;            
//                           
/*      if(l>N_cells_x/4  && l<3*N_cells_x/4  && m>N_cells_y/4 && m<3*N_cells_y/4){
      p->bc[i]   = NEUMANN;              
      u_x->bc[i] = DIRICHLET;            
      u_y->bc[i] = DIRICHLET;            
      u_z->bc[i] = DIRICHLET;            
      //   bc[i]=WALL;    
   }                                     */                   
    }                                                            
  }                                                              
  return;                                                        
}         
void set_bc(Field * phi )                                              
{                                                                  
  int i;                                                           
int N_cells = phi->N_x * phi->N_y;
int N_cells_x = phi->N_x;
int N_cells_y = phi->N_y;
  for(i=0;i<N_cells;i++){                                          
    int l = i%N_cells_x, m = (int) i/N_cells_x; // x,y coordinates 
    if(phi->bc[i] == DIRICHLET){                                        
      if(l==0)                phi->val[i] = phi->bc_val[XMIN];     
      else if(l==N_cells_x-1) phi->val[i] = phi->bc_val[XMAX];     
      else if(m==0)           phi->val[i] = phi->bc_val[YMIN];     
      else if(m==N_cells_y-1) phi->val[i] = phi->bc_val[YMAX];     
      else                    phi->val[i] = phi->bc_val[SOLID];    
    }                                                              
  }                                                                
  return;                                                                       
}        
