#include<stdio.h>                                                               
#include<stdlib.h>                                                              
#include<math.h>                                                                
#include<string.h>  
#include<stdbool.h>                                                             
#include "fvmvof.h"

void advection( Field * phi, Field * u_x, 
    Field * u_y, Field * u_z,
    Field * tmp,
    Constant constant
    )
{
 
  double dx = constant.dx, dy = constant.dy, dz =constant.dz;
  int i, l, m;
  int N = phi->N;
  int N_x = phi->N_x;
  //int N_y = u_x->N_y;

  double u_x_e, u_x_w; //u_x_s, u_x_n, 
  double u_y_s, u_y_n; // u_y_e, u_y_w;
  double phi_s, phi_n , phi_e, phi_w;

  for(i = 0;i<N;i++){
    if(phi->bc[i] == NONE ){
      l= i%N_x;                                                           
      m =(int) i/N_x;

      int south = (m-1)*N_x + l, north =(m+1)*N_x + l,              
          west =m*N_x + (l-1), east = m*N_x + (l+1);

      u_x_e = u_x->val[east];
      u_x_w = u_x->val[west];
      u_y_s = u_y->val[south];
      u_y_n = u_y->val[north];
      phi_e = phi->val[east];
      phi_w = phi->val[west];
      phi_s = phi->val[south];
      phi_n = phi->val[north];

      if(phi->bc[east] != NONE )
        phi_e = 2.0*(phi->val[east]*abs(2-phi->bc[east]) + phi->val[i]*abs(1-phi->bc[east])) - phi->val[i];
      if(phi->bc[west] != NONE )
        phi_w = 2.0*(phi->val[west]*abs(2-phi->bc[west]) + phi->val[i]*abs(1-phi->bc[west])) - phi->val[i];
      if(phi->bc[north] != NONE )
        phi_n = 2.0*(phi->val[north]*abs(2-phi->bc[north]) + phi->val[i]*abs(1-phi->bc[north])) - phi->val[i];
      if(phi->bc[south] != NONE )
        phi_s = 2.0*(phi->val[south]*abs(2-phi->bc[south]) + phi->val[i]*abs(1-phi->bc[south])) - phi->val[i];
      if(u_x->bc[east] != NONE )
        u_x_e = 2.0*(u_x->val[east]*abs(2-u_x->bc[east]) + u_x->val[i]*abs(1-u_x->bc[east])) - u_x->val[i];
      if(u_x->bc[west] != NONE )
        u_x_w = 2.0*(u_x->val[west]*abs(2-u_x->bc[west]) + u_x->val[i]*abs(1-u_x->bc[west])) - u_x->val[i];
      if(u_y->bc[north] != NONE )
        u_y_n = 2.0*(u_y->val[north]*abs(2-u_y->bc[north]) + u_y->val[i]*abs(1-u_y->bc[north])) - u_y->val[i];
      if(u_y->bc[south] != NONE )
        u_y_s = 2.0*(u_y->val[south]*abs(2-u_y->bc[south]) + u_y->val[i]*abs(1-u_y->bc[south])) - u_y->val[i];

     tmp->val[i] += - 0.5* ((phi_e*u_x_e - phi_w*u_x_w)*dy +(phi_n*u_y_n -phi_s*u_y_s)*dx);
       
      // must be plus equal to 
    } else 
      tmp->val[i] = 0.0;

  }

//        printf("%lf %lf %lf %lf %lf\n", tmp[12*N_x + 12], phi->val[12*N_x + 12], dy, u_x->val[12*N_x + 12], u_y->val[12*N_x+12]);
return;
}
/*
m = N_x * mm
mpl = N_x *(mm+1)
  m_min = N_x*(mm-1)
for(l = 0;l<N_x;l++)
temp[m+ l] = u [ m+ l -1 ] + u [m+l +1] - u [mpl + 1]  
 Technique for vectorization */

void divergence(Field * div, Field * u_x, Field * u_y, Constant constant)
{
  int i,l,m;
  double dx = constant.dx, dy = constant.dy, dz = constant.dz;
  double dt = constant.dt;
  int N= div->N;
  int N_x = div->N_x;

  double u_x_e, u_x_w, u_y_n, u_y_s;
  for(i = 0;i<N;i++){                                                           
    if(u_x->bc[i] == NONE ){                                                    
      l= i%N_x;                                                                 
      m =(int) i/N_x; 
      int south = (m-1)*N_x + l, north =(m+1)*N_x + l,                          
          west =m*N_x + (l-1), east = m*N_x + (l+1);
      u_x_e = u_x->val[east];
      u_x_w = u_x->val[west];
      u_y_n = u_y->val[north];
      u_y_s = u_y->val[south];

      if(u_x->bc[east] != NONE )                                                
        u_x_e = 2.0*(u_x->val[east]*abs(2-u_x->bc[east]) + u_x->val[i]*abs(1-u_x->bc[east])) - u_x->val[i];
      if(u_x->bc[west] != NONE )                                                
        u_x_w = 2.0*(u_x->val[west]*abs(2-u_x->bc[west]) + u_x->val[i]*abs(1-u_x->bc[west])) - u_x->val[i];
      if(u_y->bc[north] != NONE )                                               
        u_y_n = 2.0*(u_y->val[north]*abs(2-u_y->bc[north]) + u_y->val[i]*abs(1-u_y->bc[north])) - u_y->val[i];
      if(u_y->bc[south] != NONE )                                               
        u_y_s = 2.0*(u_y->val[south]*abs(2-u_y->bc[south]) + u_y->val[i]*abs(1-u_y->bc[south])) - u_y->val[i];

      div->val[i] = 0.5*((u_x_e - u_x_w )*dy  + (u_y_n-u_y_s)*dx)/dt;

    }else{
      div->val[i] = 0.0;

    }

  }
  return;
}
void gradient(Field * p, Field * temp_x, Field * temp_y, Constant constant)
{
  int i,l,m;
  double dx = constant.dx, dy = constant.dy, dz = constant.dz;
  double dt = constant.dt;
  int N= p->N;
  int N_x = p->N_x;

  double p_e, p_w, p_n, p_s;
  for(i = 0;i<N;i++){                                                           
    if(p->bc[i] == NONE ){                                                    
      l= i%N_x;                                                                 
      m =(int) i/N_x; 
      int south = (m-1)*N_x + l, north =(m+1)*N_x + l,                          
          west =m*N_x + (l-1), east = m*N_x + (l+1);
      p_e = p->val[east];
      p_w = p->val[west];
      p_n = p->val[north];
      p_s = p->val[south];

      if(p->bc[east] != NONE )                                                
        p_e = 2.0*(p->val[east]*abs(2-p->bc[east]) + p->val[i]*abs(1-p->bc[east])) - p->val[i];
      if(p->bc[west] != NONE )                                                
        p_w = 2.0*(p->val[west]*abs(2-p->bc[west]) + p->val[i]*abs(1-p->bc[west])) - p->val[i];
      if(p->bc[north] != NONE )                                               
        p_n = 2.0*(p->val[north]*abs(2-p->bc[north]) + p->val[i]*abs(1-p->bc[north])) - p->val[i];
      if(p->bc[south] != NONE )                                               
        p_s = 2.0*(p->val[south]*abs(2-p->bc[south]) + p->val[i]*abs(1-p->bc[south])) - p->val[i];

      temp_x->val[i] = 0.5*(p_e - p_w )*dy  ;
      temp_y->val[i] =  0.5* (p_n-p_s)*dx;

    }else{
      temp_x->val[i] = 0.0;
      temp_y->val[i] = 0.0;

    }

  }
  return;
}
