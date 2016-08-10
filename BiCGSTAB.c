#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

#include "fvm.h"

// Solves Ax = b using BiCGSTAB Alogorithm in Alogorithm.txt
// Program takes Ax vector and b  as inputs 


int solve_Pressure_Poisson_BiCGSTAB(Field *x, Field *rhs, Constant constant, Domain domain) 
{
  int    i ;
  double *r0_cap, *sj, *rj, *pj, *pcap, *scap ; 
  double *Ax_vector, *h, *vj  ;
  double rhoj_Minus=0.0, alphaj=0.0, omegaj =0.0, rhoj, betaj, H1, H2 ;
  double norm, BICGEPS = 1.0E-12; 
  int    BICG_ITER ;
  bool   STOP = false ;
  int    N =  x->N;
 
  double * b = rhs->val;

  r0_cap     = malloc(N*sizeof(double));
  sj         = malloc(N*sizeof(double));
  rj         = malloc(N*sizeof(double));
  pj         = malloc(N*sizeof(double));
  pcap       = malloc(N*sizeof(double));
  scap       = malloc(N*sizeof(double));
  Ax_vector  = malloc(N*sizeof(double));
  h          = malloc(N*sizeof(double));
  vj         = malloc(N*sizeof(double));

  
  for(i = 0 ; i < N ; i++) 
  {
    Ax_vector[i] = 0.0;
  }

  Pressure_Poisson(x,constant,Ax_vector,domain) ;

  for(i = 0 ; i < N ; i++) 
  {
     h[i] = x->val[i] ;
     rj[i] = b[i] - Ax_vector[i] ;            // Step 1
     r0_cap[i] = rj[i] ;                      // Step 2
  }

  BICG_ITER = 0 ; norm = 0.0 ;

  do 
  {
    rhoj = 0.0 ;
    
    for(i = 0 ; i < N ; i++) 
        rhoj += rj[i]*r0_cap[i] ;             // Step 5

    if( sqrt(rhoj/((double)N)) < BICGEPS ) 
        STOP = true ;
    else
    {
        if( BICG_ITER == 0 ) 
        {   
          for(i = 0 ; i < N ; i++) 
            pj[i] = rj[i];                    // Step 7 for j = 1 
        } 
        else 
        {
          betaj = (rhoj/rhoj_Minus)*(alphaj/omegaj) ;  // Step 6 
          for(i = 0 ; i < N ; i++) 
            pj[i] = rj[i] + betaj*(pj[i] - omegaj*vj[i]);  // Step 7 
        } 
        
        for(i = 0 ; i < N ; i++)
          pcap[i] = pj[i] ;
        for(i = 0 ; i < N ; i++) 
          x->val[i] = pcap[i] ;
        
        Pressure_Poisson(x,constant,Ax_vector,domain) ;
        
        for(i = 0 ; i < N ; i++) 
          vj[i] = Ax_vector[i] ;                     // Step 8
        
        H1 = 0.0 ;
        
        for(i = 0 ; i < N ; i++) 
          H1 += vj[i]*r0_cap[i] ;
        
        alphaj = rhoj/H1 ;                           // Step 9 
        
        for(i = 0 ; i < N ; i++) 
          sj[i] = rj[i] - alphaj*vj[i] ;             // Step 12
      
        for(i = 0 ; i < N ; i++) 
          scap[i] = sj[i] ;
          
        norm = 0.0 ;
      
        for(i = 0 ; i < N ; i++) 
          norm += scap[i]*scap[i] ;
      
        norm = sqrt(norm/((double)N)) ;
        
        if( norm < BICGEPS) 
        {
          STOP = true ; 
          for(i = 0 ; i < N ; i++) 
              h[i] += alphaj*pcap[i] ;               // Step 10
        } 
        else 
        {
           for(i = 0 ; i < N ; i++) 
              x->val[i] = scap[i] ;                          //
                                                             // Together Step 13
           Pressure_Poisson(x,constant,Ax_vector,domain) ;   // Now Ax_vector is t
           
           H1 = H2 = 0.0 ;
        
           for(i =0 ; i < N ; i++) 
           {
               H1 += Ax_vector[i]*sj[i] ;
               H2 += Ax_vector[i]*Ax_vector[i] ;
           }	
           omegaj = H1/H2;                                // Step 14
           
           norm = 0.0 ;
           for(i = 0 ; i < N ; i++) 
           { 
               H1 = (alphaj*pcap[i] + omegaj*scap[i]) ;
               h[i] += H1 ;                                // Step 15                     
               norm += H1*H1 ;
           }  
           norm = sqrt(norm/((double)N)) ;
           if(norm < BICGEPS) 
              STOP = true ;
          
           for(i = 0 ; i < N ; i++) 
              rj[i] = sj[i] - omegaj*Ax_vector[i] ;          // Step 17
            
           rhoj_Minus = rhoj ;
        } 
    }
 
    BICG_ITER++;

  } while( (BICG_ITER < 10000) && (!STOP) ) ;

  for(i = 0 ; i < N ; i++)
    x->val[i] = h[i] ;
  
  free(r0_cap);
  free(sj);
  free(rj);
  free(pj);
  free(pcap);
  free(scap);
  free(Ax_vector);
  free(h);
  free(vj);

  return 0;
}
