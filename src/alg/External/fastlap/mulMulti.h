#ifndef MUL_MULTI_H
#define MUL_MULTI_H

   #include <stdio.h>
   #include <math.h>
   #include "mulStruct.h"
   #include "mulGlobal.h"
   
   /*cdel extern int nr_mulLocal_calls, nr_mulMulti_calls;*/
   
   /* 
      Used various places.  Returns number of coefficients in the multipole 
      expansion. 
   */
   int multerms(int order);
   
   /*
     returns number of cos(m*phi)-weighted terms in the real (not cmpx) multi exp
   */
   int costerms(int order);
   /*
     returns number of sin(m*phi)-weighted terms in the real (not cmpx) multi exp
   */
   int sinterms(int order);
   
   
   /*
     takes two sets of cartesian absolute coordinates; finds rel. spherical coor.
   */
   void xyz2sphere(double x, double y, double z, double x0, double y0, double z0, double *rho, double *cosA, double *beta);
   
   /*
     gives the linear index into vector from n and m used by all routines dealing
     with moments (cosine parts) and Leg. function evals, e.g. Mn^m and Pn^m(cosA)
     used for all cases except for the sine part of arrays (use sindex() instead)
     assumed entry order: (n,m) = (0,0) (1,0) (1,1) (2,0) (2,1) (2,2) (3,0)...
   */
   int index(int n, int m);
   
   /*
     gives the linear index into vector from n and m used by all routines dealing
     with moments (sine parts), e.g. Mn^m 
     assumes an array with all m = 0 (Mn^0) entries omitted to save space
     assumed entry order: (n,m) = (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1)...
   */
   int sindex(int n, int m, int cterms);
   /*
     returns i = sqrt(-1) to the power of the argument
   */
   double iPwr(int e);
   
   /*
     returns factorial of the argument (x!)
   */
   double fact(int x);
   
   /*
     produces factorial factor array for mulMulti2P
   */
   void evalFactFac(double** array, int order);
   
   /*
     Allocates space for temporary vectors.
   */
   void mulMultiAlloc(int maxsngs, int order, int depth);
   
   /*
    * evalLegendre returns a vector of Legendre function evaluations of the 
    * form Pn^m(cosA) n and m have maximum value = order.
    * Vector entries correspond to (n,m) = (0,0) (1,0) (1,1) (2,0) (2,1)...
    */
   void evalLegendre(double cosA, double* vector, int order);

   /* 
     Used for the upward pass. 
   */
   double **mulQ2Multi(snglrty** sngs, int numsngs, double x, double y, double z, int didthis, int order, double** mat);
     
   
   double **mulMulti2Multi(double x, double y, double z, double xp, double yp, double zp, int order, double** mat);
   
   /* 
     builds multipole evaluation matrix; used only for fake downward pass
     x,y,z is the multipole coord, fpts[]->x,y,z are the evaluation points.
   */
   double **mulMulti2P(double x, double y, double z, fieldpt** fpts, int numfpts, int order);  

#endif