#include <float.h>                          // required for DBL_EPSILON
#include <math.h>                           // required for fabs(), sqrt()

//                         Internally Defined Routines
static int Stopping_Rule(double x0, double x1, double tolerance);

#define sqrt5 2.236067977499789696

void Min_Search_Golden_Section( double (*f)(double), double* a, double *fa,
                                     double* b, double* fb, double tolerance)
{
   static const double lambda = 0.5 * (sqrt5 - 1.0);
   static const double mu = 0.5 * (3.0 - sqrt5);         // = 1 - lambda
   double x1;
   double x2;
   double fx1;
   double fx2;

   x1 = *b - lambda * (*b - *a);
   x2 = *a + lambda * (*b - *a);
   fx1 = f(x1);
   fx2 = f(x2);


   if (tolerance <= 0.0) tolerance = sqrt(DBL_EPSILON) * (*b - *a);

   while ( ! Stopping_Rule( *a, *b, tolerance) ) {
      if (fx1 > fx2) {
         *a = x1;
         *fa = fx1;
         if ( Stopping_Rule( *a, *b, tolerance) ) break;
         x1 = x2;
         fx1 = fx2;
         x2 = *b - mu * (*b - *a);
         fx2 = f(x2);
      } else {
         *b = x2;
         *fb = fx2;
         if ( Stopping_Rule( *a, *b, tolerance) ) break;
         x2 = x1;
         fx2 = fx1;
         x1 = *a + mu * (*b - *a);
         fx1 = f(x1);
      }
   }
   return;
}

static int Stopping_Rule(double x0, double x1, double tolerance)
{
   double xm = 0.5 * fabs( x1 + x0 );

   if ( xm <= 1.0 ) return ( fabs( x1 - x0 ) < tolerance ) ? 1 : 0;
   return ( fabs( x1 - x0 ) < tolerance * xm ) ? 1 : 0;
}
