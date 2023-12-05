#include <nonlinfunc.h>
#include <ode.h>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};


Matrix<> Gauss2a { { 0.25, 0.25 - sqrt(3)/6 }, { 0.25 + sqrt(3)/6, 0.25 } };
Vector<> Gauss2b { 0.5, 0.5 };
Vector<> Gauss2c { 0.5 - sqrt(3)/6, 0.5 + sqrt(3)/6 };

Vector<> Gauss3c { 0.5 - sqrt(15)/10, 0.5, 0.5+sqrt(15)/10 };


// codes from Numerical Recipes, https://numerical.recipes/book.html

// Gauss integration on [0,1]
void GaussLegendre(VectorView<> x, VectorView<> w)
// Given the lower and upper limits of integration x1 and x2, this routine returns arrays x[0..n-1]
// and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-Legendre n-point
// quadrature formula.
  {
    double x1 = 0;
    double x2 = 1;
    const double EPS=1.0e-14;  // EPS is the relative precision.
    double z1,z,xm,xl,pp,p3,p2,p1;
    int n=x.Size();
    int m=(n+1)/2;  // The roots are symmetric in the interval, so
    xm=0.5*(x2+x1); // we only have to find half of them.
    xl=0.5*(x2-x1);
    for (int i=0;i<m;i++) {  // Loop over the desired roots.
      z=cos(3.141592654*(i+0.75)/(n+0.5));
      // Starting with this approximation to the ith root, we enter the main loop of refinement
      // by Newton’s method.
      do {
        p1=1.0;
        p2=0.0;
        for (int j=0;j<n;j++) { // Loop up the recurrence relation to get the
            p3=p2;  // Legendre polynomial evaluated at z.
          p2=p1;
          p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
        }
        // p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
        // by a standard relation involving also p2, the polynomial of one lower order.
        pp=n*(z*p1-p2)/(z*z-1.0);
        z1=z;
        z=z1-p1/pp;   // Newton’s method.
      } while (abs(z-z1) > EPS);
      x[i]=xm-xl*z;      // Scale the root to the desired interval,
      x[n-1-i]=xm+xl*z;  //  and put in its symmetric counterpart.
      w[i]=2.0*xl/((1.0-z*z)*pp*pp);  // Compute the weight
      w[n-1-i]=w[i];     // and its symmetric counterpart.
    }
  }


void GaussJacobi (VectorView<> x, VectorView<> w, const double alf, const double bet)
// Given alf and bet, the parameters ˛ and ˇ of the Jacobi polynomials, this routine returns
// arrays x[0..n-1] and w[0..n-1] containing the abscissas and weights of the n-point GaussJacobi quadrature formula. The largest abscissa is returned in x[0], the smallest in x[n-1].
{
  const int MAXIT=10;
  const double EPS=1.0e-14; // EPS is the relative precision.
  int i,its,j;
  double alfbet,an,bn,r1,r2,r3;
  double a,b,c,p1,p2,p3,pp,temp,z,z1;
  int n=x.Size();
  for (i=0;i<n;i++) { // Loop over the desired roots.
    if (i == 0) {  // Initial guess for the largest root.
      an=alf;
      bn=bet/n;
      r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
      r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
      z=1.0-r1/r2;
    } else if (i == 1) { // Initial guess for the second largest root.
      r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
      r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
      r3=1.0+0.012*bet*(1.0+0.25*abs(alf))/n;
      z -= (1.0-z)*r1*r2*r3;
    } else if (i == 2) { // Initial guess for the third largest root.
      r1=(1.67+0.28*alf)/(1.0+0.37*alf);
      r2=1.0+0.22*(n-8.0)/n;
      r3=1.0+8.0*bet/((6.28+bet)*n*n);
      z -= (x[0]-z)*r1*r2*r3;
    } else if (i == n-2) { // Initial guess for the second smallest root.
      r1=(1.0+0.235*bet)/(0.766+0.119*bet);
      r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
      r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
      z += (z-x[n-4])*r1*r2*r3;
    } else if (i == n-1) { // Initial guess for the smallest root.
      r1=(1.0+0.37*bet)/(1.67+0.28*bet);
      r2=1.0/(1.0+0.22*(n-8.0)/n);
      r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
      z += (z-x[n-3])*r1*r2*r3;
    } else { // Initial guess for the other roots.
      z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
    }
    alfbet=alf+bet;
    for (its=1;its<=MAXIT;its++) { // Refinement by Newton’s method.
      temp=2.0+alfbet; // Start the recurrence with P0 and P1 to avoid
      // a division by zero when ˛ C ˇ D 0 or 1.
      p1=(alf-bet+temp*z)/2.0;
      p2=1.0;
      for (j=2;j<=n;j++) { // Loop up the recurrence relation to get the
        p3=p2;    // Jacobi polynomial evaluated at z.
        p2=p1;
        temp=2*j+alfbet;
        a=2*j*(j+alfbet)*(temp-2.0);
        b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
        c=2.0*(j-1+alf)*(j-1+bet)*temp;
        p1=(b*p2-c*p3)/a;
      }
      pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
      // p1 is now the desired Jacobi polynomial. We next compute pp, its derivative, by
      //  a standard relation involving also p2, the polynomial of one lower order.
      z1=z;
      z=z1-p1/pp; // Newton’s formula.
      if (abs(z-z1) <= EPS) break;
    }
    if (its > MAXIT) throw("too many iterations in gaujac");
    x[i]=z;    // Store the root and the weight.
    w[i]=exp(std::lgamma(alf+n)+std::lgamma(bet+n)-std::lgamma(n+1.0)-
             std::lgamma(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
  }
}




auto ComputeABfromC (const Vector<> & c)
{
  int s = c.Size();
  Matrix M(s, s);
  Vector tmp(s);
  
  for (int i = 0; i < s; i++)
    for (int j = 0; j < s; j++)
      M(i,j) = std::pow(c(j), i);
  CalcInverse(M);
  
  for (int i = 0; i < s; i++)
    tmp(i) = 1.0 / (i+1);

  Vector b = M * tmp;
  cout << "b = " << b << endl;

  Matrix a(s,s);

  for (int j = 0; j < s; j++)
    {
      for (int i = 0; i < s; i++)
        tmp(i) = std::pow(c(j),i+1) / (i+1);
      a.Row(j) = M * tmp;
    }
  /*
  cout << "b = " << b << endl;
  cout << "a = " << a << endl;
  */
  return tuple { a, b };
}
  

void GaussRadau (VectorView<> x, VectorView<> w)
{
  GaussJacobi (x.Range(0, x.Size()-1),
               w.Range(0, w.Size()-1), 1, 0);
  for (int i = 0; i < x.Size()-1; i++)
    {
      x(i) = 0.5*(x(i)+1);
      w(i) *= 0.5;
    }
  x(x.Size()-1) = 1.0;
  double sum = 0;
  for (int i = 0; i < x.Size()-1; i++)
    sum += w(i);
  w(x.Size()-1) = 1-sum;
}

int main()
{
  double tend = 4*M_PI;
  int steps = 200;
  Vector<> y { 1, 0 };
  auto rhs = make_shared<MassSpring>();

  /*
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });
  */

  auto [a,b] = ComputeABfromC (Gauss3c);
  SolveODE_RK(tend, steps, a, b, Gauss3c, y, rhs, 
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });


  cout << "Gauss3c = " << Gauss3c << endl;
  cout << "weights = " << b << endl;
  GaussLegendre (Gauss3c, b);
  cout << "with generic function, c = " << Gauss3c << ", weights = " << b << endl;


  Vector<> Radau(3), RadauWeight(3);
  GaussRadau (Radau, RadauWeight);
  // not sure about weights, comput them via ComputeABfromC
  cout << "Radau = " << Radau << ", weight = " << RadauWeight <<  endl;
  
  /*
  SolveODE_RK(tend, steps, Gauss2a, Gauss2b, Gauss2c, y, rhs, 
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });
  */
}
