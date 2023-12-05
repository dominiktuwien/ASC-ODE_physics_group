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

  
  /*
  SolveODE_RK(tend, steps, Gauss2a, Gauss2b, Gauss2c, y, rhs, 
              [](double t, VectorView<double> y) { cout << t << "  " << y(0) << " " << y(1) << endl; });
  */
}
