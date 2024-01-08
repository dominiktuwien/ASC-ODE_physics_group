#define _USE_MATH_DEFINES
#include <cmath>
#include <nonlinfunc.h>
#include <ode.h>
#include <fstream>


using namespace ASC_ode;

class RHS : public NonlinearFunction
{
  size_t DimX() const override { return 1; }
  size_t DimF() const override { return 1; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = -x(0);
  }
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df(0,0) = -1;
  }
};


int main()
{
  double tend = 2*M_PI;
  int steps = 100;
  Vector<> x { 1, };
  Vector<> dx { 0, };
  auto rhs = make_shared<RHS>();
  auto mass = make_shared<IdentityFunction>(1);

  /*
  SolveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [](double t, VectorView<double> x) { std::cout << "t = " << t << ", x = " << x(0) << std::endl; }
                   );
  */

  //write output on txt file for plotting in ODE_plot.ipynb
  std::ofstream outf("test_ode_newmark.txt");
  SolveODE_Newmark(tend, steps, x, dx, rhs, mass,
                   [&outf](double t, VectorView<double> x) { outf << t << " " << x(0) << std::endl; }
                   );
  outf.close();
}
