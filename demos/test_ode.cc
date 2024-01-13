#define _USE_MATH_DEFINES
#include <cmath>
#include <nonlinfunc.h>
#include <ode.h>
#include <fstream>

using namespace ASC_ode;


class MassSpring : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f(0) = x(1);
    f(1) = -x(0);
    //std::cout << "here" << f << std::endl;
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,1) = 1;
    df(1,0) = -1;
  }
};


int main()
{
  double tend = 2*4*M_PI;
  int steps = 100;
  Vector<double> y(2);
  y(0) = 1;
  y(1) = 0;
  auto rhs = make_shared<MassSpring>();
  
  
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });
  
  /*SolveODE_EE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });
  
  SolveODE_CN(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });
  */

  
  //write output on txt file for plotting in ODE_plot.ipynb
  /*std::ofstream outf1("results_ode_ie.txt"), outf2("results_ode_ee.txt"), outf3("results_ode_cn.txt");

  SolveODE_IE(tend, steps, y, rhs,
              [&outf1](double t, VectorView<double> y) { outf1 << t << "  " << y(0) << " " << y(1) << std::endl; });
  outf1.close();

  SolveODE_EE(tend, steps, y, rhs,
              [&outf2](double t, VectorView<double> y) { outf2 << t << "  " << y(0) << " " << y(1) << std::endl; });
  outf2.close();
  
  SolveODE_CN(tend, steps, y, rhs,
              [&outf3](double t, VectorView<double> y) { outf3 << t << "  " << y(0) << " " << y(1) << std::endl; });
  outf3.close();*/

  
}
