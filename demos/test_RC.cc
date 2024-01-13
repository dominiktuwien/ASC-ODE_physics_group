#define _USE_MATH_DEFINES
#include <cmath>
#include <nonlinfunc.h>
#include <ode.h>
#include <fstream>

using namespace ASC_ode;

// TODO: Newton does not converge, results of explicit Euler don't seem right -> fix


// ODE: U_C + RC * U_C'(t) = U_0(t), with U_0(t) = cos(100 * pi * t)
// -> U_C'(t) = (cos(100*pi*t) - U_C(t)) / RC
// autonomous form: new variable y2(t) := t where y2'(t)=1
// -> U_C'(t) = (cos(100*pi*y2(t)) - U_C(t)) / RC
class ElectricalCircuit : public NonlinearFunction
{
  size_t DimX() const override { return 2; }
  size_t DimF() const override { return 2; }

  double R = 1.0; // 1 or 100
  double C = 1.0; // 1 or 10^-6
  
  void Evaluate (VectorView<double> x, VectorView<double> f) const override
  { 
    f(0) = (cos(100 * M_PI * x(0)) - x(1)) / (R * C); // U_C'(t) = (cos(100*pi*y2(t)) - U_C(t)) / RC
    f(1) = 1; // y2'(t) = 1
  }
  
  void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
    df = 0.0;
    df(0,0) = -1 / (R*C);
    df(0,1) = - (100 * M_PI * sin(100 * M_PI * x(0))) / (R*C);
  }
};


int main()
{
  double tend = 2*M_PI;
  int steps = 400;
  Vector<double> y(2);
  y(0) = 0; // starting values
  y(1) = 0;
  auto rhs = make_shared<ElectricalCircuit>();
  
  
  SolveODE_IE(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << t << " " << y(1) << std::endl; });
  
  /*
  std::ofstream outf("RC_EE.txt");
  SolveODE_EE(tend, steps, y, rhs,
              [&outf](double t, VectorView<double> y) { outf << t << "  " << y(0) << " " << y(1) << std::endl; });
  outf.close();
  */
  
  /*
  SolveODE_CN(tend, steps, y, rhs,
              [](double t, VectorView<double> y) { std::cout << t << "  " << y(0) << " " << y(1) << std::endl; });
  */
}