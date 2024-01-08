#ifndef Newton_h
#define Newton_h

#include "nonlinfunc.h"
//#include <lapack_interface.h> // vorübergehend für Inverse

namespace ASC_ode
{

  void NewtonSolver (shared_ptr<NonlinearFunction> func, VectorView<double> x,
                     double tol = 1e-10, int maxsteps = 10,
                     std::function<void(int,double,VectorView<double>)> callback = nullptr)
  {
    Vector<> res(func->DimF());
    //std::cout << "res= " << res << std::endl;
    //Matrix<> fprime(func->DimF(), func->DimX());
    Matrix<> fprime(func->DimF(), func->DimX());
    //std::cout << "x=" << x << std::endl;
    for (int i = 0; i < maxsteps; i++)
      {
        func->Evaluate(x, res);
        //std::cout << "res= " << res << std::endl;
        // cout << "|res| = " << L2Norm(res) << endl;
        //std::cout << "fprime bf evalderiv= "<< fprime << std::endl;
        func->EvaluateDeriv(x, fprime);
        //std::cout << "fprime after evalderiv= "<< fprime << std::endl;
        fprime = fprime.inverse();

        x -= fprime*res;
        //std::cout << "newx= "<< x << std::endl;

        double err= L2Norm(res);
        //std::cout << "err= " << err << std::endl;
        if (callback)
          callback(i, err, x);
        if (err < tol) return;
      }

    throw std::domain_error("Newton did not converge");
  }

}

#endif
