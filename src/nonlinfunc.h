#ifndef NONLINFUNC_H
#define NONLINFUNC_H

#include <vector.h>
#include <Matrix-Neu-using_expressions.h>

#include <memory> // fix?



namespace ASC_ode
{
  using namespace ASC_bla;
  using std::shared_ptr;
  using std::make_shared;

  class NonlinearFunction
  {
  public:
    virtual ~NonlinearFunction() = default;
    virtual size_t DimX() const = 0;
    virtual size_t DimF() const = 0;
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const = 0;
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const = 0;
  };


  class IdentityFunction : public NonlinearFunction
  {
    size_t n;
  public:
    IdentityFunction (size_t _n) : n(_n) {  } 
    size_t DimX() const override { return n; }
    size_t DimF() const override { return n; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = x;
    }
    
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      //disabled, because it caused a wrong result when compared to an online newton solver for our problem
      df.Diag() = 1.0;
      //std::cout << "Diag - Matr = " << df << std::endl;
    }
  };



  class ConstantFunction : public NonlinearFunction
  {
    Vector<double> val;
  public:
    ConstantFunction (VectorView<double> _val) : val(_val) {  }
    void Set(VectorView<double> _val) { val = _val; }
    VectorView<double> Get() const { return val.View(); }
    size_t DimX() const override { return val.Size(); }
    size_t DimF() const override { return val.Size(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = val;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      //added as test to see if newton conv.
      //df.Diag() = 1.0;
    }
  };

  
  
  class SumFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa, fb;
    double faca, facb;
  public:
    SumFunction (shared_ptr<NonlinearFunction> _fa,
                 shared_ptr<NonlinearFunction> _fb,
                 double _faca, double _facb)
      : fa(_fa), fb(_fb), faca(_faca), facb(_facb) { } 
    
    size_t DimX() const override { return fa->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      
      //std::cout << "x bf sum eval= " << x << std::endl;
      fa->Evaluate(x, f);
      //std::cout << "f=" << f << std::endl;
      //std::cout << "x=" << x << std::endl;
      f *= faca;
      //std::cout << "f*faca" << f << std::endl;
      Vector<> tmp(DimF());
      fb->Evaluate(x, tmp);
      f += facb*tmp;
      //std::cout << "f += facb*tmp" << f << std::endl;
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa->EvaluateDeriv(x, df);
      df *= faca;
      //std::cout << "faca*df = " << df << std::endl;
      Matrix<double> tmp(DimF(), DimX());      
      fb->EvaluateDeriv(x, tmp);
      //std::cout << "df bf in Sumfunc: " << df << std::endl;
      df += facb*tmp;
      //df = tmp;
      //std::cout << "df += facb*tmp " << df << std::endl;
    }
  };


  inline auto operator- (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<SumFunction>(fa, fb, 1, -1);
  }

  inline auto operator+ (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<SumFunction>(fa, fb, 1, 1);
  }

  
  class ScaleFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa;
    double fac;
  public:
    ScaleFunction (shared_ptr<NonlinearFunction> _fa,
                   double _fac)
      : fa(_fa), fac(_fac) { } 
    
    size_t DimX() const override { return fa->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      fa->Evaluate(x, f);
      f *= fac;

    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      fa->EvaluateDeriv(x, df);
      df *= fac;
    }
  };

  inline auto operator* (double a, shared_ptr<NonlinearFunction> f)
  {
    return make_shared<ScaleFunction>(f, a);
  }




  // fa(fb)
  class ComposeFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa, fb;
  public:
    ComposeFunction (shared_ptr<NonlinearFunction> _fa,
                     shared_ptr<NonlinearFunction> _fb)
      : fa(_fa), fb(_fb) { } 
    
    size_t DimX() const override { return fb->DimX(); }
    size_t DimF() const override { return fa->DimF(); }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      Vector<> tmp(fb->DimF());
      fb->Evaluate (x, tmp);
      fa->Evaluate (tmp, f);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      Vector<> tmp(fb->DimF());
      fb->Evaluate (x, tmp);
      
      Matrix<> jaca(fa->DimF(), fa->DimX());
      Matrix<> jacb(fb->DimF(), fb->DimX());
      // Matrix<double> jaca(fa->DimF(), fa->DimX());
      // Matrix<double> jacb(fb->DimF(), fb->DimX());
      
      fb->EvaluateDeriv(x, jacb);
      fa->EvaluateDeriv(tmp, jaca);

      df = jaca*jacb;
    }
  };
  
  
  inline auto Compose (shared_ptr<NonlinearFunction> fa, shared_ptr<NonlinearFunction> fb)
  {
    return make_shared<ComposeFunction> (fa, fb);
  }
  
  class EmbedFunction : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> fa;
    size_t firstx, dimx, firstf, dimf;
    size_t nextx, nextf;
  public:
    EmbedFunction (shared_ptr<NonlinearFunction> _fa,
                   size_t _firstx, size_t _dimx,
                   size_t _firstf, size_t _dimf)
      : fa(_fa),
        firstx(_firstx), dimx(_dimx), firstf(_firstf), dimf(_dimf),
        nextx(_firstx+_fa->DimX()), nextf(_firstf+_fa->DimF())
    { }
    
    size_t DimX() const override { return dimx; }
    size_t DimF() const override { return dimf; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      fa->Evaluate(x.Range(firstx, nextx), f.Range(firstf, nextf));
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0;
      fa->EvaluateDeriv(x.Range(firstx, nextx),
                        df.Rows(firstf, nextf).Columns(firstx, nextx));      
    }
  };

  
  class Projector : public NonlinearFunction
  {
    size_t size, first, next;
  public:
    Projector (size_t _size, 
               size_t _first, size_t _next)
      : size(_size), first(_first), next(_next) { }
    
    size_t DimX() const override { return size; }
    size_t DimF() const override { return size; }
    void Evaluate (VectorView<double> x, VectorView<double> f) const override
    {
      f = 0.0;
      f.Range(first, next) = x.Range(first, next);
    }
    void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
    {
      df = 0.0;
      df.Diag().Range(first, next) = 1;
    }
  };


  class MultipleFunc : public NonlinearFunction
  {
    shared_ptr<NonlinearFunction> func;
    size_t num, fdimx, fdimf;
  public:
    MultipleFunc (shared_ptr<NonlinearFunction> _func, int _num)
      : func(_func), num(_num)
    {
      fdimx = func->DimX();
      fdimf = func->DimF();
    }

    virtual size_t DimX() const { return num * fdimx; } 
    virtual size_t DimF() const { return num * fdimf; }
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const
    {
      for (size_t i = 0; i < num; i++)
        func->Evaluate(x.Range(i*fdimx, (i+1)*fdimx),
                       f.Range(i*fdimf, (i+1)*fdimf));
    }
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
    {
      df = 0.0;
      for (size_t i = 0; i < num; i++)
        func->EvaluateDeriv(x.Range(i*fdimx, (i+1)*fdimx),
                            df.Rows(i*fdimf, (i+1)*fdimf).Columns(i*fdimx, (i+1)*fdimx));
    }
  };


  class MatVecFunc : public NonlinearFunction
  {
    Matrix<> a;
    size_t n;
  public:
    MatVecFunc (Matrix<> _a, size_t _n)
      : a(_a), n(_n) { }

    virtual size_t DimX() const { return n*a.Get_height(); } 
    virtual size_t DimF() const { return n*a.Get_width(); }
    virtual void Evaluate (VectorView<double> x, VectorView<double> f) const
    {
      MatrixView<> mx(a.Get_width(), n, x.Data(), n);
      MatrixView<> mf(a.Get_height(), n, f.Data(), n);
      mf = a * mx;
    }
    virtual void EvaluateDeriv (VectorView<double> x, MatrixView<double> df) const
    {
      df = 0.0;
      for (size_t i = 0; i < a.Get_height(); i++)
        for (size_t j = 0; j < a.Get_width(); j++)
          df.Rows(i*n, (i+1)*n).Columns(j*n, (j+1)*n).Diag() = a(i,j);
    }
  };

  
}

#endif
