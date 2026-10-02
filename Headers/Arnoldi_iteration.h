#include "NLA_FUNCTIONS.h"
#include <iostream>

class Arnoldi_iteration
{
public:
      Matrix Q;
      Matrix H;
      Arnoldi_iteration(Matrix &A, Matrix &q, int k);
};

class Polynomial
{
private:
      vector<double> coeff;

public:
      Polynomial(const vector<double> &coefficients);
      Polynomial();

      void increasePower();

      void addPolynomial(
          const Polynomial &p,
          double scalar);
      int degree() const;
      void scale(double scalar);

      double evaluate(double x) const;

      vector<double> getCoefficients() const;
};

class Arnoldi_Krylov_Polynomial
{
public:
      Matrix Q;
      Matrix X;
      Matrix H;

      vector<Polynomial> candidatePolynomials;
      vector<double> coefficients;
      vector<double> residuals;
      Polynomial finalPolynomial;
      int effectiveDegree;

      Arnoldi_Krylov_Polynomial(const Matrix &data);
};

