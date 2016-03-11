/**************************************************************************
 
  The dune-fem module is a module of DUNE (see www.dune-project.org).
  It is based on the dune-grid interface library 
  extending the grid interface by a number of discretization algorithms
  for solving non-linear systems of partial differential equations.

  Copyright (C) 2003 - 2014 Robert Kloefkorn
  Copyright (C) 2003 - 2010 Mario Ohlberger 
  Copyright (C) 2004 - 2014 Andreas Dedner
  Copyright (C) 2005        Adrian Burri
  Copyright (C) 2005 - 2014 Mirko Kraenkel
  Copyright (C) 2006 - 2014 Christoph Gersbacher
  Copyright (C) 2006 - 2014 Martin Nolte
  Copyright (C) 2011 - 2014 Tobias Malkmus
  Copyright (C) 2012 - 2014 Stefan Girke
  Copyright (C) 2013 - 2014 Claus-Justus Heine
  Copyright (C) 2013 - 2014 Janick Gerstenberger
  Copyright (C) 2013        Sven Kaulman
  Copyright (C) 2013        Tom Ranner


  The dune-fem module is free software; you can redistribute it and/or 
  modify it under the terms of the GNU General Public License as 
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.

  The dune-fem module is distributed in the hope that it will be useful, 
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
**************************************************************************/
#ifndef POISSON_PROBLEMS_HH
#define POISSON_PROBLEMS_HH

#include <cassert>
#include <cmath>

#include "probleminterface.hh"

template <class FunctionSpace>
class BulkProblem : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double Xx = x[0], Yy = x[1];
    const double beta = 1.0;
#warning using beta = 1 here
    phi[ 0 ] =  beta*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(3. + 4.*Xx - 4.*pow(Xx,2) + 4.*Yy - 4.*pow(Yy,2));
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    const double beta = 1.0;
#warning using beta = 1 here
    phi[ 0 ] = beta * exp( - x[0] * (x[0] - 1.0) - x[1] * (x[1] - 1.0)  );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    const double beta = 1.0;
#warning using beta = 1 here

    const double Xx = x[0], Yy = x[1];
    ret[0][0] = beta*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. - 2.*Xx);
    ret[0][1] = beta*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. - 2.*Yy);
    ret[0][2] = 0.0;
  }
  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(1);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
		 RangeType& value) const
  {
    value = RangeType( 0 );
  }
};

template <class FunctionSpace>
class SurfaceProblem : public ProblemInterface < FunctionSpace >
{
  typedef ProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double alpha = 1.0;
    const double beta = 1.0;
#warning using alpha = 1 here
#warning using beta = 1 here

    const double Xx = x[0], Yy = x[1], Zz = x[2];
    phi[ 0 ] = beta*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1.*Xx - 2.*pow(Xx,2) + 1.*Yy - 2.*pow(Yy,2)) + alpha*pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(3. + 4.*pow(Xx,5) - 4.*pow(Xx,6) + 4.*pow(Yy,5) - 4.*pow(Yy,6) + pow(Yy,4)*(13. - 4.*pow(Zz,2)) + pow(Xx,4)*(13. + 4.*Yy - 12.*pow(Yy,2) - 4.*pow(Zz,2)) + Yy*(8. - 2.*pow(Zz,2)) + pow(Yy,3)*(-10. + 4.*pow(Zz,2)) + pow(Xx,3)*(-10. - 2.*Yy + 8.*pow(Yy,2) + 4.*pow(Zz,2)) + pow(Yy,2)*(-14. + 5.*pow(Zz,2)) + pow(Xx,2)*(-14. + 8.*pow(Yy,3) - 12.*pow(Yy,4) + 5.*pow(Zz,2) + pow(Yy,2)*(26. - 8.*pow(Zz,2)) + Yy*(-10. + 4.*pow(Zz,2))) + Xx*(8. - 2.*pow(Yy,3) + 4.*pow(Yy,4) - 2.*pow(Zz,2) + Yy*(4. - 2.*pow(Zz,2)) + pow(Yy,2)*(-10. + 4.*pow(Zz,2)))) + pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(4. + 19.*Xx - 34.*pow(Xx,2) - 12.*pow(Xx,7) + 8.*pow(Xx,8) + 19.*Yy - 34.*pow(Yy,2) - 56.*pow(Yy,3) + 60.*pow(Yy,4) + 45.*pow(Yy,5) + 12.*Xx*pow(Yy,5) - 36.*pow(Xx,2)*pow(Yy,5) - 38.*pow(Yy,6) - 12.*Xx*pow(Yy,6) + 32.*pow(Xx,2)*pow(Yy,6) - 12.*pow(Yy,7) + 8.*pow(Yy,8) - 2.*Xx*pow(Zz,2) + 8.*pow(Xx,2)*pow(Zz,2) - 2.*Yy*pow(Zz,2) + 8.*pow(Yy,2)*pow(Zz,2) + 21.*pow(Yy,3)*pow(Zz,2) - 22.*pow(Yy,4)*pow(Zz,2) - 12.*pow(Yy,5)*pow(Zz,2) + 8.*pow(Yy,6)*pow(Zz,2) + pow(Xx,2)*pow(Yy,2)*(120. - 44.*pow(Zz,2)) + pow(Xx,2)*pow(Yy,3)*(88. - 24.*pow(Zz,2)) + Xx*pow(Yy,4)*(43. - 12.*pow(Zz,2)) + pow(Xx,5)*(45. + 12.*Yy - 36.*pow(Yy,2) - 12.*pow(Zz,2)) + Xx*Yy*(24. - 8.*pow(Zz,2)) + pow(Xx,6)*(-38. - 12.*Yy + 32.*pow(Yy,2) + 8.*pow(Zz,2)) + Xx*pow(Yy,3)*(-32. + 12.*pow(Zz,2)) + pow(Xx,2)*Yy*(-52. + 19.*pow(Zz,2)) + Xx*pow(Yy,2)*(-52. + 19.*pow(Zz,2)) + pow(Xx,2)*pow(Yy,4)*(-114. + 24.*pow(Zz,2)) + pow(Xx,3)*(-56. + 24.*pow(Yy,3) - 36.*pow(Yy,4) + 21.*pow(Zz,2) + pow(Yy,2)*(88. - 24.*pow(Zz,2)) + Yy*(-32. + 12.*pow(Zz,2))) + pow(Xx,4)*(60. - 36.*pow(Yy,3) + 48.*pow(Yy,4) - 22.*pow(Zz,2) + Yy*(43. - 12.*pow(Zz,2)) + pow(Yy,2)*(-114. + 24.*pow(Zz,2))));
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    const double alpha = 1.0;
#warning using alpha = 1 here
    const double Xx = x[0], Yy = x[1];
    phi[ 0 ] = exp(-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1.*alpha + 1.*Xx - 2.*Xx*Xx + (1. - 2.*Yy)*Yy);  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    const double alpha = 1.0;
#warning using alpha = 1 here

    const double Xx = x[0], Yy = x[1];
    JacobianRangeType grad(0);
    grad[0][0] =  pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. + alpha*(1. - 2.*Xx) - 4.*pow(Xx,2) + 4.*pow(Xx,3) + 1.*Yy - 2.*pow(Yy,2) + Xx*(-3. - 2.*Yy + 4.*pow(Yy,2)));
    grad[0][1] =  pow(exp(1),-1.*(-1. + Xx)*Xx - 1.*(-1. + Yy)*Yy)*(1. + alpha*(1. - 2.*Yy) + Xx*(1. - 2.*Yy) - 3.*Yy - 4.*pow(Yy,2) + 4.*pow(Yy,3) + pow(Xx,2)*(-2. + 4.*Yy));

    double dot = 0.0;
    for( int i = 0; i < dimDomain; ++i )
      dot += grad[ 0 ][ i ] * x[ i ];

    for( int i = 0; i < 3; ++i )
      ret[ 0 ][ i ] = grad[ 0 ][ i ] - dot * x[ i ];
  }
  //! mass coefficient has to be 1 for this problem
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(1);
  }
  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x,
		 RangeType& value) const
  {
    value = RangeType( 0 );
  }
};

#endif // #ifndef POISSON_HH
