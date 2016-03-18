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

#include "temporalprobleminterface.hh"

double Power( const double y, const int a )
{
  assert( a > 0 );

  if( a == 1 )
    return y;

  return y * Power( y, a-1 );
}

template <class FunctionSpace>
class SurfaceHeatProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  SurfaceHeatProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    // define evolution of surface
    const double at = 1.0 + 0.25 * sin( time() );
    const double apt = 0.25 * cos( time() );

    // calculated surface parameters
    const double divGammaV = 0.5 * at * apt * ( x[1]*x[1] + x[2]*x[2] ) / ( x[0]*x[0] + at*at * ( x[1]*x[1] + x[2]*x[2] ) );
    const double N1 = 1/at * x[0] / sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] );
    const double N2 = x[1] / sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] );
    const double H = ( 2.0 * x[0] * x[0] + at * ( 1 + at ) * ( x[1]*x[1] + x[2]*x[2] ) )
      / ( sqrt( x[0]*x[0] / (at*at) + x[1]*x[1] + x[2]*x[2] ) * ( x[0]*x[0] + at*at * ( x[1]*x[1] + x[2]*x[2] ) ) );


    // calculate solution and derivatives
    const double ux = sin( time() ) * x[0] * x[1];
    const double mdux = ( cos( time() ) + 0.5 * sin( time() ) * apt / at ) * x[0] * x[1];
    const double mlapux = sin( time() ) * (  2.0 * N1 * N2 + H * ( x[1] * N1 + x[0] * N2 ) );

    // construct solution
    phi = mdux + divGammaV * ux + mlapux;
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0.0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to identity by default
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = 1;
  }

  //! advection coefficient (default = 0)
  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    // set to zero by default
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = RangeType(0);
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( time() ) * x[0] * x[1];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * x[1];
    grad[ 0 ][ 1 ] = sin( time() ) * x[0];
    grad[ 0 ][ 2 ] = 0.0;

    const double at = 1.0 + 0.25 * sin( time() );
    DomainType nu;
    nu[ 0 ] = 2.0 * x[ 0 ] / at;
    nu[ 1 ] = 2.0 * x[ 1 ];
    nu[ 2 ] = 2.0 * x[ 2 ];
    nu /= nu.two_norm();

    double dot = 0;
    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	dot += nu[ i ] * grad[ 0 ][ i ];
      }

    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	ret[ 0 ][ i ] = grad[ 0 ][ i ] - dot * nu[ i ];
      }
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return false ;
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  virtual bool isNeumannPoint( const DomainType& x ) const
  {
    return true ;
  }
};

template <class FunctionSpace>
class SurfaceMCProblem : public TemporalProblemInterface < FunctionSpace >
{
  typedef TemporalProblemInterface < FunctionSpace >  BaseType;
public:
  typedef typename BaseType :: RangeType            RangeType;
  typedef typename BaseType :: DomainType           DomainType;
  typedef typename BaseType :: JacobianRangeType    JacobianRangeType;
  typedef typename BaseType :: DiffusionTensorType  DiffusionTensorType;
  typedef typename BaseType :: AdvectionVectorType  AdvectionVectorType;

  enum { dimRange  = BaseType :: dimRange };
  enum { dimDomain = BaseType :: dimDomain };

  // get time function from base class
  using BaseType :: time ;
  using BaseType :: deltaT ;

  SurfaceMCProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {}

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    phi = RangeType(0);
    phi += Y1_;
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0.0;
  }

  //! diffusion coefficient (default = Id)
  virtual void D(const DomainType& x, DiffusionTensorType& D ) const
  {
    // set to identity by default
    D = 0;
    for( int i=0; i<D.rows; ++i )
      D[ i ][ i ] = 1;
  }

  //! advection coefficient (default = 0)
  virtual void b(const DomainType& x, AdvectionVectorType& b ) const
  {
    b = 0;
  }

  //! mass coefficient (default = 0)
  virtual void m(const DomainType& x, RangeType &m) const
  {
    m = 0;
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(1);
  }

  //! the exact solution
  virtual void u(const DomainType& x,
		 RangeType& phi) const
  {
    phi = sin( time() ) * x[0] * x[1];
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * x[1];
    grad[ 0 ][ 1 ] = sin( time() ) * x[0];
    grad[ 0 ][ 2 ] = 0.0;

    const double at = 1.0 + 0.25 * sin( time() );
    DomainType nu;
    nu[ 0 ] = 2.0 * x[ 0 ] / at;
    nu[ 1 ] = 2.0 * x[ 1 ];
    nu[ 2 ] = 2.0 * x[ 2 ];
    nu /= nu.two_norm();

    double dot = 0;
    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	dot += nu[ i ] * grad[ 0 ][ i ];
      }

    for( unsigned int i = 0; i < nu.size(); ++i )
      {
	ret[ 0 ][ i ] = grad[ 0 ][ i ] - dot * nu[ i ];
      }
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const
  {
    return false ;
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  virtual bool isNeumannPoint( const DomainType& x ) const
  {
    return true ;
  }

  virtual void setY1Y2( const double Y1, const double Y2 )
  {
    Y1_ = Y1;
    Y2_ = Y2;
  }

private:
  double Y1_;
  double Y2_;
};

#endif // #ifndef POISSON_HH
