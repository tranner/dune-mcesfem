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

  virtual void setYs( const std::vector<double> )
  {}
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

  double Power( const double y, const int a ) const
  {
    assert( a > 0 );

    if( a == 1 )
      return y;

    return y * Power( y, a-1 );
  }

  double Sin( const double x ) const
  {
    return sin( x );
  }

  double Cos( const double x ) const
  {
    return cos( x );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double xx = x[0];
    const double yy = x[1];
    const double zz = x[2];
    const double tt = time();

    const double Y1 = Ys_.at(0);
    const double Y2 = Ys_.at(1);

    phi = (-128*Power(xx,7)*Y1*yy*Power(Sin(tt),2)*(48 + Sin(tt)) - 1024*Power(xx,8)*Y1*(-20*Y2 + (6*Y1 - Y2)*Sin(tt))*Sin(2*tt) -
	   (Power(xx,3)*yy*Sin(tt)*Power(4 + Sin(tt),2)*(-2656 + 6168*Y1 + 160*Cos(tt) + 96*Cos(2*tt) - 24*Y1*Cos(2*tt) + 32*Y1*Power(zz,2)*(-227 + 3*Cos(2*tt) - 104*Sin(tt)) -
							 8*(Y1 + Y2)*Power(zz,4)*(-5 + 5*Cos(2*tt) - 96*Sin(tt)) + (-1155 + 2048*Y1)*Sin(tt) + 32*Sin(2*tt) + Sin(3*tt)))/2. +
	   Power(xx,5)*yy*Sin(tt)*(-1091 + 31040*Y1 + 16*Cos(tt) + 1092*Cos(2*tt) - 2368*Y1*Cos(2*tt) - 16*Cos(3*tt) - Cos(4*tt) + 8*(-542 + 3081*Y1)*Sin(tt) + 96*Sin(2*tt) + 80*Sin(3*tt) - 24*Y1*Sin(3*tt) +
				   48*Y1*Power(zz,2)*(-40 + 40*Cos(2*tt) - 259*Sin(tt) + Sin(3*tt))) + (xx*yy*Power(4 + Sin(tt),3)*
													(16 + 80*Y1 + 66*Cos(tt) - 16*Cos(2*tt) - 80*Y1*Cos(2*tt) - 2*Cos(3*tt) + (253 + 643*Y1)*Sin(tt) + 24*Sin(2*tt) + Sin(3*tt) - Y1*Sin(3*tt) +
													 2*Y1*Power(zz,2)*(-112 + 112*Cos(2*tt) - 777*Sin(tt) + 3*Sin(3*tt)) - (Y1 + Y2)*Power(zz,4)*(-144 + 144*Cos(2*tt) - 911*Sin(tt) + 5*Sin(3*tt))))/2. +
	   (Power(4 + Sin(tt),4)*(-4*(1 + Y1)*(Y1 - Y2)*Cos(tt) + 32*Y2*Cos(2*tt) + 4*Y1*Cos(3*tt) + 4*Power(Y1,2)*Cos(3*tt) - 4*Y2*Cos(3*tt) - 4*Y1*Y2*Cos(3*tt) - 3*Y2*Sin(tt) - 32*Y1*Sin(2*tt) - 32*Power(Y1,2)*Sin(2*tt) +
				  64*Y2*Sin(2*tt) + 64*Y1*Y2*Sin(2*tt) - 8*Power(zz,4)*(4*(Power(Y1,2) - 15*Y1*Y2 - 6*Power(Y2,2)) + (Power(Y1,2) - 12*Y1*Y2 - 5*Power(Y2,2))*Sin(tt))*Sin(2*tt) -
				  8*Y2*(Y1 + Y2)*Power(zz,6)*(3*Cos(tt) - 3*Cos(3*tt) + 28*Sin(2*tt)) + 5*Y2*Sin(3*tt) +
				  Power(zz,2)*(8*(Power(Y1,2) - Y2 - 4*Y1*Y2)*Cos(tt) - 32*Y2*Cos(2*tt) - 8*Power(Y1,2)*Cos(3*tt) + 8*Y2*Cos(3*tt) + 32*Y1*Y2*Cos(3*tt) + 3*Y2*Sin(tt) + 64*Power(Y1,2)*Sin(2*tt) - 96*Y2*Sin(2*tt) -
					       352*Y1*Y2*Sin(2*tt) - 5*Y2*Sin(3*tt))))/2. + Power(xx,4)*(-1 + Cos(2*tt) - 8*Sin(tt))*
	   (31*Y1 - 32*Y2 - 5392*Y1*Cos(tt) + 15696*Power(Y1,2)*Cos(tt) + 5280*Y2*Cos(tt) - 16480*Y1*Y2*Cos(tt) + 296*Y1*Cos(2*tt) - 288*Y2*Cos(2*tt) + 272*Y1*Cos(3*tt) - 336*Power(Y1,2)*Cos(3*tt) - 160*Y2*Cos(3*tt) +
	    96*Y1*Y2*Cos(3*tt) - 7*Y1*Cos(4*tt) - 36*(Y1 - Y2)*Sin(tt) + 4*Power(zz,2)*(-8*(461*Power(Y1,2) + 3*Y2 - 684*Y1*Y2)*Cos(tt) + 8*(13*Power(Y1,2) + 3*Y2 - 12*Y1*Y2)*Cos(3*tt) +
											2*((-1280*Power(Y1,2) - 129*Y2 + 1760*Y1*Y2)*Cos(tt) + Y2*(1 + 5*Cos(2*tt) + Cos(3*tt)))*Sin(tt)) - 1948*Y1*Sin(2*tt) + 4608*Power(Y1,2)*Sin(2*tt) + 1544*Y2*Sin(2*tt) - 3456*Y1*Y2*Sin(2*tt) +
	    32*Power(zz,4)*(24*(Power(Y1,2) - 2*Y1*Y2 - Power(Y2,2)) + (6*Power(Y1,2) - 7*Y1*Y2 - 5*Power(Y2,2))*Sin(tt))*Sin(2*tt) + 92*Y1*Sin(3*tt) - 60*Y2*Sin(3*tt) + 14*Y1*Sin(4*tt) - 4*Y2*Sin(4*tt)) -
	   4*Power(xx,6)*Sin(tt)*(Y1 + 256*Y1*Cos(tt) - 14944*Power(Y1,2)*Cos(tt) - 160*Y2*Cos(tt) + 24704*Y1*Y2*Cos(tt) - 10240*Y1*Y2*Power(yy,2)*Cos(tt) - 8*Y1*Cos(2*tt) - 256*Y1*Cos(3*tt) + 608*Power(Y1,2)*Cos(3*tt) +
				  160*Y2*Cos(3*tt) - 128*Y1*Y2*Cos(3*tt) + 7*Y1*Cos(4*tt) + 24*(Y1 - Y2)*Sin(tt) - 256*Y1*Power(zz,2)*Cos(tt)*(-3*Y1 + 42*Y2 + (3*Y1 - 2*Y2)*Cos(2*tt) + (-24*Y1 + 36*Y2)*Sin(tt)) + 1056*Y1*Sin(2*tt) -
				  6656*Power(Y1,2)*Sin(2*tt) - 1032*Y2*Sin(2*tt) + 5376*Y1*Y2*Sin(2*tt) - 40*Y1*Sin(3*tt) + 40*Y2*Sin(3*tt) - 16*Y1*Sin(4*tt) + 4*Y2*Sin(4*tt)) -
	   (Power(xx,2)*Power(4 + Sin(tt),2)*(48*Y1 - 72*Y2 - 24*(-33*Y1 + 192*Power(Y1,2) + 43*Y2 - 160*Y1*Y2)*Cos(tt) - 2240*Y1*Cos(2*tt) + 2240*Y2*Cos(2*tt) - 804*Y1*Cos(3*tt) + 4608*Power(Y1,2)*Cos(3*tt) +
					      1036*Y2*Cos(3*tt) - 3840*Y1*Y2*Cos(3*tt) + 144*Y1*Cos(4*tt) - 120*Y2*Cos(4*tt) + 12*Y1*Cos(5*tt) - 4*Y2*Cos(5*tt) + (586*Y1 - 896*Y2)*Sin(tt) + 384*Y1*Sin(2*tt) - 19008*Power(Y1,2)*Sin(2*tt) -
					      2304*Y2*Sin(2*tt) + 22784*Y1*Y2*Sin(2*tt) - 768*Y2*(Y1 + Y2)*Power(zz,6)*Power(Sin(tt),2)*(8*Cos(tt) + Sin(2*tt)) -
					      64*Power(zz,4)*(28*(Power(Y1,2) - 2*Y1*Y2 - Power(Y2,2)) + (7*Power(Y1,2) - 19*Y1*Y2 - 10*Power(Y2,2))*Sin(tt))*(Cos(tt) - Cos(3*tt) + 8*Sin(2*tt)) - 975*Y1*Sin(3*tt) + 1152*Y2*Sin(3*tt) - 192*Y1*Sin(4*tt) +
					      288*Power(Y1,2)*Sin(4*tt) + 128*Y2*Sin(4*tt) - 128*Y1*Y2*Sin(4*tt) - 16*Power(zz,2)*Sin(tt)*
					      (8*Y2 - 4*(1040*Power(Y1,2) + 327*Y2 - 2072*Y1*Y2)*Cos(tt) + 72*Y2*Cos(2*tt) + 64*Power(Y1,2)*Cos(3*tt) + 28*Y2*Cos(3*tt) - 96*Y1*Y2*Cos(3*tt) - 6*Y2*Sin(tt) - 1024*Power(Y1,2)*Sin(2*tt) -
					       322*Y2*Sin(2*tt) + 1792*Y1*Y2*Sin(2*tt) + 10*Y2*Sin(3*tt) + Y2*Sin(4*tt)) + 7*Y1*Sin(5*tt)))/8.)/(2.*(4 + Sin(tt))*Power(-4*Power(xx,2)*Sin(tt) + Power(4 + Sin(tt),2),2));
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

    D *= ( 1 + x[0]*x[0] + Ys_.at(0) * x[1]*x[1]*x[1]*x[1] + Ys_.at(1) * x[2]*x[2]*x[2]*x[2] );
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

  virtual void setYs( const std::vector<double> Ys )
  {
    Ys_ = Ys;
  }

private:
  std::vector<double> Ys_;
};

template <class FunctionSpace>
class SurfaceMCStationaryProblem : public TemporalProblemInterface < FunctionSpace >
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

  SurfaceMCStationaryProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {
    assert( dimDomain == 2 );
  }

  double Cos( const double a ) const { return cos( a ); }
  double Sin( const double a ) const { return sin( a ); }
  double Power( const double a, int n ) const
  {
    if( n == 1 )
      return a;
    else
      return a * Power( a, n-1 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double t = time();
    const double xx = x[0];
    const double yy = x[1];

    const double Y1 = Ys_.at(0);
    const double Y2 = Ys_.at(1);

    phi = (4*Cos(t)*Cos(3*xx) + 4*Y1*Cos(t)*Cos(5*xx) + 4*Cos(t)*Cos(3*yy) + 4*Y2*Cos(t)*Cos(5*yy) + 36*Cos(3*yy)*Sin(t) + 100*Y2*Cos(5*yy)*Sin(t) + 9*Y1*Cos(3*yy)*Sin(t)*Sin(2*xx) + 25*Y1*Y2*Cos(5*yy)*Sin(t)*Sin(2*xx) + 9*Y2*Cos(3*yy)*Sin(t)*Sin(2*yy) + 25*Power(Y2,2)*Cos(5*yy)*Sin(t)*Sin(2*yy) - xx*Sin(t)*(3*Sin(3*xx) + 5*Y1*Sin(5*xx))*(4 + Y1*Sin(2*xx) + Y2*Sin(2*yy)) + 3*Y2*Sin(t)*(Sin(yy) + Sin(5*yy)) - yy*Sin(t)*(4 + Y1*Sin(2*xx) + Y2*Sin(2*yy))*(3*Sin(3*yy) + 5*Y2*Sin(5*yy)) - 2*xx*yy*Sin(t)*(Y2*Cos(2*yy)*(3*Sin(3*xx) + 5*Y1*Sin(5*xx)) + Y1*Cos(2*xx)*(3*Sin(3*yy) + 5*Y2*Sin(5*yy))) + 5*Power(Y2,2)*Sin(t)*(Sin(3*yy) + Sin(7*yy)) + Power(yy,2)*Sin(t)*(-36*Cos(3*yy) - 100*Y2*Cos(5*yy) + 3*Y1*Sin(xx) - 9*Y1*Cos(3*yy)*Sin(2*xx) - 25*Y1*Y2*Cos(5*yy)*Sin(2*xx) + 5*Power(Y1,2)*Sin(3*xx) + 3*Y1*Sin(5*xx) + 5*Power(Y1,2)*Sin(7*xx) - 3*Y2*Sin(yy) - 9*Y2*Cos(3*yy)*Sin(2*yy) - 25*Power(Y2,2)*Cos(5*yy)*Sin(2*yy) + 9*Cos(3*xx)*(4 + Y1*Sin(2*xx) + Y2*Sin(2*yy)) + 25*Y1*Cos(5*xx)*(4 + Y1*Sin(2*xx) + Y2*Sin(2*yy)) - 5*Power(Y2,2)*Sin(3*yy) - 3*Y2*Sin(5*yy) - 5*Power(Y2,2)*Sin(7*yy)))/4;
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

    D *= ( 1 + 0.25 * Ys_.at(0) * sin( 2 * x[0] ) + 0.25 * Ys_.at(1) * sin( 2 * x[1] ) );
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
    phi = sin( time() ) * ( cos( 3*x[0] ) + cos( 3 * x[1] ) );
    // + sin( time() ) * ( Ys_.at(0) * cos( 5 * x[0] ) + Ys_.at(1) * cos( 5 * x[1] ) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * ( -3*sin( 3*x[0] ) );
    grad[ 0 ][ 1 ] = sin( time() ) * ( -3*sin( 3*x[1] ) );

    DomainType nu = x;

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

  virtual void setYs( const std::vector<double> Ys )
  {
    Ys_ = Ys;
  }

private:
  std::vector<double> Ys_;
};

template <class FunctionSpace>
class CurveMCProblem : public TemporalProblemInterface < FunctionSpace >
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

  CurveMCProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {
    assert( dimDomain == 2 );
  }

  double Cos( const double a ) const { return cos( a ); }
  double Sin( const double a ) const { return sin( a ); }
  double Power( const double a, int n ) const
  {
    if( n == 1 )
      return a;
    else
      return a * Power( a, n-1 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& z,
		 RangeType& phi) const
  {
    const double t = time();
    const double x = z[0];
    const double y = z[1];
    const double Y1 = Ys_.at(0);
    const double Y2 = Ys_.at(1);

    phi = (8*Cos(t)*(Cos(3*x) + Y1*Cos(5*x) + Cos(3*y) + Y2*Cos(5*y)) + 
     (4*(Cos(3*x) + Y1*Cos(5*x) + Cos(3*y) + Y2*Cos(5*y))*Sin(t)*
        (-4*Power(x,2)*Sin(t) + 
          Cos(t)*(4*Power(y,2) + (-Power(x,2) + Power(y,2))*Sin(t))))/
      (16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
        Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
        Power(y,2)*Power(Sin(t),2)) - 
     (4*x*Cos(t)*Sin(t)*(3*Sin(3*x) + 5*Y1*Sin(5*x)))/(4 + Sin(t)) + 
     (x*(4 + Cos(t))*Sin(t)*(9*x*(4 + Cos(t))*Cos(3*y)*
           (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
             (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t)) + 
          25*x*Y2*(4 + Cos(t))*Cos(5*y)*
           (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
             (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t)) - 
          (4 + Sin(t))*(3*(33*Power(x,2) - 33*Power(y,2) + 
                16*Power(x,2)*Cos(t) + (Power(x,2) + Power(y,2))*Cos(2*t) - 
                16*Power(y,2)*Sin(t))*Sin(3*x) + 
             5*Y1*(33*Power(x,2) - 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
                (Power(x,2) + Power(y,2))*Cos(2*t) - 16*Power(y,2)*Sin(t))*
              Sin(5*x) + 4*x*y*(4 + Cos(t))*
              (3 + 5*Y2 + 2*(3 + 5*Y2)*Cos(2*y) + 10*Y2*Cos(4*y))*
              (4 + Sin(t))*Sin(y)))*(4 + Y1*Sin(2*x) + Y2*Sin(2*y)))/
      Power(16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
        Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
        Power(y,2)*Power(Sin(t),2),2) + 
     (4*y*Power(Sin(t),2)*(3*Sin(3*y) + 5*Y2*Sin(5*y)))/(4 + Cos(t)) + 
     (2*y*Sin(t)*(4 + Sin(t))*(4 + Y1*Sin(2*x) + Y2*Sin(2*y))*
        ((9*y*Cos(3*x)*(4 + Sin(t))*
             (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
               (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t)))/
           2. + (25*y*Y1*Cos(5*x)*(4 + Sin(t))*
             (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
               (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t)))/
           2. + (4 + Cos(t))*(-3*x*y*(32 + 8*Cos(t) + 8*Sin(t) + Sin(2*t))*
              Sin(3*x) - 5*x*y*Y1*(32 + 8*Cos(t) + 8*Sin(t) + Sin(2*t))*
              Sin(5*x) + (16*(Power(x,2) - Power(y,2)) + 
                8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) - 
                8*Power(y,2)*Sin(t) - Power(y,2)*Power(Sin(t),2))*
              (3*Sin(3*y) + 5*Y2*Sin(5*y)))))/
      Power(16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
        Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
        Power(y,2)*Power(Sin(t),2),2) + 
     (y*Y1*Cos(2*x)*(-1 + Cos(2*t) - 8*Sin(t))*
        (-5*y*Y1*Cos(t - 5*x) - 3*y*Cos(t - 3*x) + 3*y*Cos(t + 3*x) + 
          5*y*Y1*Cos(t + 5*x) - 24*y*Sin(3*x) - 40*y*Y1*Sin(5*x) - 
          5*x*Y2*Sin(t - 5*y) - 3*x*Sin(t - 3*y) + 24*x*Sin(3*y) + 
          40*x*Y2*Sin(5*y) + 3*x*Sin(t + 3*y) + 5*x*Y2*Sin(t + 5*y)))/
      (16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
        Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
        Power(y,2)*Power(Sin(t),2)) + 
     (x*Y2*Cos(2*y)*(8*Sin(t) + Sin(2*t))*
        (-5*y*Y1*Cos(t - 5*x) - 3*y*Cos(t - 3*x) + 3*y*Cos(t + 3*x) + 
          5*y*Y1*Cos(t + 5*x) - 24*y*Sin(3*x) - 40*y*Y1*Sin(5*x) - 
          5*x*Y2*Sin(t - 5*y) - 3*x*Sin(t - 3*y) + 24*x*Sin(3*y) + 
          40*x*Y2*Sin(5*y) + 3*x*Sin(t + 3*y) + 5*x*Y2*Sin(t + 5*y)))/
      (16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
        Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
        Power(y,2)*Power(Sin(t),2)) + 
     (Power(x,2)*y*(-2*Sin(t)*(4 + Sin(t))*(4 + Y1*Sin(2*x) + Y2*Sin(2*y))*
           ((9*y*Cos(3*x)*(4 + Sin(t))*
                (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
                  (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t))
                )/2. + (25*y*Y1*Cos(5*x)*(4 + Sin(t))*
                (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
                  (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t))
                )/2. + (4 + Cos(t))*
              (-3*x*y*(32 + 8*Cos(t) + 8*Sin(t) + Sin(2*t))*Sin(3*x) - 
                5*x*y*Y1*(32 + 8*Cos(t) + 8*Sin(t) + Sin(2*t))*Sin(5*x) + 
                (16*(Power(x,2) - Power(y,2)) + 8*Power(x,2)*Cos(t) + 
                   Power(x,2)*Power(Cos(t),2) - 8*Power(y,2)*Sin(t) - 
                   Power(y,2)*Power(Sin(t),2))*(3*Sin(3*y) + 5*Y2*Sin(5*y)))
             ) + Y1*Cos(2*x)*(-1 + Cos(2*t) - 8*Sin(t))*
           (16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
             Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
             Power(y,2)*Power(Sin(t),2))*
           (5*y*Y1*Cos(t - 5*x) + 3*y*Cos(t - 3*x) - 3*y*Cos(t + 3*x) - 
             5*y*Y1*Cos(t + 5*x) + 24*y*Sin(3*x) + 40*y*Y1*Sin(5*x) + 
             5*x*Y2*Sin(t - 5*y) + 3*x*Sin(t - 3*y) - 24*x*Sin(3*y) - 
             40*x*Y2*Sin(5*y) - 3*x*Sin(t + 3*y) - 5*x*Y2*Sin(t + 5*y))))/
      (Power(4 + Sin(t),2)*Power(16*(Power(x,2) + Power(y,2)) + 
          8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) + 
          8*Power(y,2)*Sin(t) + Power(y,2)*Power(Sin(t),2),2)*
        (Power(y,2)/Power(4 + Cos(t),2) + Power(x,2)/Power(4 + Sin(t),2)))\
      + (x*y*(2*x*(4 + Cos(t))*Sin(t)*(4 + Sin(t))*
           (4 + Y1*Sin(2*x) + Y2*Sin(2*y))*
           (9*y*Cos(3*y)*(16*(Power(x,2) + Power(y,2)) + 
                8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) + 
                8*Power(y,2)*Sin(t) + Power(y,2)*Power(Sin(t),2)) + 
             25*y*Y2*Cos(5*y)*(16*(Power(x,2) + Power(y,2)) + 
                8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) + 
                8*Power(y,2)*Sin(t) + Power(y,2)*Power(Sin(t),2)) - 
             96*x*y*Sin(3*x) - 24*x*y*Cos(t)*Sin(3*x) - 
             24*x*y*Sin(t)*Sin(3*x) - 3*x*y*Sin(2*t)*Sin(3*x) - 
             160*x*y*Y1*Sin(5*x) - 40*x*y*Y1*Cos(t)*Sin(5*x) - 
             40*x*y*Y1*Sin(t)*Sin(5*x) - 5*x*y*Y1*Sin(2*t)*Sin(5*x) + 
             48*Power(x,2)*Sin(3*y) - 48*Power(y,2)*Sin(3*y) + 
             24*Power(x,2)*Cos(t)*Sin(3*y) + 
             3*Power(x,2)*Power(Cos(t),2)*Sin(3*y) - 
             24*Power(y,2)*Sin(t)*Sin(3*y) - 
             3*Power(y,2)*Power(Sin(t),2)*Sin(3*y) + 
             80*Power(x,2)*Y2*Sin(5*y) - 80*Power(y,2)*Y2*Sin(5*y) + 
             40*Power(x,2)*Y2*Cos(t)*Sin(5*y) + 
             5*Power(x,2)*Y2*Power(Cos(t),2)*Sin(5*y) - 
             40*Power(y,2)*Y2*Sin(t)*Sin(5*y) - 
             5*Power(y,2)*Y2*Power(Sin(t),2)*Sin(5*y)) + 
          y*Y2*Cos(2*y)*(-1 + Cos(2*t) - 8*Sin(t))*
           (16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
             Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
             Power(y,2)*Power(Sin(t),2))*
           (5*y*Y1*Cos(t - 5*x) + 3*y*Cos(t - 3*x) - 3*y*Cos(t + 3*x) - 
             5*y*Y1*Cos(t + 5*x) + 24*y*Sin(3*x) + 40*y*Y1*Sin(5*x) + 
             5*x*Y2*Sin(t - 5*y) + 3*x*Sin(t - 3*y) - 24*x*Sin(3*y) - 
             40*x*Y2*Sin(5*y) - 3*x*Sin(t + 3*y) - 5*x*Y2*Sin(t + 5*y))))/
      ((4 + Cos(t))*(4 + Sin(t))*
        Power(16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
          Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
          Power(y,2)*Power(Sin(t),2),2)*
        (Power(y,2)/Power(4 + Cos(t),2) + Power(x,2)/Power(4 + Sin(t),2)))\
      + (x*y*(y*(4 + Cos(t))*Sin(t)*(4 + Sin(t))*
           (25*x*Y1*Cos(5*x)*(33*Power(x,2) + 33*Power(y,2) + 
                16*Power(x,2)*Cos(t) + (Power(x,2) - Power(y,2))*Cos(2*t) + 
                16*Power(y,2)*Sin(t)) + 
             18*x*Cos(3*x)*(16*(Power(x,2) + Power(y,2)) + 
                8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) + 
                8*Power(y,2)*Sin(t) + Power(y,2)*Power(Sin(t),2)) - 
             2*(3*(16*(Power(x,2) - Power(y,2)) + 8*Power(x,2)*Cos(t) + 
                   Power(x,2)*Power(Cos(t),2) - 8*Power(y,2)*Sin(t) - 
                   Power(y,2)*Power(Sin(t),2))*Sin(3*x) + 
                (5*Y1*(33*Power(x,2) - 33*Power(y,2) + 
                     16*Power(x,2)*Cos(t) + 
                     (Power(x,2) + Power(y,2))*Cos(2*t) - 
                     16*Power(y,2)*Sin(t))*Sin(5*x))/2. + 
                x*y*(3 + 5*Y2 + 2*(3 + 5*Y2)*Cos(2*y) + 10*Y2*Cos(4*y))*
                 (32 + 8*Cos(t) + 8*Sin(t) + Sin(2*t))*Sin(y)))*
           (4 + Y1*Sin(2*x) + Y2*Sin(2*y)) + 
          x*Y1*Cos(2*x)*(16*(Power(x,2) + Power(y,2)) + 
             8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) + 
             8*Power(y,2)*Sin(t) + Power(y,2)*Power(Sin(t),2))*
           (8*Sin(t) + Sin(2*t))*
           (5*y*Y1*Cos(t - 5*x) + 3*y*Cos(t - 3*x) - 3*y*Cos(t + 3*x) - 
             5*y*Y1*Cos(t + 5*x) + 24*y*Sin(3*x) + 40*y*Y1*Sin(5*x) + 
             5*x*Y2*Sin(t - 5*y) + 3*x*Sin(t - 3*y) - 24*x*Sin(3*y) - 
             40*x*Y2*Sin(5*y) - 3*x*Sin(t + 3*y) - 5*x*Y2*Sin(t + 5*y))))/
      ((4 + Cos(t))*(4 + Sin(t))*
        Power(16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
          Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
          Power(y,2)*Power(Sin(t),2),2)*
        (Power(y,2)/Power(4 + Cos(t),2) + Power(x,2)/Power(4 + Sin(t),2)))\
      + (x*Power(y,2)*(-((4 + Cos(t))*Sin(t)*
             (9*x*(4 + Cos(t))*Cos(3*y)*
                (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
                  (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t))
                 + 25*x*Y2*(4 + Cos(t))*Cos(5*y)*
                (33*Power(x,2) + 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
                  (Power(x,2) - Power(y,2))*Cos(2*t) + 16*Power(y,2)*Sin(t))
                 - (4 + Sin(t))*
                (3*(33*Power(x,2) - 33*Power(y,2) + 16*Power(x,2)*Cos(t) + 
                     (Power(x,2) + Power(y,2))*Cos(2*t) - 
                     16*Power(y,2)*Sin(t))*Sin(3*x) + 
                  5*Y1*(33*Power(x,2) - 33*Power(y,2) + 
                     16*Power(x,2)*Cos(t) + 
                     (Power(x,2) + Power(y,2))*Cos(2*t) - 
                     16*Power(y,2)*Sin(t))*Sin(5*x) + 
                  4*x*y*(4 + Cos(t))*
                   (3 + 5*Y2 + 2*(3 + 5*Y2)*Cos(2*y) + 10*Y2*Cos(4*y))*
                   (4 + Sin(t))*Sin(y)))*(4 + Y1*Sin(2*x) + Y2*Sin(2*y))) + 
          Y2*Cos(2*y)*(16*(Power(x,2) + Power(y,2)) + 8*Power(x,2)*Cos(t) + 
             Power(x,2)*Power(Cos(t),2) + 8*Power(y,2)*Sin(t) + 
             Power(y,2)*Power(Sin(t),2))*(8*Sin(t) + Sin(2*t))*
           (5*y*Y1*Cos(t - 5*x) + 3*y*Cos(t - 3*x) - 3*y*Cos(t + 3*x) - 
             5*y*Y1*Cos(t + 5*x) + 24*y*Sin(3*x) + 40*y*Y1*Sin(5*x) + 
             5*x*Y2*Sin(t - 5*y) + 3*x*Sin(t - 3*y) - 24*x*Sin(3*y) - 
             40*x*Y2*Sin(5*y) - 3*x*Sin(t + 3*y) - 5*x*Y2*Sin(t + 5*y))))/
      (Power(4 + Cos(t),2)*Power(16*(Power(x,2) + Power(y,2)) + 
          8*Power(x,2)*Cos(t) + Power(x,2)*Power(Cos(t),2) + 
          8*Power(y,2)*Sin(t) + Power(y,2)*Power(Sin(t),2),2)*
        (Power(y,2)/Power(4 + Cos(t),2) + Power(x,2)/Power(4 + Sin(t),2))))/
      8.; 
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

    D *= ( 1 + 0.25 * Ys_.at(0) * sin( 2 * x[0] ) + 0.25 * Ys_.at(1) * sin( 2 * x[1] ) );
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
    phi = sin( time() ) * ( cos( 3*x[0] ) + cos( 3 * x[1] ) );
    // phi += sin( time() ) * ( Ys_.at(0) * cos( 5 * x[0] ) + Ys_.at(1) * cos( 5 * x[1] ) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * ( -3*sin( 3*x[0] ) );
    grad[ 0 ][ 1 ] = sin( time() ) * ( -3*sin( 3*x[1] ) );
    // grad[ 0 ][ 0 ] += Ys_.at(0) * sin( time() ) * ( -5*sin( 5*x[0] ) );
    // grad[ 0 ][ 1 ] += Ys_.at(1) * sin( time() ) * ( -5*sin( 5*x[1] ) );

    DomainType nu = x;

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

  virtual void setYs( const std::vector<double> Ys )
  {
    Ys_ = Ys;
  }

private:
  std::vector<double> Ys_;
};

template <class FunctionSpace>
class CurveMCLessSmoothProblem : public TemporalProblemInterface < FunctionSpace >
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

  CurveMCLessSmoothProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {
    assert( dimDomain == 2 );
  }

  double Cos( const double a ) const { return cos( a ); }
  double Sin( const double a ) const { return sin( a ); }
  double Power( const double a, int n ) const
  {
    if( n == 1 )
      return a;
    else
      return a * Power( a, n-1 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double t = time();
    const double Y1 = Ys_.at(0);
    const double Y2 = Ys_.at(1);

    const auto eta = []( const double s ) { return s*std::abs(s); };
    const auto deta = []( const double s ){ return 2.0*std::abs(s); };

    RangeType f_md = Cos(t)*(Cos(3*x[0]) + Y1*Cos(5*x[0]) + Cos(3*x[1]) + Y2*Cos(5*x[1])) - (Cos(t)*x[0]*Sin(t)*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0])))/(2.*(4 + Sin(t))) + (x[1]*Power(Sin(t),2)*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])))/(2.*(4 + Cos(t)));
    RangeType f_divVu = -((Cos(3*x[0]) + Y1*Cos(5*x[0]) + Cos(3*x[1]) + Y2*Cos(5*x[1]))*Sin(t)*(8*Sin(t) - 8*Power(x[1],2)*(Cos(t) + Sin(t)) + Sin(2*t)))/(4.*(Power(4 + Cos(t),2) - 4*Power(x[1],2)*(Cos(t) - Sin(t))));

    RangeType f_lapu = (Sin(t)*(9216*Cos(3*x[1]) + 25600*Y2*Cos(5*x[1]) + 2304*Y2*Cos(3*x[1])*eta(x[1]) + 6400*Power(Y2,2)*Cos(5*x[1])*eta(x[1]) - 256*(4 + Y2*eta(x[1]))*x[0]*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0])) - 192*Cos(t)*(4 + Y2*eta(x[1]))*x[0]*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0])) - 48*Power(Cos(t),2)*(4 + Y2*eta(x[1]))*x[0]*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0])) - 4*Power(Cos(t),3)*(4 + Y2*eta(x[1]))*x[0]*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0])) - 256*(4 + Y2*eta(x[1]))*x[1]*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])) - 128*Cos(t)*(4 + Y2*eta(x[1]))*x[1]*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])) - 16*Power(Cos(t),2)*(4 + Y2*eta(x[1]))*x[1]*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])) - 64*(4 + Y2*eta(x[1]))*x[1]*Sin(t)*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])) - 32*Cos(t)*(4 + Y2*eta(x[1]))*x[1]*Sin(t)*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])) - 4*Power(Cos(t),2)*(4 + Y2*eta(x[1]))*x[1]*Sin(t)*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1])) - Y1*eta(x[0])*(-(Power(4 + Cos(t),3)*((4 + Cos(t))*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) - 4*x[0]*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0])))) + 8*Power(x[1],4)*(Cos(t) - Sin(t))*(72*Cos(3*x[0]) + 200*Y1*Cos(5*x[0]) - 25*Y2*Cos(t - 5*x[1]) - 9*Cos(t - 3*x[1]) - 72*Cos(3*x[1]) - 200*Y2*Cos(5*x[1]) - 9*Cos(t + 3*x[1]) - 25*Y2*Cos(t + 5*x[1]) + 25*Y1*Sin(t - 5*x[0]) + 9*Sin(t - 3*x[0]) + 9*Sin(t + 3*x[0]) + 25*Y1*Sin(t + 5*x[0])) + 4*Power(4 + Cos(t),2)*(3 + 5*Y2 + 2*(3 + 5*Y2)*Cos(2*x[1]) + 10*Y2*Cos(4*x[1]))*x[1]*(4 + Sin(t))*Sin(x[1]) - 2*Power(4 + Cos(t),2)*Power(x[1],2)*(72*Cos(3*x[0]) + 200*Y1*Cos(5*x[0]) - 50*Y2*Cos(t - 5*x[1]) - 18*Cos(t - 3*x[1]) - 72*Cos(3*x[1]) - 200*Y2*Cos(5*x[1]) - 18*Cos(t + 3*x[1]) - 50*Y2*Cos(t + 5*x[1]) + 25*Y1*Sin(t - 5*x[0]) + 9*Sin(t - 3*x[0]) + 9*Sin(t + 3*x[0]) + 25*Y1*Sin(t + 5*x[0]) + 25*Y2*Sin(t - 5*x[1]) + 9*Sin(t - 3*x[1]) + 9*Sin(t + 3*x[1]) + 25*Y2*Sin(t + 5*x[1]))) + 16*Power(x[1],4)*Power(Sin(t),2)*(4*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0])) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0])) + 768*Y2*Sin(3*x[1])*deta(x[1]) + 1280*Power(Y2,2)*Sin(5*x[1])*deta(x[1]) - 256*x[0]*x[1]*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) - 192*Cos(t)*x[0]*x[1]*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) - 48*Power(Cos(t),2)*x[0]*x[1]*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) - 4*Power(Cos(t),3)*x[0]*x[1]*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) + 64*Cos(t)*x[0]*Power(x[1],3)*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) + 16*Power(Cos(t),2)*x[0]*Power(x[1],3)*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) - 64*x[0]*Power(x[1],3)*Sin(t)*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) - 8*x[0]*Power(x[1],3)*Sin(2*t)*(Y1*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[0]) + Y2*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[1])) + 16*Power(Cos(t),2)*Power(x[1],2)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) - 180*Cos(3*x[1]) - 500*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) - 5*(-5*Y1*Cos(5*x[0]) + 9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) - 15*Y2*Sin(3*x[1])*deta(x[1]) - 25*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) + 256*Power(x[1],2)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) - 36*Cos(3*x[1]) - 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) - 9*Cos(3*x[1]) - 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) - 3*Y2*Sin(3*x[1])*deta(x[1]) - 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) - 64*Cos(t)*Power(x[1],4)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) - 36*Cos(3*x[1]) - 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) - 9*Cos(3*x[1]) - 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) - 3*Y2*Sin(3*x[1])*deta(x[1]) - 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) + 64*Power(x[1],4)*Sin(t)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) - 36*Cos(3*x[1]) - 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) - 9*Cos(3*x[1]) - 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) - 3*Y2*Sin(3*x[1])*deta(x[1]) - 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) + 64*Power(x[1],2)*Sin(t)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) + 36*Cos(3*x[1]) + 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) + 9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) + 3*Y2*Sin(3*x[1])*deta(x[1]) + 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) + 32*Cos(t)*Power(x[1],2)*Sin(t)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) + 36*Cos(3*x[1]) + 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) + 9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) + 3*Y2*Sin(3*x[1])*deta(x[1]) + 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) + 4*Power(Cos(t),2)*Power(x[1],2)*Sin(t)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) + 36*Cos(3*x[1]) + 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) + 9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) + 3*Y2*Sin(3*x[1])*deta(x[1]) + 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) - 8*Power(x[1],4)*Sin(2*t)*(36*Cos(3*x[0]) + 100*Y1*Cos(5*x[0]) + 36*Cos(3*x[1]) + 100*Y2*Cos(5*x[1]) + Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) + 9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) + 3*Y2*Sin(3*x[1])*deta(x[1]) + 5*Power(Y2,2)*Sin(5*x[1])*deta(x[1])) + 256*Cos(t)*(4*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) + Y2*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])) + 96*Power(Cos(t),2)*(4*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) + Y2*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])) + 16*Power(Cos(t),3)*(4*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) + Y2*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])) + Power(Cos(t),4)*(4*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) + Y2*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])) - 8*Power(Cos(t),3)*Power(x[1],2)*(4*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) + Y2*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])) + 16*Power(Cos(t),2)*Power(x[1],4)*(4*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1])) + Y2*(9*Cos(3*x[1]) + 25*Y2*Cos(5*x[1]))*eta(x[1]) + Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])) + 128*Cos(t)*Power(x[1],2)*(Y2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) - 18*Cos(3*x[1]) - 50*Y2*Cos(5*x[1]))*eta(x[1]) + Y1*(3*Sin(3*x[0]) + 5*Y1*Sin(5*x[0]))*deta(x[0]) + 2*(2*(9*Cos(3*x[0]) + 25*Y1*Cos(5*x[0]) - 18*Cos(3*x[1]) - 50*Y2*Cos(5*x[1])) - Y2*(3*Sin(3*x[1]) + 5*Y2*Sin(5*x[1]))*deta(x[1])))))/Power(2*Power(4 + Cos(t),2) - 8*Power(x[1],2)*(Cos(t) - Sin(t)),2);

    phi = f_md + f_divVu + f_lapu;
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

    D *= ( 1.0 + 0.25 * Ys_.at(0) * x[0] * std::abs(x[0])
	   + 0.25 * Ys_.at(1) * x[1] * std::abs(x[1]) );
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
    phi = sin( time() ) * ( cos( 3*x[0] ) + cos( 3 * x[1] ) );
    // phi += sin( time() ) * ( Ys_.at(0) * cos( 5 * x[0] ) + Ys_.at(1) * cos( 5 * x[1] ) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * ( -3*sin( 3*x[0] ) );
    grad[ 0 ][ 1 ] = sin( time() ) * ( -3*sin( 3*x[1] ) );
    // grad[ 0 ][ 0 ] += Ys_.at(0) * sin( time() ) * ( -5*sin( 5*x[0] ) );
    // grad[ 0 ][ 1 ] += Ys_.at(1) * sin( time() ) * ( -5*sin( 5*x[1] ) );

    DomainType nu = x;

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

  virtual void setYs( const std::vector<double> Ys )
  {
    Ys_ = Ys;
  }

private:
  std::vector<double> Ys_;
};

template <class FunctionSpace>
class CurveMCMoreNonlinearProblem : public TemporalProblemInterface < FunctionSpace >
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

  CurveMCMoreNonlinearProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {
    assert( dimDomain == 2 );
  }

  double Cos( const double a ) const { return cos( a ); }
  double Sin( const double a ) const { return sin( a ); }
  double Power( const double a, int n ) const
  {
    if( n == 1 )
      return a;
    else
      return a * Power( a, n-1 );
  }
  double Power( const double a, const double n ) const
  {
    return std::pow( a, n );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x,
		 RangeType& phi) const
  {
    const double t = time();
    const double xx = x[0];
    const double yy = x[1];
    const double Y1 = Ys_.at(0);
    const double Y2 = Ys_.at(1);
    const double p4 = 4.0 * M_PI;

    phi = (4*Cos(t)*(Cos(3*xx) + Y1*Cos(5*xx) + Cos(3*yy) + Y2*Cos(5*yy)) - 
     ((Cos(3*xx) + Y1*Cos(5*xx) + Cos(3*yy) + Y2*Cos(5*yy))*Sin(t)*
        (8*Sin(t) - 8*Power(yy,2)*(Cos(t) + Sin(t)) + Sin(2*t)))/
      (Power(4 + Cos(t),2) - 4*Power(yy,2)*(Cos(t) - Sin(t))) - 
     (2*xx*Cos(t)*Sin(t)*(3*Sin(3*xx) + 5*Y1*Sin(5*xx)))/(4 + Sin(t)) + 
     (2*yy*Power(Sin(t),2)*(3*Sin(3*yy) + 5*Y2*Sin(5*yy)))/(4 + Cos(t)) + 
     (Sin(t)*(8*p4*xx*Power(yy,3)*Cos(p4*(xx*Y1 + Y2*yy))*
           (1 + 8*Cos(t) + Cos(2*t) - 8*Sin(t) - Sin(2*t))*
           (3*Y2*Sin(3*xx) + Y1*(5*Y2*Sin(5*xx) + 3*Sin(3*yy) + 5*Y2*Sin(5*yy)))\
           - 4*yy*Power(4 + Cos(t),2)*
           (p4*xx*(4 + Cos(t))*Cos(p4*(xx*Y1 + Y2*yy))*
              (3*Y2*Sin(3*xx) + Y1*
                 (5*Y2*Sin(5*xx) + 3*Sin(3*yy) + 5*Y2*Sin(5*yy))) + 
             (4 + Sin(t))*(3*Sin(3*yy) + 5*Y2*Sin(5*yy))*
              (4 + Sin(p4*(xx*Y1 + Y2*yy)))) + 
          4*Power(yy,2)*Power(4 + Cos(t),2)*
           (-144*Cos(3*yy) - 400*Y2*Cos(5*yy) + 
             12*p4*Y1*Cos(p4*(xx*Y1 + Y2*yy))*Sin(3*xx) + 
             20*p4*Power(Y1,2)*Cos(p4*(xx*Y1 + Y2*yy))*Sin(5*xx) - 
             12*p4*Y2*Cos(p4*(xx*Y1 + Y2*yy))*Sin(3*yy) - 
             20*p4*Power(Y2,2)*Cos(p4*(xx*Y1 + Y2*yy))*Sin(5*yy) - 
             36*Cos(3*yy)*Sin(p4*(xx*Y1 + Y2*yy)) - 
             100*Y2*Cos(5*yy)*Sin(p4*(xx*Y1 + Y2*yy)) + 
             9*Cos(3*xx)*(4 + Sin(t))*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
             25*Y1*Cos(5*xx)*(4 + Sin(t))*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
             Sin(t)*(p4*Cos(p4*(xx*Y1 + Y2*yy))*
                 (3*Y1*Sin(3*xx) + 5*Power(Y1,2)*Sin(5*xx) + 3*Y2*Sin(3*yy) + 
                   5*Power(Y2,2)*Sin(5*yy)) + 
                9*Cos(3*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
                25*Y2*Cos(5*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy)))) - 
             2*Cos(t)*(9*Cos(3*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
                Y2*(p4*Cos(p4*(xx*Y1 + Y2*yy))*(3*Sin(3*yy) + 5*Y2*Sin(5*yy)) + 
                   25*Cos(5*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy)))))) + 
          16*Power(yy,4)*(Cos(t) - Sin(t))*
           (144*Cos(3*yy) + 400*Y2*Cos(5*yy) - 
             12*p4*Y1*Cos(p4*(xx*Y1 + Y2*yy))*Sin(3*xx) - 
             20*p4*Power(Y1,2)*Cos(p4*(xx*Y1 + Y2*yy))*Sin(5*xx) - 
             p4*Y1*Cos(p4*(xx*Y1 + Y2*yy))*Sin(t)*
              (3*Sin(3*xx) + 5*Y1*Sin(5*xx)) + 
             12*p4*Y2*Cos(p4*(xx*Y1 + Y2*yy))*Sin(3*yy) + 
             20*p4*Power(Y2,2)*Cos(p4*(xx*Y1 + Y2*yy))*Sin(5*yy) + 
             36*Cos(3*yy)*Sin(p4*(xx*Y1 + Y2*yy)) + 
             100*Y2*Cos(5*yy)*Sin(p4*(xx*Y1 + Y2*yy)) - 
             9*Cos(3*xx)*(4 + Sin(t))*(4 + Sin(p4*(xx*Y1 + Y2*yy))) - 
             25*Y1*Cos(5*xx)*(4 + Sin(t))*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
             Cos(t)*(9*Cos(3*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
                Y2*(p4*Cos(p4*(xx*Y1 + Y2*yy))*(3*Sin(3*yy) + 5*Y2*Sin(5*yy)) + 
                   25*Cos(5*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy)))))) + 
          Power(4 + Cos(t),3)*(-4*xx*
              (3 + 5*Y1 + 2*(3 + 5*Y1)*Cos(2*xx) + 10*Y1*Cos(4*xx))*Sin(xx)*
              (4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
             (4 + Cos(t))*(9*Cos(3*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy))) + 
                Y2*(p4*Cos(p4*(xx*Y1 + Y2*yy))*(3*Sin(3*yy) + 5*Y2*Sin(5*yy)) + 
                   25*Cos(5*yy)*(4 + Sin(p4*(xx*Y1 + Y2*yy))))))))/
	   Power(Power(4 + Cos(t),2) - 4*Power(yy,2)*(Cos(t) - Sin(t)),2))/4;
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

    D *= 1 + 0.25 * sin( 4.0 * M_PI * Ys_.at(0) * x[0] + 4.0 * M_PI * Ys_.at(1) * x[1] );
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
    phi = sin( time() ) * ( cos( 3*x[0] ) + cos( 3 * x[1] ) );
    // phi += sin( time() ) * ( Ys_.at(0) * cos( 5 * x[0] ) + Ys_.at(1) * cos( 5 * x[1] ) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * ( -3*sin( 3*x[0] ) );
    grad[ 0 ][ 1 ] = sin( time() ) * ( -3*sin( 3*x[1] ) );
    // grad[ 0 ][ 0 ] += Ys_.at(0) * sin( time() ) * ( -5*sin( 5*x[0] ) );
    // grad[ 0 ][ 1 ] += Y2_ * sin( time() ) * ( -5*sin( 5*x[1] ) );

    DomainType nu = x;

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

  virtual void setYs( const std::vector<double> Ys )
  {
    Ys_ = Ys;
  }

private:
  std::vector<double> Ys_;
};


template <class FunctionSpace>
class CurveMCNonUniformProblem : public TemporalProblemInterface < FunctionSpace >
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

  CurveMCNonUniformProblem( const Dune::Fem::TimeProviderBase &timeProvider )
    : BaseType( timeProvider )
  {
    assert( dimDomain == 2 );
  }

  double Cos( const double a ) const { return cos( a ); }
  double Sin( const double a ) const { return sin( a ); }
  double Power( const double a, int n ) const
  {
    if( n == 1 )
      return a;
    else
      return a * Power( a, n-1 );
  }

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& z,
		 RangeType& phi) const
  {
    throw "not implemented yet";
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

    D *= ( 1 + exp( -x[0]*x[0] / ( 0.5 * Ys_.at(0) + 0.5 ) )
	   + exp( -x[1]*x[1] / ( 0.5 * Ys_.at(1) + 0.5 ) ) );
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
    phi = sin( time() ) * ( cos( 3*x[0] ) + cos( 3 * x[1] ) );
    // + sin( time() ) * ( Ys_.at(0) * cos( 5 * x[0] ) + Ys_.at(1) * cos( 5 * x[1] ) );
  }

  //! the jacobian of the exact solution
  virtual void uJacobian(const DomainType& x,
			 JacobianRangeType& ret) const
  {
    JacobianRangeType grad;
    grad[ 0 ][ 0 ] = sin( time() ) * ( -3*sin( 3*x[0] ) );
    grad[ 0 ][ 1 ] = sin( time() ) * ( -3*sin( 3*x[1] ) );

    DomainType nu = x;

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

  virtual void setYs( const std::vector<double> Ys )
  {
    Ys_ = Ys;
  }

private:
  std::vector<double> Ys_;
};


#endif // #ifndef POISSON_HH
