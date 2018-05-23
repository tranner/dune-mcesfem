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
#ifndef POISSON_PROBLEMINTERFACE_HH
#define POISSON_PROBLEMINTERFACE_HH

#include <cassert>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/fem/function/common/function.hh>

/** \brief problem interface class for problem descriptions, i.e. right hand side,
 *         boudnary data, and, if exsistent, an exact solution.
 */
template <class FunctionSpace>
class ProblemInterface : public Dune::Fem::Function< FunctionSpace, ProblemInterface<FunctionSpace> >
{
public:  
  // type of function space 
  typedef FunctionSpace  FunctionSpaceType;  

  enum { dimRange  = FunctionSpaceType :: dimRange  };
  enum { dimDomain = FunctionSpaceType :: dimDomain };

  typedef typename FunctionSpaceType :: RangeFieldType   RangeFieldType;

  typedef typename FunctionSpaceType :: RangeType   RangeType;
  typedef typename FunctionSpaceType :: DomainType  DomainType;

  typedef typename FunctionSpaceType :: JacobianRangeType  JacobianRangeType;

  typedef Dune::FieldMatrix< RangeFieldType, dimDomain, dimDomain > DiffusionTensorType;
  typedef Dune::FieldVector< RangeFieldType, dimDomain > AdvectionVectorType;

  //! the right hand side data (default = 0)
  virtual void f(const DomainType& x, 
                 RangeType& value) const
  {
    value = 0; 
  }

  virtual void boundaryRhs( const DomainType& x,
			    RangeType& value ) const
  {
    value = 0;
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

  //! robin coefficient (default = 0)
  virtual void a(const DomainType& x, RangeType &a) const
  {
    a = RangeType(0);
  }

  //! capacity coefficient (default = 1)
  virtual void d(const DomainType& x, RangeType &d) const
  {
    d = RangeType(0);
  }

  //! the exact solution (default = 0)
  virtual void u(const DomainType& x, 
                 RangeType& value) const 
  {
    value = 0; 
  }

  //! the jacobian of the exact solution (default = 0)
  virtual void uJacobian(const DomainType& x, 
                         JacobianRangeType& value) const 
  {
    value = 0; 
  }

  //! return true if Dirichlet boundary is present (default is true)
  virtual bool hasDirichletBoundary () const 
  {
    return false ;
  }

  //! return true if given point belongs to the Dirichlet boundary (default is true)
  virtual bool isDirichletPoint( const DomainType& x ) const 
  {
    return true ;
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  virtual bool isNeumannPoint( const DomainType& x ) const
  {
    return false ;
  }

  //! the Dirichlet boundary data (default calls u)
  virtual void g(const DomainType& x, 
                 RangeType& value) const 
  {
    u( x, value );
  }

  //! make this into a fem function for the exact solution
  void evaluate( const DomainType& x, RangeType& ret ) const 
  {
    // call exact solution of implementation 
    u( x, ret );
  }
  //! also need the jacobian of the exact solution
  void jacobian( const DomainType& x, JacobianRangeType& jac ) const 
  {
    uJacobian( x, jac );
  }

  virtual void setYs( const std::vector<double> Ys )
  {
    assert(0);
  }
};

#endif // #ifndef ELLIPTC_PROBLEMINTERFACE_HH

