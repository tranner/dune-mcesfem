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
#ifndef ELLIPTC_MODEL_HH
#define ELLIPTC_MODEL_HH

#include <cassert>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/parameter.hh>

#include "probleminterface.hh"

/********************************************
  Full non-linar elliptic model

  struct EllipticModel
  {
    //! [Methods used for operator application]
    template< class Entity, class Point >
    void source ( const Entity &entity, const Point &x,
                  const RangeType &value, 
                  RangeType &flux ) const;
    template< class Entity, class Point >
    void diffusiveFlux ( const Entity &entity, const Point &x,
                         const RangeType &value, const JacobianRangeType &gradient,
                         JacobianRangeType &flux ) const;
    //! [Methods used for operator application]

    //! [Methods used to assemble linearized operator]
    template< class Entity, class Point >
    void linSource ( const RangeType& uBar, 
                     const Entity &entity, const Point &x,
                     const RangeType &value, 
                     RangeType &flux ) const;
    template< class Entity, class Point >
    void linDiffusiveFlux ( const RangeType& uBar, const JacobianRangeType& gradientBar,
                            const Entity &entity, const Point &x,
                            const RangeType &value, const JacobianRangeType &gradient,
                            JacobianRangeType &flux ) const;
    //! [Methods used to assemble linearized operator]

    //! [Methods used for Dirichlet constraints]
    bool hasDirichletBoundary () const;
    template <class Intersection>
    bool isDirichletIntersection( const Intersection& inter ) const 
    bool isDirichletPoint( const DomainType& x ) const;
    template< class Entity, class Point >
    void g( const RangeType& uBar, 
            const Entity &entity, const Point &x,
            RangeType &u ) const;
    //! [Methods used for Dirichlet constraints]
  };

 ********************************************/

// DiffusionModel
// --------------

template< class FunctionSpace, class GridPart >
struct DiffusionModel
{
  typedef FunctionSpace FunctionSpaceType;
  typedef GridPart GridPartType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef ProblemInterface< FunctionSpaceType > ProblemType ;

  static const bool isLinear = true;
  static const bool isSymmetric = true;

protected:
  enum FunctionId { rhs, bnd };
  template <FunctionId id>
  class FunctionWrapper;
public:
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<rhs>, GridPartType > RightHandSideType;
  typedef Dune::Fem::GridFunctionAdapter< FunctionWrapper<bnd>, GridPartType > DirichletBoundaryType;

  //! constructor taking problem reference 
  DiffusionModel( const ProblemType& problem, const GridPart &gridPart )
    : problem_( problem ),
      gridPart_(gridPart),
      rhs_(problem_),
      bnd_(problem_)
  {
  }

  template< class Entity, class Point >
  void source ( const Entity &entity, 
                const Point &x,
                const RangeType &value, 
                RangeType &flux ) const
  {
    linSource( value, entity, x, value, flux );
  }

  // the linearization of the source function
  template< class Entity, class Point >
  void linSource ( const RangeType& uBar, 
                   const Entity &entity, 
                   const Point &x,
                   const RangeType &value, 
                   RangeType &flux ) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    RangeType m;
    problem_.m(xGlobal,m);
    for (unsigned int i=0;i<flux.size();++i)
      flux[i] = m[i]*value[i];
  }
  //! return the diffusive flux 
  template< class Entity, class Point >
  void diffusiveFlux ( const Entity &entity, 
                       const Point &x,
                       const RangeType &value,
                       const JacobianRangeType &gradient,
                       JacobianRangeType &flux ) const
  {
    linDiffusiveFlux( value, gradient, entity, x, value, gradient, flux );
  }
  // linearization of diffusiveFlux
  template< class Entity, class Point >
  void linDiffusiveFlux ( const RangeType& uBar, 
                          const JacobianRangeType& gradientBar,
                          const Entity &entity, 
                          const Point &x,
                          const RangeType &value,
                          const JacobianRangeType &gradient,
                          JacobianRangeType &flux ) const
  {
    // the flux is simply the identity 
    flux = gradient;
  }

  //! exact some methods from the problem class
  bool hasDirichletBoundary () const 
  {
    return problem_.hasDirichletBoundary() ;
  }

  //! return true if given intersection belongs to the Dirichlet boundary -
  //! we test here if the center is a dirichlet point
  template <class Intersection>
  bool isDirichletIntersection( const Intersection& inter ) const 
  {
    return isDirichletPoint( inter.geometry().center() );
  }


  //! return true if given point belongs to the Dirichlet boundary (default is true)
  bool isDirichletPoint( const DomainType& x ) const 
  {
    return problem_.isDirichletPoint(x) ;
  }

  //! return true if given intersection belongs to the Neumann boundary -
  //! we test here if the center is a neumann point
  template <class Intersection>
  bool isNeumannIntersection( const Intersection& inter ) const
  {
    return isNeumannPoint( inter.geometry().center() );
  }

  //! return true if given point belongs to the Neumann boundary (default is false)
  bool isNeumannPoint( const DomainType& x ) const
  {
    return problem_.isNeumannPoint(x) ;
  }

  template< class Entity, class Point >
  void g( const RangeType& uBar, 
          const Entity &entity, 
          const Point &x,
          RangeType &u ) const
  {
    const DomainType xGlobal = entity.geometry().global( coordinate( x ) );
    problem_.g( xGlobal, u );
  }

  // return Fem :: Function for Dirichlet boundary values 
  DirichletBoundaryType dirichletBoundary( ) const 
  {
    return DirichletBoundaryType( "boundary function", bnd_, gridPart_, 5 );  
  }

  // return Fem :: Function for right hand side 
  RightHandSideType rightHandSide(  ) const 
  {
    return RightHandSideType( "right hand side", rhs_, gridPart_, 5 );  
  }
  
protected:
  template <FunctionId id>
  class FunctionWrapper : public Dune::Fem::Function< FunctionSpaceType, FunctionWrapper< id > >
  {
    const ProblemInterface<FunctionSpaceType>& impl_;
    public:   
    FunctionWrapper( const ProblemInterface<FunctionSpaceType>& impl )
    : impl_( impl ) {}
 
    //! evaluate function 
    void evaluate( const DomainType& x, RangeType& ret ) const 
    {
      if( id == rhs ) 
      {
        // call right hand side of implementation 
        impl_.f( x, ret );
      }
      else if( id == bnd ) 
      {
        // call dirichlet boudary data of implementation 
        impl_.g( x, ret );
      }
      else 
      {
        DUNE_THROW(Dune::NotImplemented,"FunctionId not implemented"); 
      }
    }
  };
   
  const ProblemType& problem_;
  const GridPart &gridPart_;
  FunctionWrapper<rhs> rhs_;
  FunctionWrapper<bnd> bnd_;
};

#endif // #ifndef ELLIPTC_MODEL_HH
