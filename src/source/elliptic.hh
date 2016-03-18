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
#ifndef ELLIPTIC_HH
#define ELLIPTIC_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

// EllipticOperator
// ----------------

//! [Class for elliptic operator]
template< class DiscreteFunction, class Model >
struct EllipticOperator
: public virtual Dune::Fem::Operator< DiscreteFunction >
//! [Class for elliptic operator]
{
  typedef DiscreteFunction DiscreteFunctionType;
  typedef Model            ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType::DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  //! type of Dirichlet constraints
  template< class MyModel, class MySpace >
  struct Constraints
  {
    Constraints( const MyModel& m, const MySpace& s )
    {}

    template< class F1, class F2 >
    void operator()( const F1& f1, F2& f2 ) const
    {
      return;
    }

    template< class Op >
    void applyToOperator( Op& op ) const
    {
      return;
    }
  };
  typedef Constraints< Model, DiscreteFunctionSpaceType > ConstraintsType;

public:
  //! contructor
  EllipticOperator ( ModelType &model, const DiscreteFunctionSpaceType &space )
  : model_( model )
  , constraints_( model, space )
  {}

  // prepare the solution vector
  template <class Function>
  void prepare( const Function &func, DiscreteFunctionType &u )
  {
    // set boundary values for solution
    constraints()( func, u );
  }

  //! application operator
  virtual void
  operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const;

  ModelType &model () { return model_; }
protected:
  const ModelType &model () const { return model_; }
  const ConstraintsType &constraints () const { return constraints_; }

private:
  ModelType model_;
  const ConstraintsType constraints_;
};

// DifferentiableEllipticOperator
// ------------------------------

//! [Class for linearizable elliptic operator]
template< class JacobianOperator, class Model >
struct DifferentiableEllipticOperator
: public EllipticOperator< typename JacobianOperator::DomainFunctionType, Model >,
  public Dune::Fem::DifferentiableOperator< JacobianOperator >
//! [Class for linearizable elliptic operator]
{
  typedef EllipticOperator< typename JacobianOperator::DomainFunctionType, Model > BaseType;

  typedef JacobianOperator JacobianOperatorType;

  typedef typename BaseType::DiscreteFunctionType DiscreteFunctionType;
  typedef typename BaseType::ModelType ModelType;

protected:
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
  typedef typename LocalFunctionType::RangeType RangeType;
  typedef typename LocalFunctionType::JacobianRangeType JacobianRangeType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity       EntityType;
  typedef typename EntityType::Geometry       GeometryType;

  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

  typedef typename DiscreteFunctionSpaceType::GridPartType  GridPartType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  typedef typename IntersectionIteratorType::Intersection IntersectionType;
  typedef typename IntersectionType::Geometry  IntersectionGeometryType;

  typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef typename BaseType::QuadratureType QuadratureType;

public:
  //! constructor
  DifferentiableEllipticOperator ( ModelType &model, const DiscreteFunctionSpaceType &space, bool sw=true )
  : BaseType( model, space )
  {}

  //! method to setup the jacobian of the operator for storage in a matrix
  void jacobian ( const DiscreteFunctionType &u, JacobianOperatorType &jOp ) const;

  using BaseType::model;
protected:
  using BaseType::constraints;
};

// Implementation of EllipticOperator
// ----------------------------------

template< class DiscreteFunction, class Model >
void EllipticOperator< DiscreteFunction, Model >
  ::operator() ( const DiscreteFunctionType &u, DiscreteFunctionType &w ) const
{
  // clear destination
  w.clear();

  // get discrete function space
  const DiscreteFunctionSpaceType &dfSpace = w.space();

  // iterate over grid
  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    // get entity (here element)
    const EntityType &entity = *it;
    // get elements geometry
    const GeometryType &geometry = entity.geometry();

    // get local representation of the discrete functions
    const LocalFunctionType uLocal = u.localFunction( entity );
    LocalFunctionType wLocal = w.localFunction( entity );

    // obtain quadrature order
    const int quadOrder = 2*( uLocal.order() + wLocal.order() );

    { // element integral
      QuadratureType quadrature( entity, quadOrder );
      const size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
      {
        //! [Compute local contribution of operator]
        const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
        const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

        RangeType vu;
        uLocal.evaluate( quadrature[ pt ], vu );
        JacobianRangeType du;
        uLocal.jacobian( quadrature[ pt ], du );

        // compute mass contribution (studying linear case so linearizing around zero)
        RangeType avu( 0 );
        model().source( entity, x, vu, avu );
        avu *= weight;
        // add to local functional wLocal.axpy( quadrature[ pt ], avu );

        JacobianRangeType adu( 0 );
        // apply diffusive flux
        model().diffusiveFlux( entity, x, vu, du, adu );
        adu *= weight;

        // add to local function
        wLocal.axpy( quadrature[ pt ], avu, adu );
        //! [Compute local contribution of operator]
      }

      // boundary integral
      const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
      for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity );
	     iit != iitend; ++iit ) // iterate over intersections
	{
	  const IntersectionType &intersection = *iit;

	  if( intersection.boundary() )
	    {
	      if( not model().isNeumannIntersection( intersection ) )
		{
		  continue;
		}

	      // extract intersection geometry
	      const auto& intersectionGeometry = intersection.geometry();

	      // find quadrature rule on intersection
	      FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, quadOrder, FaceQuadratureType :: INSIDE );
	      const size_t numQuadraturePoints = quadInside.nop();

	      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
		{
		  // find quadrature
		  const typename FaceQuadratureType::LocalCoordinateType &xLocal = quadInside.localPoint( pt );
		  double weight = quadInside.weight( pt ) * intersectionGeometry.integrationElement( xLocal );

		  // evaluate boundary flux
		  RangeType vu, avu;
		  uLocal.evaluate( quadInside[ pt ], vu );
		  model().boundaryFlux( entity, coordinate( quadInside[ pt ] ), vu, avu );

		  // multiply by quadrature weight
		  avu *= weight;

		  // add to local function
		  wLocal.axpy( quadInside[ pt ], avu );
		}
	    }
	}
    }
  }

  // communicate data (in parallel runs)
  w.communicate();

  // apply constraints, e.g. Dirichlet contraints, to the result
  constraints()( u, w );
}

// Implementation of DifferentiableEllipticOperator
// ------------------------------------------------

template< class JacobianOperator, class Model >
void DifferentiableEllipticOperator< JacobianOperator, Model >
  ::jacobian ( const DiscreteFunctionType &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;
  typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;

  const DiscreteFunctionSpaceType &dfSpace = u.space();

  Dune::Fem::DiagonalStencil<DiscreteFunctionSpaceType,DiscreteFunctionSpaceType> stencil(dfSpace,dfSpace);
  jOp.reserve(stencil);
  jOp.clear();

  const int blockSize = dfSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename LocalFunctionType::RangeType > phi( dfSpace.blockMapper().maxNumDofs()*blockSize );
  std::vector< typename LocalFunctionType::JacobianRangeType > dphi( dfSpace.blockMapper().maxNumDofs()*blockSize );

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    const LocalFunctionType uLocal = u.localFunction( entity );
    LocalMatrixType jLocal = jOp.localMatrix( entity, entity );

    const BasisFunctionSetType &baseSet = jLocal.domainBasisFunctionSet();
    const unsigned int numBasisFunctions = baseSet.size();

    QuadratureType quadrature( entity, 4*dfSpace.order() );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      //! [Assembling the local matrix]
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );
      const double weight = quadrature.weight( pt ) * geometry.integrationElement( x );

      // evaluate all basis functions at given quadrature point
      baseSet.evaluateAll( quadrature[ pt ], phi );

      // evaluate jacobians of all basis functions at given quadrature point
      baseSet.jacobianAll( quadrature[ pt ], dphi );

      // get value for linearization
      RangeType u0;
      JacobianRangeType jacU0;
      uLocal.evaluate( quadrature[ pt ], u0 );
      uLocal.jacobian( quadrature[ pt ], jacU0 );

      RangeType aphi( 0 );
      JacobianRangeType adphi( 0 );
      for( unsigned int localCol = 0; localCol < numBasisFunctions; ++localCol )
      {
        // if mass terms or right hand side is present
        model().linSource( u0, entity, x, phi[ localCol ], aphi );

        // if gradient term is present
        model().linDiffusiveFlux( u0, jacU0, entity, x, phi[ localCol ], dphi[ localCol ], adphi );

        // get column object and call axpy method
        jLocal.column( localCol ).axpy( phi, dphi, aphi, adphi, weight );
      }
      //! [Assembling the local matrix]
    }

      // boundary integral
      const IntersectionIteratorType iitend = dfSpace.gridPart().iend( entity );
      for( IntersectionIteratorType iit = dfSpace.gridPart().ibegin( entity );
	     iit != iitend; ++iit ) // iterate over intersections
	{
	  const IntersectionType &intersection = *iit;

	  if( intersection.boundary() )
	    {
	      if( not model().isNeumannIntersection( intersection ) )
		{
		  continue;
		}

	      // extract intersection geometry
	      const auto& intersectionGeometry = intersection.geometry();

	      // find quadrature rule on intersection
	      FaceQuadratureType quadInside( dfSpace.gridPart(), intersection, 2*dfSpace.order(),
					     FaceQuadratureType :: INSIDE );
	      const size_t numQuadraturePoints = quadInside.nop();

	      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
		{
		  // find quadrature
		  const typename FaceQuadratureType::LocalCoordinateType &xLocal = quadInside.localPoint( pt );
		  double weight = quadInside.weight( pt ) * intersectionGeometry.integrationElement( xLocal );

		  // evaluate basis functions at quadrature point
		  baseSet.evaluateAll( quadInside[ pt ], phi );

		  // get value for linearization
		  RangeType u0;
		  uLocal.evaluate( quadInside[ pt ], u0 );

		  RangeType aphi( 0 );
		  for( unsigned int localCol = 0; localCol < numBasisFunctions; ++localCol )
		    {
		      // if boundary flux term is present
		      model().linBoundaryFlux( u0, entity, coordinate( quadInside[ pt ] ), phi[ localCol ], aphi );

		      // get column object and call axpy method
		      jLocal.column( localCol ).axpy( phi, aphi, weight );
		    }
		}
	    }
	}
  } // end grid traversal

  // apply constraints to matrix operator
  constraints().applyToOperator( jOp );
  jOp.communicate();
}

#endif // #ifndef ELLIPTIC_HH
