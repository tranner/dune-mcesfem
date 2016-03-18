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
#ifndef DEFORMATION_HH
#define DEFORMATION_HH

// include discrete function space
#include <dune/fem/space/lagrange.hh>
// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>
// lagrange interpolation 
#include <dune/fem/operator/lagrangeinterpolation.hh>

// DeformationCoordFunction
// ------------------------

template< class Deformation, class BoundaryProjection, class GridPart,
	  const unsigned int codim, const unsigned int polorder = 1 >
class DiscreteDeformationCoordHolder;

template< class Deformation, class BoundaryProjection, class GridPart, const unsigned int polorder >
class DiscreteDeformationCoordHolder< Deformation, BoundaryProjection, GridPart, 0, polorder >
{
  typedef Deformation DeformationType;
  typedef BoundaryProjection BoundaryProjectionType;
  typedef GridPart GridPartType;

public:
  typedef typename DeformationType :: FunctionSpaceType FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
  // choose type of discrete function
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  DiscreteDeformationCoordHolder( Deformation& deformation, BoundaryProjectionType& boundaryProjection,
				  GridPart& gridPart )
    : deformation_( deformation ),
      boundaryProjection_( boundaryProjection ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart ),
      coordFunction_( "deformation", discreteSpace_ )
  {
    interpolate();
  }

  const DiscreteFunctionType& coordFunction() const
  {
    return coordFunction_;
  }

  int dofs() const
  {
    return discreteSpace_.size();
  }

  int elements() const
  {
    int n = 0;
    for( auto it = discreteSpace_.begin(); it != discreteSpace_.end(); ++it )
      ++n;
    return n;
  }

  void interpolate()
  {
    typedef typename DiscreteFunctionType::DofType DofType;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    // set all DoFs to infinity
    const DofIteratorType dend = coordFunction_.dend();
    for( DofIteratorType dit = coordFunction_.dbegin(); dit != dend; ++dit )
      *dit = std::numeric_limits< DofType >::infinity();

    for( const auto& entity : discreteSpace_ )
      {
	const auto& lagrangePointSet
	  = discreteSpace_.lagrangePointSet( entity );
	const int nop = lagrangePointSet.nop();

	auto df_local = coordFunction_.localFunction( entity );

	// does element contain a boundary segment?
	const bool boundary = entity.hasBoundaryIntersections();

	// if not on boundary, map is identity
	if( not boundary )
	  {
	    const auto geometry = entity.geometry();

	    // assume point based local dofs
	    int k = 0;
	    for( int qp = 0; qp < nop; ++qp )
	      {
		// if the first DoF for this point is already valid, continue
		if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
		  {
		    // find local coordinate
		    const auto hatx = coordinate( lagrangePointSet[ qp ] );

		    // find global coordinate
		    typename DiscreteFunctionType::DomainType x = geometry.global( hatx );

		    // evaluate the function in the Lagrange point
		    typename DiscreteFunctionType::RangeType phi;
		    deformation_.evaluate( x, phi );

		    // assign the appropriate values to the DoFs
		    for( int i = 0; i < dimRange; ++i, ++k )
		      df_local[ k ] = phi[ i ];
		  }
		else
		  k += dimRange;
	      }
	  }
	else
	  {
	    const auto geometry = entity.geometry();
	    std::vector< bool > isVertexOnBoundary( entity.subEntities( dimRange ) );
	    for( int i = 0; i < geometry.corners(); ++i )
	      {
		isVertexOnBoundary[i] = boundaryProjection_.onBoundary( geometry.corner(i) );
	      }

	    // assume point based local dofs
	    const int nop = lagrangePointSet.nop();
	    int k = 0;
	    for( int qp = 0; qp < nop; ++qp )
	      {
		// if the first DoF for this point is already valid, continue
		if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
		  {
		    // find local coordinate
		    const auto hatx = coordinate( lagrangePointSet[ qp ] );

		    // find barycentric coordinates
		    Dune::FieldVector< double, dimRange+1 > lambda;
		    lambda[0] = 1.0;
		    for( int i = 0; i < dimRange; ++i )
		      {
			lambda[i+1] = hatx[i];
			lambda[0] -= hatx[i];
		      }

		    double lambdaStar = 0.0;
		    for( int i = 0; i < dimRange+1; ++i )
		      {
			if( isVertexOnBoundary[i] )
			  {
			    lambdaStar += lambda[i];
			  }
		      }

		    // find global coordinate
		    typename DiscreteFunctionType::DomainType x = geometry.global( hatx );

		    // find perturbation
		    if( lambdaStar > 1.0e-8 )
		      {
			Dune::FieldVector< double, dimRange > y(0);
			for( int i = 0; i < dimRange+1; ++i )
			  {
			    if( isVertexOnBoundary[i] )
			      y.axpy(  lambda[i] / lambdaStar, geometry.corner(i) );
			  }

			Dune::FieldVector< double, dimRange > p;
			boundaryProjection_.evaluate( y, p );

			x.axpy( std::pow( lambdaStar, df_local.order()+2 ), ( p - y ) );
		      }

		    // evaluate deformation
		    typename DiscreteFunctionType::RangeType phi;
		    deformation_.evaluate( x, phi );

		    // assign the appropriate values to the DoFs
		    for( int i = 0; i < dimRange; ++i, ++k )
		      df_local[ k ] = phi[ i ];
		  }
		else
		  k += dimRange;
	      }
	  }
      }
  }

private:
  Deformation& deformation_;
  BoundaryProjectionType& boundaryProjection_;
  GridPart& gridPart_;
  DiscreteFunctionSpaceType discreteSpace_;
  DiscreteFunctionType coordFunction_;
};

template< class Deformation, class BoundaryProjection, class GridPart, const unsigned int polorder >
class DiscreteDeformationCoordHolder< Deformation, BoundaryProjection, GridPart, 1, polorder >
{
  typedef Deformation DeformationType;
  typedef BoundaryProjection BoundaryProjectionType;
  typedef GridPart GridPartType;

public:
  typedef typename DeformationType :: FunctionSpaceType FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polorder > DiscreteFunctionSpaceType;
  // choose type of discrete function
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  DiscreteDeformationCoordHolder( Deformation& deformation, BoundaryProjectionType& boundaryProjection, GridPart& gridPart )
    : deformation_( deformation ),
      boundaryProjection_( boundaryProjection ),
      gridPart_( gridPart ),
      discreteSpace_( gridPart ),
      coordFunction_( "deformation", discreteSpace_ )
  {
    interpolate();
  }

  const DiscreteFunctionType& coordFunction() const
  {
    return coordFunction_;
  }

  int dofs() const
  {
    return discreteSpace_.size();
  }

  int elements() const
  {
    int n = 0;
    for( auto it = discreteSpace_.begin(); it != discreteSpace_.end(); ++it )
      ++n;
    return n;
  }

  void interpolate()
  {
    typedef typename DiscreteFunctionType::DofType DofType;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    static const int dimRange = DiscreteFunctionSpaceType::dimRange;

    // set all DoFs to infinity
    const DofIteratorType dend = coordFunction_.dend();
    for( DofIteratorType dit = coordFunction_.dbegin(); dit != dend; ++dit )
      *dit = std::numeric_limits< DofType >::infinity();

    for( const typename DiscreteFunctionSpaceType::EntityType &entity : discreteSpace_ )
      {
	const auto geometry = entity.geometry();

        const auto &lagrangePointSet
          = discreteSpace_.lagrangePointSet( entity );

        auto df_local = coordFunction_.localFunction( entity );

        // assume point based local dofs
        const int nop = lagrangePointSet.nop();
        int k = 0;
        for( int qp = 0; qp < nop; ++qp )
	  {
	    // if the first DoF for this point is already valid, continue
	    if( df_local[ k ] == std::numeric_limits< DofType >::infinity() )
	      {
		// evaluate the function in the Lagrange point
		const auto x = geometry.global( coordinate( lagrangePointSet[qp] ) );
		typename BoundaryProjectionType::DomainType p;
		boundaryProjection_.evaluate( x, p );

		typename Deformation::RangeType phi;
		deformation_.evaluate( p, phi );

		// assign the appropriate values to the DoFs
		for( int i = 0; i < dimRange; ++i, ++k )
		  df_local[ k ] = phi[ i ];
	      }
	    else
	      k += dimRange;
	  }
      }
  }

private:
  Deformation& deformation_;
  BoundaryProjectionType& boundaryProjection_;
  GridPart& gridPart_;
  DiscreteFunctionSpaceType discreteSpace_;
  DiscreteFunctionType coordFunction_;
};

#endif // #ifndef DEFORMATION_HH
