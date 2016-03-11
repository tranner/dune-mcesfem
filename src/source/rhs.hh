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
#ifndef RHS_HH
#define RHS_HH

#include <dune/fem/quadrature/cachingquadrature.hh>


// assembleRHS
// -----------

template< class Function, class DiscreteFunction >
void assembleRHS ( const Function &function, DiscreteFunction &rhs )
{
  typedef typename DiscreteFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunction::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;
  typedef typename EntityType::Geometry GeometryType;

  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef Dune::Fem::CachingQuadrature< GridPartType, 0 > QuadratureType;

  rhs.clear();

  const DiscreteFunctionSpaceType &dfSpace = rhs.space();

  const IteratorType end = dfSpace.end();
  for( IteratorType it = dfSpace.begin(); it != end; ++it )
  {
    const EntityType &entity = *it;
    const GeometryType &geometry = entity.geometry();

    typename Function::LocalFunctionType localFunction = 
             function.localFunction( entity);
    LocalFunctionType rhsLocal = rhs.localFunction( entity );
    
    QuadratureType quadrature( entity, 2*dfSpace.order()+1 );
    const size_t numQuadraturePoints = quadrature.nop();
    for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
    {
      // obtain quadrature point
      const typename QuadratureType::CoordinateType &x = quadrature.point( pt );

      // evaluate f
      typename Function::RangeType f;
      localFunction.evaluate( quadrature[ pt ], f );

      // multiply by quadrature weight
      f *= quadrature.weight( pt ) * geometry.integrationElement( x );

      // add f * phi_i to rhsLocal[ i ]
      rhsLocal.axpy( quadrature[ pt ], f );
    }
  }
  rhs.communicate();
}

#endif // #ifndef RHS_HH
