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
#ifndef ELLIPT_MEANSCHEME_HH
#define ELLIPT_MEANSCHEME_HH

// iostream includes
#include <iostream>

// include discrete function space
#include <dune/fem/space/lagrange.hh>

// adaptation ...
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/common/adaptmanager.hh>

// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>

// include linear operators
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/diagonalpreconditioner.hh>

#include <dune/fem/operator/linear/istloperator.hh>
#include <dune/fem/solver/istlsolver.hh>
#include <dune/fem/solver/cginverseoperator.hh>

// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

/*********************************************************/

// include norms
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include parameter handling
#include <dune/fem/io/parameter.hh>

// local includes
#include "probleminterface.hh"

#include "model.hh"

#include "rhs.hh"
#include "elliptic.hh"

#include "femscheme.hh"

// FemScheme
//----------

/*******************************************************************************
 * template arguments are:
 * - GridPsrt: the part of the grid used to tesselate the
 *             computational domain
 * - Model: description of the data functions and methods required for the
 *          elliptic operator (massFlux, diffusionFlux)
 *     Model::ProblemType boundary data, exact solution,
 *                        and the type of the function space
 *******************************************************************************/
template < class Scheme >
class MeanScheme
{
public:
  typedef Scheme SchemeType;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename SchemeType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename SchemeType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef typename SchemeType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
  typedef typename SchemeType :: DiscreteFunctionType DiscreteFunctionType;

  // norm types
  typedef Dune::Fem::L2Norm< GridPartType > L2NormType;
  typedef Dune::Fem::H1Norm< GridPartType > H1NormType;

  MeanScheme( GridPartType &gridPart )
    : gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      linftyl2Error_( 0 ),
      l2h1Error_( 0 )
  {
    // set all DoF to zero
    solution_.clear();
  }

  void computeMean( const std::vector< SchemeType >& schemeVector )
  {
    solution_.clear();

    for( const auto& scheme : schemeVector )
      {
	solution_ += scheme.solution();
      }

    solution_ /= (double)schemeVector.size();
  }

  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  template< class GridExactSolution >
  void closeTimestep( const GridExactSolution &exact, const double deltaT )
  {
    // find l2 error
    L2NormType l2norm( gridPart_ );
    const double l2Error = l2norm.distance( exact, solution() );
    linftyl2Error_ = std::max( linftyl2Error_, l2Error );

    // find h1 error
    H1NormType h1norm( gridPart_ );
    const double h1Error = h1norm.distance( exact, solution() );
    l2h1Error_ = std::sqrt( l2h1Error_*l2h1Error_ + deltaT * h1Error * h1Error );
  }

  double linftyl2Error() const
  {
    return linftyl2Error_;
  }

  double l2h1Error() const
  {
    return l2h1Error_;
  }

  const int dofs() const
  {
    int tmp = discreteSpace_.size();
    return Dune::Fem::MPIManager::comm().sum( tmp );
  }

protected:
  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown

  double linftyl2Error_;
  double l2h1Error_;
};

#endif // end #if ELLIPT_MEANSCHEME_HH
