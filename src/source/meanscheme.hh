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
#ifndef ELLIPT_FEMSCHEME_HH
#define ELLIPT_FEMSCHEME_HH

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

// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "poisson-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};

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
template < class Model >
class MeanScheme
{
public:
  //! type of the mathematical model
  typedef Model ModelType ;

  //! grid view (e.g. leaf grid view) provided in the template argument list
  typedef typename ModelType::GridPartType GridPartType;

  //! type of underyling hierarchical grid needed for data output
  typedef typename GridPartType::GridType GridType;

  //! type of function space (scalar functions, \f$ f: \Omega -> R) \f$
  typedef typename ModelType :: FunctionSpaceType   FunctionSpaceType;

  //! choose type of discrete function space
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;

  // choose type of discrete function, Matrix implementation and solver implementation
#if HAVE_DUNE_ISTL && WANT_ISTL
  typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::ISTLLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  // typedef Dune::Fem::ISTLCGOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
  typedef Dune::Fem::ISTLGMResOp< DiscreteFunctionType, LinearOperatorType > LinearInverseOperatorType;
#else
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
  typedef Dune::Fem::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;
#endif

  /*********************************************************/

  //! define Laplace operator
  typedef DifferentiableEllipticOperator< LinearOperatorType, ModelType > EllipticOperatorType;

  typedef FemScheme< ModelType > FemSchemeType;

  MeanScheme( GridPartType &gridPart,
	      const vector< FemSchemeType >& schemeVector )
    : gridPart_( gridPart ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      schemeVector_( schemeVector ),
      linftyl2Error_( 0 ),
      l2h1Error_( 0 )
  {
    // set all DoF to zero
    solution_.clear();
  }

  const DiscreteFunctionType &solution() const
  {
    solution_.clear();

    for( auto scheme: schemeVector )
      solution += scheme.solution();

    solution /= (double)schemeVector.size();

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

  const int dofs() const
  {
    int tmp = discreteSpace_.size();
    return Dune::Fem::MPIManager::comm().sum( tmp );
  }

protected:
  GridPartType  &gridPart_;         // grid part(view), e.g. here the leaf grid the discrete space is build with

  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown

  const vector< FemSchemeType >& schemeVector;

  double linftyl2Error_;
  double l2h1Error_;
};

#endif // end #if ELLIPT_FEMSCHEME_HH
