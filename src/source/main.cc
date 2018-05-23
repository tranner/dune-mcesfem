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
#include <config.h>

// iostream includes
#include <iostream>

// include for rng
#include <random>

// dune grid includes
#include <dune/grid/albertagrid/agrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>

// include grid part
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

// include output
#include <dune/fem/io/file/dataoutput.hh>

// include eoc output
#include <dune/fem/misc/femeoc.hh>
#include "gridwidth.hh"

// include geometrty grid part
#include <dune/fem/gridpart/geogridpart.hh>

// include discrete function space
#include <dune/fem/space/lagrange.hh>
// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>
// lagrange interpolation
#include <dune/fem/operator/lagrangeinterpolation.hh>

#include "deformation.hh"

// include header for heat model
#include "heat.hh"

#include "heatmodel.hh"
#include "heatscheme.hh"

#include "meanscheme.hh"

#include <dune/fem/space/common/functionspace.hh>
template< int dimWorld >
struct BoundaryProjection
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
    y /= y.two_norm();
  }

  bool onBoundary( const DomainType& x )
  {
    return true;
  }
};

template< int dimWorld >
struct DeformationCoordFunctionBase
{
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  explicit DeformationCoordFunctionBase ( const double time = 0.0 )
  : time_( time )
  {}

  virtual void evaluate ( const DomainType &x, RangeType &y ) const
  {
    y = x;
  }

  double time() const { return time_; }
  void setTime ( const double time ) { time_ = time; }

private:
  double time_;
};

template< int dimWorld >
struct DeformationCoordFunction
  : public DeformationCoordFunctionBase< dimWorld >
{
  typedef DeformationCoordFunctionBase< dimWorld > BaseType;
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  explicit DeformationCoordFunction ( const double time = 0.0 )
    : BaseType( time )
  {}

  virtual void evaluate ( const DomainType &x, RangeType &y ) const
  {
    const double at = 1.0 + 0.25 * sin( time() );

    y[ 0 ] = x[ 0 ] * sqrt(at);
    y[ 1 ] = x[ 1 ];
    y[ 2 ] = x[ 2 ];
  }


private:
  using BaseType :: time;
};

template< int dimWorld >
struct DeformationCoordFunction2
  : public DeformationCoordFunctionBase< dimWorld >
{
  typedef DeformationCoordFunctionBase< dimWorld > BaseType;
  typedef Dune::Fem::FunctionSpace< double, double, dimWorld, dimWorld > FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;

  explicit DeformationCoordFunction2 ( const double time = 0.0 )
    : BaseType( time )
  {}

  virtual void evaluate ( const DomainType &x, RangeType &y ) const
  {
    const double at = 1.0 + 0.25 * sin( time() );
    const double bt = 1.0 + 0.25 * cos( time() );

    y[ 0 ] = x[ 0 ] * sqrt(at);
    y[ 1 ] = x[ 1 ] * sqrt(bt);
  }

private:
  using BaseType :: time;
};

// assemble-solve-estimate-mark-refine-IO-error-doitagain
template <class HGridType>
void algorithm ( HGridType &grid, int step, const int eocId )
{
  typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimensionworld, 1 > FunctionSpaceType;

  // create time provider
  Dune::Fem::GridTimeProvider< HGridType > timeProvider( grid );

  // we want to solve the problem on the leaf elements of the grid
  //! [Setup the grid part for a deforming domain]
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType, Dune::InteriorBorder_Partition > HostGridPartType;
  HostGridPartType hostGridPart( grid );

  // construct deformation
  typedef BoundaryProjection< HGridType::dimensionworld > BoundaryProjectionType;
  BoundaryProjectionType boundaryProjection;
  typedef DeformationCoordFunctionBase< HGridType::dimensionworld > DeformationType;

  // choose deformation
  const std::string problemNames [] = { "surface_heat", "surface_mc", "surface_stationary_mc", "curve_mc", "curve_mc_less-smooth", "curve_mc_more-nonlinear", "curve_mc_non-uniform" };
  const int problemNumber = Dune :: Fem :: Parameter :: getEnum( "heat.problem", problemNames );
  DeformationType *deformationPtr = 0;
  switch( problemNumber )
    {
    case 0:
      deformationPtr = new DeformationCoordFunction< HGridType::dimensionworld >;
      std::cout << "using ellipsoid with varying major axis" << std::endl;
      break;
    case 1:
      deformationPtr = new DeformationCoordFunction< HGridType::dimensionworld >;
      std::cout << "using ellipsoid with varying major axis" << std::endl;
      break;
    case 2:
      deformationPtr = new DeformationCoordFunctionBase< HGridType::dimensionworld >;
      std::cout << "using stationary surface" << std::endl;
    case 3:
    case 4:
    case 5:
    case 6:
      deformationPtr = new DeformationCoordFunction2< HGridType::dimensionworld >;
      std::cout << "using ellipse with both axes varying" << std::endl;
      break;

      break;

    default:
      std::cerr << "unrecognised problem name" << std::endl;
      assert(0);
    }

  // recover deformation
  assert( deformationPtr );
  DeformationType &deformation = *deformationPtr;

  typedef DiscreteDeformationCoordHolder< DeformationType, BoundaryProjectionType,
					  HostGridPartType, 1, 1 > DiscreteDeformationCoordHolderType;
  typedef typename DiscreteDeformationCoordHolderType :: DiscreteFunctionType CoordinateFunctionType;
  DiscreteDeformationCoordHolderType discreteDeformation( deformation, boundaryProjection, hostGridPart );

  typedef Dune::Fem::GeoGridPart< CoordinateFunctionType > GridPartType;
  GridPartType gridPart( discreteDeformation.coordFunction() );
  //! [Setup the grid part for a deforming domain]

  // type of the mathematical model used
  typedef HeatModel< FunctionSpaceType, GridPartType > ModelType;
  typedef typename ModelType :: ProblemType ProblemType;

  ProblemType* problemPtr = 0;
  switch( problemNumber )
    {
    case 0:
      problemPtr = new SurfaceHeatProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using surface heat problem" << std::endl;
      break;
    case 1:
      problemPtr = new SurfaceMCProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using surface mc problem" << std::endl;
      break;
    case 2:
      problemPtr = new SurfaceMCStationaryProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using stationary surface mc problem" << std::endl;
      break;
    case 3:
      problemPtr = new CurveMCProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using curve mc problem" << std::endl;
      break;
    case 4:
      problemPtr = new CurveMCLessSmoothProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using less smooth curve problem" << std::endl;
      break;
    case 5:
      problemPtr = new CurveMCMoreNonlinearProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using more nonlinear problem" << std::endl;
      break;
    case 6:
      problemPtr = new CurveMCNonUniformProblem< FunctionSpaceType > ( timeProvider );
      std::cout << "using non-uniform curve problem" << std::endl;
      break;

    default:
      std::cerr << "unrecognised problem name" << std::endl;
      assert(0);
    }

  // recover problem
  assert( problemPtr );
  ProblemType& problem = *problemPtr;

  // implicit model for left hand side
  ModelType implicitModel( problem, gridPart, true );

  // explicit model for right hand side
  ModelType explicitModel( problem, gridPart, false );

  // decide how many runs I should do
  const int globalM = Dune::Fem::Parameter::getValue< double >( "mcesfem.M", 1 );
  assert( globalM > 0 );
  const int M = globalM / Dune::Fem::MPIManager::size()
    + ( Dune::Fem::MPIManager::rank() < ( globalM % Dune::Fem::MPIManager::size() ) ? 1 : 0 );

  for( int p = 0; p < Dune::Fem::MPIManager::size(); ++p )
    {
      if( p == Dune::Fem::MPIManager::rank() )
	std::cout << "[" << p << "] " << "my M: " << M << std::endl;
      Dune::Fem::MPIManager::comm().barrier();
    }

  // create heat schemes
  typedef HeatScheme< ModelType, ModelType > SchemeType;
  SchemeType scheme( gridPart, implicitModel, explicitModel, step );
  std::vector< SchemeType > schemeVector( M, scheme );
  typedef MeanScheme< SchemeType > MeanSchemeType;
  MeanSchemeType meanScheme( gridPart, timeProvider, step );

  // initialise random components of scheme
  const unsigned int nYs = Dune::Fem::Parameter::getValue< double >( "mcesfem.number_of_Ys");
  const int seed = Dune::Fem::Parameter::getValue< double >( "mcesfem.rng.seed" );
  const double center = Dune::Fem::Parameter::getValue< double >( "mcesfem.rng.center" );
  const double range = Dune::Fem::Parameter::getValue< double >( "mcesfem.rng.range" );
  static std::mt19937 mt( seed + Dune::Fem::MPIManager::rank() );
  std::uniform_real_distribution<> dist( center-range/2.0, center+range/2.0 );
  for( auto& scheme : schemeVector )
    {
      const std::vector< double > Ys( nYs, dist(mt) );
      scheme.setYs( Ys );

      for( int p = 0; p < Dune::Fem::MPIManager::size(); ++p )
	{
	  if( p == Dune::Fem::MPIManager::rank() )
	    {
	      std::cout << "[" << p << "] Ys:";
	      for( auto Y : Ys )
		{
		  std::cout << " " << Y;
		}
	      std::cout << std::endl;
	    }
	  Dune::Fem::MPIManager::comm().barrier();
	}
    }
  problem.setYs( std::vector<double>( nYs, center ) );

  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename MeanSchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(meanScheme.solution()), &gridExactSolution) ; // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  const double endTime  = Dune::Fem::Parameter::getValue< double >( "heat.endtime", 2.0 );
  const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("heat.reducetimestepfactor", 1 );
  double timeStep = Dune::Fem::Parameter::getValue< double >( "heat.timestep", 0.125 );

  timeStep *= pow(dtreducefactor,step);

  //! [time loop]
  // initialize with fixed time step
  timeProvider.init( timeStep ) ;

  // initialize scheme and output initial data
  for( auto& scheme : schemeVector )
    scheme.initialize();
  meanScheme.computeMean( schemeVector );

  // write initial solve
  dataOutput.write( timeProvider );
  // finalise (compute errors)
  for( auto& scheme : schemeVector )
    scheme.closeTimestep( gridExactSolution, timeProvider.deltaT() );
  meanScheme.closeTimestep( gridExactSolution, timeProvider.deltaT() );

  // increment time
  timeProvider.next( timeStep );

  // time loop, increment with fixed time step
  for( ; timeProvider.time() <= endTime; timeProvider.next( timeStep ) )
  //! [time loop]
  {
    // assemble explicit pare
    for( auto& scheme : schemeVector )
      scheme.prepare();
    //! [Set the new time to move to new surface]
    deformation.setTime( timeProvider.time() + timeProvider.deltaT() );
    discreteDeformation.interpolate();
    // solve once - but now we need to reassmble
    for( auto& scheme : schemeVector )
      scheme.solve(true);
    meanScheme.computeMean( schemeVector );

    //! [Set the new time to move to new surface]
    dataOutput.write( timeProvider );
    // finalise (compute errors)
    for( auto& scheme : schemeVector )
      scheme.closeTimestep( gridExactSolution, timeProvider.deltaT() );
    meanScheme.closeTimestep( gridExactSolution, timeProvider.deltaT() );
  }

  // output final solution
  dataOutput.write( timeProvider );

  // get errors
  std::vector< double > store;
  store.push_back( meanScheme.linftyl2Error() );
  store.push_back( meanScheme.l2h1Error() );
  Dune :: Fem :: FemEoc :: setErrors( eocId, store );

  // write to file / output
  const double h = EvolvingDomain :: GridWidth :: gridWidth( gridPart );
  const int dofs = meanScheme.dofs();
  Dune::Fem::FemEoc::write( h, dofs, 0.0, 0.0, std::cout );
}

// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );

  // append default parameter file
  Dune::Fem::Parameter::append( "../data/parameter" );

  // initialzie eoc file
  std::string eocOutPath = Dune::Fem::Parameter::getValue<std::string>("fem.eocOutPath", std::string("."));
  Dune::Fem::FemEoc::initialize( eocOutPath, "eoc", "surface only" );

  // add entries to eoc calculation
  std::vector<std::string> femEocHeaders;
  femEocHeaders.push_back("$L^\\infty( L^2 )$ error");
  femEocHeaders.push_back("$L^2( H^1 )$ error");

  // get eoc id
  const int eocId = Dune::Fem::FemEoc::addEntry( femEocHeaders );

  // type of hierarchical grid
  typedef Dune :: AlbertaGrid< GRIDDIM, WORLDDIM > HGridType;
  static_assert( HGridType :: dimension == HGridType :: dimensionworld - 1, "this code is written with the assumption grid dim = world dim -1" );

  // create grid from DGF file
  const std::string gridkey = Dune::Fem::IOInterface::defaultGridKey( HGridType::dimension );
  const std::string gridfile = Dune::Fem::Parameter::getValue< std::string >( gridkey );

  // the method rank and size from MPIManager are static
  if( Dune::Fem::MPIManager::rank() == 0 )
    std::cout << "Loading macro grid: " << gridfile << std::endl;

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  // do initial load balance
  grid.loadBalance();

  // setup EOC loop
  const int repeats = Dune::Fem::Parameter::getValue< int >( "heat.repeats", 0 );

  // initial grid refinement
  const int level = Dune::Fem::Parameter::getValue< int >( "heat.level" );

  // number of global refinements to bisect grid width
  const int refineStepsForHalf = Dune::DGFGridInfo< HGridType >::refineStepsForHalf();

  // refine grid
  Dune::Fem::GlobalRefine::apply( grid, level * refineStepsForHalf );

  // calculate first step
  algorithm( grid, (repeats > 0) ? 0 : -1, eocId );

  for( int step = 1; step <= repeats; ++step )
  {
    // refine globally such that grid with is bisected
    // and all memory is adjusted correctly
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

    algorithm( grid, step, eocId );
  }

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}

