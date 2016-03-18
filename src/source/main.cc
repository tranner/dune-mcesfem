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
// include description of surface deformation
#include "deformation.hh"

// include discrete function space
#include <dune/fem/space/lagrange.hh>
// include discrete function
#include <dune/fem/function/blockvectorfunction.hh>
// lagrange interpolation 
#include <dune/fem/operator/lagrangeinterpolation.hh>

// include header for heat model
#include "heat.hh"

#include "heatmodel.hh"
#include "heatscheme.hh"

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
  typedef DeformationCoordFunction< HGridType::dimensionworld > DeformationType;
  DeformationType deformation;

  typedef DiscreteDeformationCoordHolder< DeformationType, HostGridPartType, 1 > DiscreteDeformationCoordHolderType;
  typedef typename DiscreteDeformationCoordHolderType :: DiscreteFunctionType CoordinateFunctionType;
  DiscreteDeformationCoordHolderType discreteDeformation( deformation, hostGridPart );

  typedef Dune::Fem::GeoGridPart< CoordinateFunctionType > GridPartType;
  GridPartType gridPart( discreteDeformation.coordFunction() );
  //! [Setup the grid part for a deforming domain]

  // type of the mathematical model used
  typedef HeatModel< FunctionSpaceType, GridPartType > ModelType;
  typedef typename ModelType :: ProblemType ProblemType;

  // choose problem
  ProblemType* problemPtr = 0;
  const std::string problemNames [] = { "surface_heat", "surface_mc" };
  const int problemNumber = Dune :: Fem :: Parameter :: getEnum( "heat.problem", problemNames );
  switch( problemNumber )
    {
    case 0:
      problemPtr = new SurfaceHeatProblem< FunctionSpaceType > ( timeProvider );
      break;
    case 1:
      problemPtr = new SurfaceMCProblem< FunctionSpaceType > ( timeProvider );
      break;

    default:
      std::cerr << "unrecognised problem name" << std::endl;
    }

  // recover problem
  assert( problemPtr );
  ProblemType& problem = *problemPtr;

  // implicit model for left hand side
  ModelType implicitModel( problem, gridPart, true );

  // explicit model for right hand side
  ModelType explicitModel( problem, gridPart, false );

  // create heat schemes
  const int M = Dune::Fem::Parameter::getValue< double >( "mcesfem.M", 1 );
  assert( M > 0 );
  typedef HeatScheme< ModelType, ModelType > SchemeType;
  SchemeType scheme( gridPart, implicitModel, explicitModel, step );
  std::vector< SchemeType > schemeVector( M, scheme );
  SchemeType meanScheme( gridPart, implicitModel, explicitModel, step );

  // initialise random components of scheme
  const int seed = Dune::Fem::Parameter::getValue< double >( "mcesfem.seed" );
  std::mt19937 mt( 1729 );
  std::uniform_real_distribution<> dist(-1, 1);
  for( auto& scheme : schemeVector )
    {
      const double Y1 = dist(mt);
      const double Y2 = dist(mt);
      scheme.setY1Y2( Y1, Y2 );
      std::cout << "Y1: " << Y1 << " Y2: " << Y2 << std::endl;
    }

  typedef Dune::Fem::GridFunctionAdapter< ProblemType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", problem, gridPart, 5 );
  //! input/output tuple and setup datawritter
  typedef Dune::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(meanScheme.solution()), &gridExactSolution) ; // tuple with pointers
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );

  const double endTime  = Dune::Fem::Parameter::getValue< double >( "heat.endtime", 2.0 );
  const double dtreducefactor = Dune::Fem::Parameter::getValue< double >("heat.reducetimestepfactor", 1 );
  double timeStep = Dune::Fem::Parameter::getValue< double >( "heat.timestep", 0.125 );

  // timeStep *= pow(dtreducefactor,step);

  //! [time loop]
  // initialize with fixed time step
  timeProvider.init( timeStep ) ;

  // initialize scheme and output initial data
  for( auto& scheme : schemeVector )
    scheme.initialize();

  // write initial solve
  dataOutput.write( timeProvider );

  // time loop, increment with fixed time step
  for( ; timeProvider.time() < endTime; timeProvider.next( timeStep ) )
  //! [time loop]
  {
    // assemble explicit pare
    for( auto& scheme : schemeVector )
      scheme.prepare();
    //! [Set the new time to move to new surface]
    discreteDeformation.setTime( timeProvider.time() + timeProvider.deltaT() );
    // solve once - but now we need to reassmble
    for( auto& scheme : schemeVector )
      scheme.solve(true);
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
  std::vector< double > store = { meanScheme.linftyl2Error(), meanScheme.l2h1Error() };
  Dune :: Fem :: FemEoc :: setErrors( eocId, store );

  for( auto& scheme : schemeVector )
    std::cout << "L2: " << scheme.linftyl2Error() << " H1: " << scheme.l2h1Error() << std::endl;

  // write to file / output
  const double h = EvolvingDomain :: GridWidth :: gridWidth( gridPart );
  const int dofs = 0;//scheme.dofs();
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
  typedef Dune :: AlbertaGrid< 2, 3 > HGridType;
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

