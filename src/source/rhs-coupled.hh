#ifndef RHSBULK_HH
#define RHSBULK_HH

#include <dune/fem/quadrature/cachingquadrature.hh>

// assembleRHS
// -----------

template< class CoupledModel,
	  class BulkDiscreteFunction, class SurfaceDiscreteFunction,
	  class CoupledGrid >
void assembleRHS ( const CoupledModel &coupledModel,
		   const BulkDiscreteFunction &bulkSolution,
		   const SurfaceDiscreteFunction &surfaceSolution,
		   const CoupledGrid &coupledGrid,
		   BulkDiscreteFunction &bulkRhs, SurfaceDiscreteFunction &surfaceRhs )
{
  typedef typename CoupledModel :: BulkModelType BulkModelType;
  typedef typename CoupledModel :: SurfaceModelType SurfaceModelType;

  typedef typename BulkModelType :: TimeProviderType TimeProviderType;

  typedef typename BulkModelType :: RightHandSideType BulkRightHandSideType;
  typedef typename SurfaceModelType :: RightHandSideType SurfaceRightHandSideType;

  typedef typename BulkDiscreteFunction::DiscreteFunctionSpaceType BulkDiscreteFunctionSpaceType;
  typedef typename BulkDiscreteFunction::LocalFunctionType BulkLocalFunctionType;
  typedef typename SurfaceDiscreteFunction::DiscreteFunctionSpaceType SurfaceDiscreteFunctionSpaceType;
  typedef typename SurfaceDiscreteFunction::LocalFunctionType SurfaceLocalFunctionType;

  typedef typename BulkDiscreteFunctionSpaceType::IteratorType BulkIteratorType;
  typedef typename BulkIteratorType::Entity BulkEntityType;
  typedef typename BulkEntityType::Geometry BulkGeometryType;

  typedef typename SurfaceDiscreteFunctionSpaceType::IteratorType SurfaceIteratorType;
  typedef typename SurfaceIteratorType::Entity SurfaceEntityType;
  typedef typename SurfaceEntityType::Geometry SurfaceGeometryType;

  typedef typename BulkDiscreteFunctionSpaceType::GridPartType BulkGridPartType;
  typedef typename SurfaceDiscreteFunctionSpaceType::GridPartType SurfaceGridPartType;

  typedef typename BulkDiscreteFunction::DomainType GlobalCoordType;
  typedef typename BulkLocalFunctionType::DomainType LocalCoordType;
  typedef typename SurfaceLocalFunctionType::DomainType SurfaceLocalCoordType;

  typedef Dune::Fem::CachingQuadrature< BulkGridPartType, 0 > BulkQuadratureType;
  typedef Dune::Fem::CachingQuadrature< SurfaceGridPartType, 0 > SurfaceQuadratureType;

  const BulkModelType &bulkModel = coupledModel.bulkModel();
  const SurfaceModelType &surfaceModel = coupledModel.surfaceModel();
  const TimeProviderType &timeProvider = bulkModel.timeProvider();

  const BulkDiscreteFunctionSpaceType &bulkDfSpace = bulkRhs.space();

  const BulkIteratorType bend = bulkDfSpace.end();
  for( BulkIteratorType bit = bulkDfSpace.begin(); bit != bend; ++bit )
    {
      const BulkEntityType &entity = *bit;
      const BulkGeometryType &geometry = entity.geometry();

      // bulk forcing term f_h
      {
	const BulkLocalFunctionType bulkLocal = bulkSolution.localFunction( entity );
	BulkLocalFunctionType bulkRhsLocal = bulkRhs.localFunction( entity );

	BulkQuadratureType quadrature( entity, 2*bulkDfSpace.order() + 1 );
	size_t numQuadraturePoints = quadrature.nop();

	for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	  {
	    // obatain quadrature points in quadrature coords
	    const LocalCoordType &xLocal = quadrature.point( pt );
	    const double weight = quadrature.weight( pt ) *
	      geometry.integrationElement( xLocal );

	    // map to global coord
	    GlobalCoordType xGlobal = geometry.global( xLocal );

	    // evaluate f
	    typename BulkRightHandSideType :: RangeType fx;
	    bulkModel.rightHandSide().evaluate( xGlobal, fx );

	    // multiply by quadrature weight (and deltaT)
	    fx *= weight * timeProvider.deltaT();

	    // add fx * phi_i to rhsLocal[ i ]
	    bulkRhsLocal.axpy( xLocal, fx );
	  }
      }
    }

  const SurfaceDiscreteFunctionSpaceType &surfaceDfSpace = surfaceRhs.space();

  const SurfaceIteratorType send = surfaceDfSpace.end();
  for( SurfaceIteratorType sit = surfaceDfSpace.begin(); sit != send; ++sit )
    {
      const SurfaceEntityType &entity = *sit;
      const SurfaceGeometryType &geometry = entity.geometry();

      const typename BulkEntityType :: EntitySeed bulkSeed
	= coupledGrid.surfaceBulkMap( entity );
      const BulkEntityType &bulkEntity = bulkDfSpace.gridPart().entity( bulkSeed );
      const BulkGeometryType &bulkGeometry = bulkEntity.geometry();

      // surface forcing term g_h
      {
	const SurfaceLocalFunctionType surfaceLocal = surfaceSolution.localFunction( entity );
	const BulkLocalFunctionType bulkLocal = bulkSolution.localFunction( bulkEntity );
	SurfaceLocalFunctionType surfaceRhsLocal = surfaceRhs.localFunction( entity );
	BulkLocalFunctionType bulkRhsLocal = bulkRhs.localFunction( bulkEntity );

	SurfaceQuadratureType quadrature( entity, 2*surfaceDfSpace.order() + 1 );
	size_t numQuadraturePoints = quadrature.nop();

	for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	  {
	    // obatain quadrature points in quadrature coords
	    const typename SurfaceQuadratureType :: CoordinateType &xLocal
	      = quadrature.point( pt );
	    const double weight = quadrature.weight( pt ) *
	      geometry.integrationElement( xLocal );

	    // map to global coord
	    GlobalCoordType xGlobal = geometry.global( xLocal );

	    // evaluate f
	    typename SurfaceRightHandSideType :: RangeType fx;
	    surfaceModel.rightHandSide().evaluate( xGlobal, fx );

	    // multiply by quadrature weight (and deltaT)
	    fx *= weight * timeProvider.deltaT();

	    // add fx * phi_i to rhsLocal[ i ]
	    surfaceRhsLocal.axpy( xLocal, fx );

	    // map global coord to local coord in bulk entity
	    typename BulkGeometryType :: LocalCoordinate xBulkLocal
	      = bulkGeometry.local( xGlobal );

	    // add bulk to surface coupling
	    typename BulkLocalFunctionType :: RangeType uhx;
	    bulkLocal.evaluate( xBulkLocal, uhx );
	    uhx *= weight * coupledModel.alpha();
	    surfaceRhsLocal.axpy( xLocal, uhx );

	    // add surface to bulk coupling
	    typename SurfaceLocalFunctionType :: RangeType vhx;
	    surfaceLocal.evaluate( xLocal, vhx );
	    vhx *= weight * coupledModel.beta();
	    bulkRhsLocal.axpy( xBulkLocal, vhx );
	  }
      }
    }

  bulkRhs.communicate();
  surfaceRhs.communicate();
}

#endif // #ifndef RHS_HH
