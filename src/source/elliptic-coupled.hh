#ifndef ELLIPTIC_COUPLED_HH
#define ELLIPTIC_COUPLED_HH

#include <dune/common/fmatrix.hh>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/grid/common/entitypointer.hh>

// EllipticOperator
// ----------------

//! [Class for mixing operator]
template< class DomainFunction, class RangeFunction, class Model, class CoupledGrid,
	  int dimDiff = DomainFunction :: GridType :: dimension
	  - RangeFunction :: GridType :: dimension >
  struct MixingOperator;

template< class DomainFunction, class RangeFunction, class Model, class CoupledGrid >
struct MixingOperator< DomainFunction, RangeFunction, Model, CoupledGrid, 1 >
  : public virtual Dune :: Fem :: Operator< DomainFunction, RangeFunction >
{
  using BulkDiscreteFunctionType = DomainFunction;
  using SurfaceDiscreteFunctionType = RangeFunction;
  using ModelType = Model;
  using CoupledGridType = CoupledGrid;

protected:

  // typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< typename CoupledGridType ::  SurfaceGeoGridPartType, 0 > QuadratureType;
  // typedef typename QuadratureType::CoordinateType LocalCoordType;

public:
  MixingOperator( const ModelType &model, const CoupledGridType &coupledGrid )
    : model_( model ), coupledGrid_( coupledGrid )
  {}

  //! application operator
  virtual void operator() ( const BulkDiscreteFunctionType &u,
			    SurfaceDiscreteFunctionType &w ) const;

protected:
  const ModelType& model() const { return model_; }
  const CoupledGridType& coupledGrid() const { return coupledGrid_; }

private:
  const ModelType &model_;
  const CoupledGridType &coupledGrid_;
};

template< class DomainFunction, class RangeFunction, class Model, class CoupledGrid >
struct MixingOperator< DomainFunction, RangeFunction, Model, CoupledGrid, -1 >
  : public virtual Dune :: Fem :: Operator< DomainFunction, RangeFunction >
{
  using SurfaceDiscreteFunctionType = DomainFunction;
  using BulkDiscreteFunctionType = RangeFunction;
  using ModelType = Model;
  using CoupledGridType = CoupledGrid;

protected:

  // typedef Dune::Fem::ElementQuadrature< GridPartType, 1 > FaceQuadratureType;
  typedef Dune::Fem::CachingQuadrature< typename CoupledGridType ::  SurfaceGeoGridPartType, 0 > QuadratureType;
  // typedef typename QuadratureType::CoordinateType LocalCoordType;

public:
  MixingOperator( const ModelType &model, const CoupledGridType &coupledGrid )
    : model_( model ), coupledGrid_( coupledGrid )
  {}

  //! application operator
  virtual void operator() ( const SurfaceDiscreteFunctionType &u,
			    BulkDiscreteFunctionType &w ) const;

protected:
  const ModelType& model() const { return model_; }
  const CoupledGridType& coupledGrid() const { return coupledGrid_; }

private:
  const ModelType &model_;
  const CoupledGridType &coupledGrid_;
};

// Implementation of EllipticOperator
// ----------------------------------

template< class DomainFunction, class RangeFunction,
	  class Model, class CoupledGrid >
void MixingOperator< DomainFunction, RangeFunction, Model, CoupledGrid, 1 >
::operator() ( const BulkDiscreteFunctionType &uBulk,
	       SurfaceDiscreteFunctionType &wSurf ) const
{
  // clear destination
  wSurf.clear();

  // extract grid parts
  const auto& surfaceGridPart = coupledGrid().surfaceGridPart();
  const auto& bulkGridPart = coupledGrid().bulkGridPart();

  // iterate over grid
  using SurfaceGridView = typename CoupledGrid :: SurfaceGeoGridPartType :: GridViewType;
  for( const auto& entity : elements( static_cast< SurfaceGridView >( surfaceGridPart ) ) )
  {
    // extract geometry
    const auto& geometry = entity.geometry();

    // find bulk entity
    const auto seed = coupledGrid().surfaceBulkMap( entity );
    const auto bulkEntity = bulkGridPart.entity( seed );
    const auto bulkGeometry = bulkEntity.geometry();

    // get local representation of the discrete functions
    const auto uBulkLocal = uBulk.localFunction( bulkEntity );
    auto wSurfLocal = wSurf.localFunction( entity );

    // construct quadrature
    const int quadOrder = uBulkLocal.order() + wSurfLocal.order();
    QuadratureType quadrature( entity, quadOrder );
    unsigned int numQuadraturePoints = quadrature.nop();

    // quadrature loop
    for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
	// extract quadrature rule
	const auto xLocal = quadrature.point( pt );
	const double weight = quadrature.weight( pt ) * geometry.integrationElement( xLocal );

	// convert quadrature point to local coordinate on bulk element
	const auto xGlobal = geometry.global( xLocal );
	const auto xBulkLocal = bulkGeometry.local( xGlobal );

	// evaluate local functions
	typename SurfaceDiscreteFunctionType :: RangeType cSurfx;
	typename BulkDiscreteFunctionType :: RangeType uBulkx;
	uBulkLocal.evaluate( xBulkLocal, uBulkx );

	// evaluate flux
	model().boundaryFlux( bulkEntity, xBulkLocal, uBulkx, cSurfx );

	// multiply by quadrature weight
	cSurfx *= weight;

	// add to rhs
	wSurfLocal.axpy( quadrature[ pt ], cSurfx );
      }
  }

  // communicate in parallel runs
  wSurf.communicate();
}

template< class DomainFunction, class RangeFunction,
	  class Model, class CoupledGrid >
void MixingOperator< DomainFunction, RangeFunction, Model, CoupledGrid, -1 >
::operator() ( const SurfaceDiscreteFunctionType &uSurf,
	       BulkDiscreteFunctionType &wBulk ) const
{
  // clear destination
  wBulk.clear();

  // extract grid parts
  const auto& surfaceGridPart = coupledGrid().surfaceGridPart();
  const auto& bulkGridPart = coupledGrid().bulkGridPart();

  // iterate over grid
  using SurfaceGridView = typename CoupledGrid :: SurfaceGeoGridPartType :: GridViewType;
  for( const auto& entity : elements( static_cast< SurfaceGridView >( surfaceGridPart ) ) )
  {
    // exact geometry
    const auto& geometry = entity.geometry();

    // find bulk entity
    const auto seed = coupledGrid().surfaceBulkMap( entity );
    const auto bulkEntity = bulkGridPart.entity( seed );
    const auto bulkGeometry = bulkEntity.geometry();

    // get local representation of the discrete functions
    const auto uSurfLocal = uSurf.localFunction( entity );
    auto wBulkLocal = wBulk.localFunction( bulkEntity );

    // construct quadrature
    const int quadOrder = uSurfLocal.order() + wBulkLocal.order();
    QuadratureType quadrature( entity, quadOrder );
    unsigned int numQuadraturePoints = quadrature.nop();

    // quadrature loop
    for( unsigned int pt = 0; pt < numQuadraturePoints; ++pt )
      {
	// extract quadrature rule
	const auto xLocal = quadrature.point( pt );
	const double weight = quadrature.weight( pt ) * geometry.integrationElement( xLocal );

	// convert quadrature point to local coordinate on bulk element
	const auto xGlobal = geometry.global( xLocal );
	const auto xBulkLocal = bulkGeometry.local( xGlobal );

	// evaluate local functions
	typename SurfaceDiscreteFunctionType :: RangeType uSurfx;
	uSurfLocal.evaluate( quadrature[ pt ], uSurfx );
	typename BulkDiscreteFunctionType :: RangeType cBulkx;

	// evaluate flux
	model().boundaryFlux( entity, xLocal, uSurfx, cBulkx );

	// multiply by quadrature weight
	cBulkx *= weight;

	// add to rhs
	wBulkLocal.axpy( xBulkLocal, cBulkx );
      }
  }

  // communicate in parallel runs
  wBulk.communicate();
}


// DifferentiableEllipticOperator
// ------------------------------

//! [Class for linearizable mixing operator]
template< class JacobianOperator, class Model, class CoupledGrid,
	  int dimDiff = JacobianOperator :: DomainFunctionType :: GridType :: dimension
	  - JacobianOperator :: RangeFunctionType :: GridType :: dimension >
  struct DifferentiableMixingOperator;

template< class JacobianOperator, class Model, class CoupledGrid >
struct DifferentiableMixingOperator< JacobianOperator, Model, CoupledGrid, 1 >
  : public MixingOperator< typename JacobianOperator :: DomainFunctionType,
			   typename JacobianOperator :: RangeFunctionType,
			   Model, CoupledGrid, 1 >,
    public Dune :: Fem :: DifferentiableOperator< JacobianOperator >
{
  using BaseType = MixingOperator< typename JacobianOperator :: DomainFunctionType,
				   typename JacobianOperator :: RangeFunctionType,
				   Model, CoupledGrid >;

  using JacobianOperatorType = JacobianOperator;

  using ModelType = typename BaseType :: ModelType;
  using CoupledGridType = typename BaseType :: CoupledGridType;

  using BulkDiscreteFunctionType = typename BaseType :: BulkDiscreteFunctionType;
  using SurfaceDiscreteFunctionType = typename BaseType :: SurfaceDiscreteFunctionType;

  using QuadratureType = typename BaseType :: QuadratureType;

public:
  //! constructor
  DifferentiableMixingOperator( const ModelType& model, const CoupledGridType& coupledGrid )
    : BaseType( model, coupledGrid )
  {}

  //! method to setup the jacobian of the operator for storage in a matrix
  void jacobian ( const BulkDiscreteFunctionType &u, JacobianOperatorType &jOp ) const;

protected:
  using BaseType::model;
  using BaseType::coupledGrid;
};

template< class JacobianOperator, class Model, class CoupledGrid >
struct DifferentiableMixingOperator< JacobianOperator, Model, CoupledGrid, -1 >
  : public MixingOperator< typename JacobianOperator :: DomainFunctionType,
			   typename JacobianOperator :: RangeFunctionType,
			   Model, CoupledGrid, -1 >,
    public Dune :: Fem :: DifferentiableOperator< JacobianOperator >
{
  using BaseType = MixingOperator< typename JacobianOperator :: DomainFunctionType,
				   typename JacobianOperator :: RangeFunctionType,
				   Model, CoupledGrid >;

  using JacobianOperatorType = JacobianOperator;

  using ModelType = typename BaseType :: ModelType;
  using CoupledGridType = typename BaseType :: CoupledGridType;

  using SurfaceDiscreteFunctionType = typename BaseType :: SurfaceDiscreteFunctionType;
  using BulkDiscreteFunctionType = typename BaseType :: BulkDiscreteFunctionType;

  using QuadratureType = typename BaseType :: QuadratureType;

public:
  //! constructor
  DifferentiableMixingOperator( const ModelType& model, const CoupledGridType& coupledGrid )
    : BaseType( model, coupledGrid )
  {}

  //! method to setup the jacobian of the operator for storage in a matrix
  void jacobian ( const SurfaceDiscreteFunctionType &u, JacobianOperatorType &jOp ) const;

protected:
  using BaseType::model;
  using BaseType::coupledGrid;
};

// Implementation of DifferentiableEllipticOperator
// ------------------------------------------------

template < class JacobianOperator, class Model, class CoupledGrid >
void DifferentiableMixingOperator< JacobianOperator, Model, CoupledGrid, 1 >
::jacobian ( const BulkDiscreteFunctionType &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;

  using BulkDiscreteFunctionSpaceType = typename BulkDiscreteFunctionType :: DiscreteFunctionSpaceType;
  using SurfaceDiscreteFunctionSpaceType = typename SurfaceDiscreteFunctionType :: DiscreteFunctionSpaceType;

  // get discrete function spaces
  const auto& domainSpace = jOp.domainSpace();
  const auto& rangeSpace = jOp.rangeSpace();

  // get grid parts
  const auto& surfaceGridPart = coupledGrid().surfaceGridPart();
  const auto& bulkGridPart = coupledGrid().bulkGridPart();
  // set up matrix stencil
  Dune :: Fem :: Stencil< BulkDiscreteFunctionSpaceType, SurfaceDiscreteFunctionSpaceType > stencil( domainSpace, rangeSpace );
  // iterate over grid
  using SurfaceGridView = typename CoupledGrid :: SurfaceGeoGridPartType :: GridViewType;
  for( const auto& entity : elements( static_cast< SurfaceGridView >( surfaceGridPart ) ) )
  {
    // find bulk entity
    const auto seed = coupledGrid().surfaceBulkMap( entity );
    const auto bulkEntity = bulkGridPart.entity( seed );

    stencil.fill( bulkEntity, entity );
  }

  jOp.reserve( stencil );
  jOp.clear();

  // set up basis function storage
  const int domainBlockSize = domainSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename BulkDiscreteFunctionType::RangeType > domainPhi( domainSpace.blockMapper().maxNumDofs()*domainBlockSize );
  const int rangeBlockSize = rangeSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename SurfaceDiscreteFunctionType::RangeType > rangePhi( rangeSpace.blockMapper().maxNumDofs()*rangeBlockSize );

  // loop over grid
  using SurfaceGridView = typename CoupledGrid :: SurfaceGeoGridPartType :: GridViewType;
  for( const auto& entity : elements( static_cast< SurfaceGridView >( surfaceGridPart ) ) )
    {
      // find geometry
      const auto& geometry = entity.geometry();

      // find bulk entity
      const auto seed = coupledGrid().surfaceBulkMap( entity );
      const auto bulkEntity = bulkGridPart.entity( seed );
      const auto bulkGeometry = bulkEntity.geometry();

      // construct local representations of data
      const auto uLocal = u.localFunction( bulkEntity );
      auto jLocal = jOp.localMatrix( bulkEntity, entity );

      // find local basis functions
      const auto& domainBasisSet = jLocal.domainBasisFunctionSet();
      const unsigned int numDomainBasisFunctions = domainBasisSet.size();
      const auto& rangeBasisSet = jLocal.rangeBasisFunctionSet();

      // perform quadrature loop
      QuadratureType quadrature( entity, domainSpace.order() + rangeSpace.order() );
      size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  // obatain quadrature points in quadrature coords
	  const auto xLocal = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) *
	    geometry.integrationElement( xLocal );

	  const auto xGlobal = geometry.global( xLocal );
	  const auto xBulkLocal = bulkGeometry.local( xGlobal );

	  // evaluate basis functions at quadrature points
	  domainBasisSet.evaluateAll( xBulkLocal, domainPhi );
	  rangeBasisSet.evaluateAll( quadrature[ pt ], rangePhi );

	  // get value for linearization
	  typename BulkDiscreteFunctionType :: RangeType ux;
	  uLocal.evaluate( xBulkLocal, ux );

	  for( unsigned int i = 0; i < numDomainBasisFunctions; ++i )
	    {
	      // evaluate fluxes
	      typename SurfaceDiscreteFunctionType :: RangeType cphi;
	      model().linBoundaryFlux( ux, bulkEntity, xBulkLocal, domainPhi[ i ], cphi );

	      // add contribution to local matrix
	      jLocal.column( i ).axpy( rangePhi, cphi, weight );
	    }
	}
    }
  jOp.communicate();
}

template < class JacobianOperator, class Model, class CoupledGrid >
void DifferentiableMixingOperator< JacobianOperator, Model, CoupledGrid, -1 >
::jacobian ( const SurfaceDiscreteFunctionType &u, JacobianOperator &jOp ) const
{
  typedef typename JacobianOperator::LocalMatrixType LocalMatrixType;

  using SurfaceDiscreteFunctionSpaceType = typename SurfaceDiscreteFunctionType :: DiscreteFunctionSpaceType;
  using BulkDiscreteFunctionSpaceType = typename BulkDiscreteFunctionType :: DiscreteFunctionSpaceType;

  // get discrete function spaces
  const auto& domainSpace = jOp.domainSpace();
  const auto& rangeSpace = jOp.rangeSpace();

  // get grid parts
  const auto& surfaceGridPart = coupledGrid().surfaceGridPart();
  const auto& bulkGridPart = coupledGrid().bulkGridPart();
  // set up matrix stencil
  Dune :: Fem :: Stencil< SurfaceDiscreteFunctionSpaceType, BulkDiscreteFunctionSpaceType > stencil( domainSpace, rangeSpace );
  // iterate over grid
  using SurfaceGridView = typename CoupledGrid :: SurfaceGeoGridPartType :: GridViewType;
  for( const auto& entity : elements( static_cast< SurfaceGridView >( surfaceGridPart ) ) )
  {
    // find bulk entity
    const auto seed = coupledGrid().surfaceBulkMap( entity );
    const auto bulkEntity = bulkGridPart.entity( seed );

    stencil.fill( entity, bulkEntity );
  }

  jOp.reserve( stencil );
  jOp.clear();

  // set up basis function storage
  const int domainBlockSize = domainSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename BulkDiscreteFunctionType::RangeType > domainPhi( domainSpace.blockMapper().maxNumDofs()*domainBlockSize );
  const int rangeBlockSize = rangeSpace.localBlockSize; // is equal to 1 for scalar functions
  std::vector< typename SurfaceDiscreteFunctionType::RangeType > rangePhi( rangeSpace.blockMapper().maxNumDofs()*rangeBlockSize );

  // loop over grid
  using SurfaceGridView = typename CoupledGrid :: SurfaceGeoGridPartType :: GridViewType;
  for( const auto& entity : elements( static_cast< SurfaceGridView >( surfaceGridPart ) ) )
    {
      // find geometry
      const auto& geometry = entity.geometry();

      // find bulk entity
      const auto seed = coupledGrid().surfaceBulkMap( entity );
      const auto bulkEntity = bulkGridPart.entity( seed );
      const auto bulkGeometry = bulkEntity.geometry();

      // construct local representations of data
      const auto uLocal = u.localFunction( entity );
      auto jLocal = jOp.localMatrix( entity, bulkEntity );

      // find local basis functions
      const auto& domainBasisSet = jLocal.domainBasisFunctionSet();
      const unsigned int numDomainBasisFunctions = domainBasisSet.size();
      const auto& rangeBasisSet = jLocal.rangeBasisFunctionSet();

      // perform quadrature loop
      QuadratureType quadrature( entity, domainSpace.order() + rangeSpace.order() );
      size_t numQuadraturePoints = quadrature.nop();
      for( size_t pt = 0; pt < numQuadraturePoints; ++pt )
	{
	  // obatain quadrature points in quadrature coords
	  const auto xLocal = quadrature.point( pt );
	  const double weight = quadrature.weight( pt ) *
	    geometry.integrationElement( xLocal );

	  const auto xGlobal = geometry.global( xLocal );
	  const auto xBulkLocal = bulkGeometry.local( xGlobal );

	  // evaluate basis functions at quadrature points
	  domainBasisSet.evaluateAll( quadrature[ pt ], domainPhi );
	  rangeBasisSet.evaluateAll( xBulkLocal, rangePhi );

	  // get value for linearization
	  typename SurfaceDiscreteFunctionType :: RangeType ux;
	  uLocal.evaluate( quadrature[ pt ], ux );

	  for( unsigned int i = 0; i < numDomainBasisFunctions; ++i )
	    {
	      // evaluate fluxes
	      typename BulkDiscreteFunctionType :: RangeType cphi;
	      model().linBoundaryFlux( ux, entity, xLocal, domainPhi[ i ], cphi );

	      // add contribution to local matrix
	      jLocal.column( i ).axpy( rangePhi, cphi, weight );
	    }
	}
    }
  jOp.communicate();
}

#endif // #ifndef ELLIPTIC_COUPLED_HH
