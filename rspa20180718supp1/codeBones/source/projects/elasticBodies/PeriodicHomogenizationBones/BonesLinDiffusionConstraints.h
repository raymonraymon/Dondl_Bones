#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESDIFFUSIONCONSTRAINTS_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESDIFFUSIONCONSTRAINTS_H


#include "BonesMatOptDefines.h"
#include "BonesOptDiffusionSolver.h"

using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{  

// Compl = E_PER + 2 E_MIXED + E_AFF 
 //   = (for optimal diffusion) E_MIXED + E_AFF
template<typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class DiffusionConstraintMultipleLoad {
  
  typedef typename MatOptConfiguratorType::RealType RealType;
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

protected:
  const OptimalDiffusionSolverMultipleLoad<MatOptConfiguratorType> & _OptDiffusionFinder;
  
public:
  DiffusionConstraintMultipleLoad( const OptimalDiffusionSolverMultipleLoad<MatOptConfiguratorType> & OptDiffusionFinder  ) : 
  _OptDiffusionFinder ( OptDiffusionFinder ) { }

  void apply( const VectorType & v, const int numLoad, RealType &dest ) const { 
     _OptDiffusionFinder.updatePhasefield( v );

       dest = 0.0;
       RealType _compliance_Mixed = -1. * ( _OptDiffusionFinder.template getRHS<MaterialType>(numLoad) ).dot( _OptDiffusionFinder.template getSolutionDiffusionAndMultiplier<MaterialType>(numLoad) );
       RealType _compliance_Affine = ( _OptDiffusionFinder.template getHessianLinDiffusionAffinePart<MaterialType>() * _OptDiffusionFinder.template getAffineDiffusion<MaterialType>(numLoad) ).dot( _OptDiffusionFinder.template getAffineDiffusion<MaterialType>(numLoad));
       dest = _compliance_Mixed + _compliance_Affine;
     

  }
  
};


//======================================================================================================================================
//================================= Derivative of Jphys  ===============================================================================
//======================================================================================================================================


template<typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType >
class DiffusionConstraintMultipleLoad_DerivativeInM
: public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, DiffusionConstraintMultipleLoad_DerivativeInM<MatOptConfiguratorType,MaterialType> > {
protected: 

    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename MatOptConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DomVecType DomVecType;
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;

    const OptimalDiffusionSolverMultipleLoad<MatOptConfiguratorType> &_OptDiffusionFinder;
    const Material<RealType> &_MaterialBone, &_MaterialPolymer;
    RealType _diffusionCoefficient;
    const RealType _factorVoidMaterial;
    const int _dimDomain;
    mutable DomVecType _Grad_diffusionAffine;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_diffusionPeriodicPtr;

public:

    DiffusionConstraintMultipleLoad_DerivativeInM ( const OptimalDiffusionSolverMultipleLoad<MatOptConfiguratorType> & OptDiffusionFinder ) 
    : QuocFENonlinOpIntegrator<ConfiguratorType, DiffusionConstraintMultipleLoad_DerivativeInM<MatOptConfiguratorType,MaterialType> > (
    OptDiffusionFinder.getMatOptConfigurator()._conf ),
    _OptDiffusionFinder ( OptDiffusionFinder ),
    _MaterialBone ( OptDiffusionFinder.getMatOptConfigurator()._materialInfo._MaterialBone ),
    _MaterialPolymer ( OptDiffusionFinder.getMatOptConfigurator()._materialInfo._MaterialPolymer ),
    _factorVoidMaterial( OptDiffusionFinder.getMatOptConfigurator()._factorVoidMaterial ),
    _dimDomain( OptDiffusionFinder.getMatOptConfigurator()._conf.dimDomain ),
    _PfPtr ( NULL ), _diffusionPeriodicPtr ( NULL ) {
      switch( MaterialType ){
        case BONE :
            _diffusionCoefficient = _MaterialBone.getDiffusionCoefficient(); 
            break;
        case POLYMER :
            _diffusionCoefficient = _MaterialPolymer.getDiffusionCoefficient();
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      } 
    }

    ~DiffusionConstraintMultipleLoad_DerivativeInM() {
        delete _PfPtr;
        delete _diffusionPeriodicPtr;
    };

    void apply ( const VectorType &Pf, const int numLoad, VectorType &Dest ) const {
        Dest.setZero();
        _OptDiffusionFinder.updatePhasefield( Pf ); // This computes opt diffusion u(Pf) and by reference this updated in adjoint problem
             
        delete _PfPtr;
        _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _OptDiffusionFinder.getMatOptConfigurator()._conf, _OptDiffusionFinder.getPhaseFieldPeriodicallyExtended() );
        

        delete _diffusionPeriodicPtr;
        _diffusionPeriodicPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _OptDiffusionFinder.getMatOptConfigurator()._conf, _OptDiffusionFinder.template getSolutionDiffusionPeriodicallyExtended<MaterialType>(numLoad) );

//          QuocDiscreteFunctionDefaultAffine<DataTypeContainer, ConfiguratorType::dimDomain> diffusionAffine ( _OptDiffusionFinder.template getAffineDiffusion<MaterialType>(i) );
            _Grad_diffusionAffine = _OptDiffusionFinder.template getAffineDiffusion<MaterialType> ( numLoad );
            
        QuocFENonlinOpIntegrator< ConfiguratorType, DiffusionConstraintMultipleLoad_DerivativeInM<MatOptConfiguratorType,MaterialType> >::assembleAdd( Dest );
        
    }

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {

        DomVecType Grad_diffusion;
        _diffusionPeriodicPtr->evaluateGradientAtQuadPoint( El, QuadPoint, Grad_diffusion );

        const RealType Dchi = _OptDiffusionFinder.getMatOptConfigurator().template approxCharFct_material_Derivative<MaterialType>( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
        const RealType materialfactor = Dchi * ( 1. - _factorVoidMaterial );
        
        RealType aux = 0.0;
        RealType auxPeriodic = 0.0, auxMixed = 0.0, auxAffine = 0.0;

        auxPeriodic =  materialfactor *  _diffusionCoefficient * Grad_diffusion.dot( Grad_diffusion );
        auxMixed =   materialfactor * _diffusionCoefficient * Grad_diffusion.dot( _Grad_diffusionAffine );
        auxAffine =  materialfactor * _Grad_diffusionAffine.dot( _Grad_diffusionAffine );

        aux = auxPeriodic + 2. * auxMixed + auxAffine;
        
        return aux; 
    }

};




template <typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class DiffusionOpInMaterialPeriodicBC : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
protected:
  const MatOptConfigurator &_matOptConf;
  const QuocHandler<ConfiguratorType> & _quocHandler;
  const OptimalDiffusionSolverMultipleLoad<MatOptConfigurator> & _OptDiffusionFinder;
  const int _numConstraint;

public:
  DiffusionOpInMaterialPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler, 
                                               const OptimalDiffusionSolverMultipleLoad<MatOptConfigurator> & OptDiffusionFinder,
                                               const int numConstraint   )  :
//   aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer > ( ),
  _matOptConf ( matOptConf ), _quocHandler ( quocHandler ), 
  _OptDiffusionFinder ( OptDiffusionFinder ), _numConstraint( numConstraint ) { } 
      
    void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
        
        DiffusionConstraintMultipleLoad<MatOptConfigurator,MaterialType>( _OptDiffusionFinder ).apply( vCollabsedAndExtended, _numConstraint, energy );
    }
    
    void evaluateJacobian( const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
        
        DiffusionConstraintMultipleLoad_DerivativeInM<MatOptConfigurator,MaterialType>( _OptDiffusionFinder ).apply( vCollabsedAndExtended, _numConstraint, Deriv );
        _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
        throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
      
};



template <typename MatOptConfigurator> //, MaterialTypeBonePolymer MaterialType>
class DiffusionAndBarycenterConstraintPeriodicBC : public aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
protected:
  const MatOptConfigurator &_matOptConf;
  const QuocHandler<ConfiguratorType> & _quocHandler;
  const OptimalDiffusionSolverMultipleLoad<MatOptConfigurator> & _OptDiffusionFinder;
  const PointType &_c;

  public:
  DiffusionAndBarycenterConstraintPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler, 
                                                          const OptimalDiffusionSolverMultipleLoad<MatOptConfigurator> & OptDiffusionFinder,
                                                          const PointType &c )  :
  aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  > ( OptDiffusionFinder.template getNumLoads<BONE>() + OptDiffusionFinder.template getNumLoads<POLYMER>() + c.size() ),
  _matOptConf ( matOptConf ), _quocHandler ( quocHandler ), 
  _OptDiffusionFinder ( OptDiffusionFinder ),
  _c ( c ) { } 
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
        
        if( numConstraint < _OptDiffusionFinder.template getNumLoads<BONE>() ){
                const int numConstraintBone = numConstraint;
                DiffusionConstraintMultipleLoad<MatOptConfigurator,BONE>( _OptDiffusionFinder ).apply( vCollabsedAndExtended, numConstraintBone, energy );
        }else{ 
            if( numConstraint < _OptDiffusionFinder.template getNumLoads<BONE>() + _OptDiffusionFinder.template getNumLoads<POLYMER>() ){
                 const int numConstraintPolymer = numConstraint - _OptDiffusionFinder.template getNumLoads<BONE>();
                 DiffusionConstraintMultipleLoad<MatOptConfigurator,POLYMER>( _OptDiffusionFinder ).apply( vCollabsedAndExtended, numConstraintPolymer, energy );
            } else {
                const int numConstraintBarycenter = numConstraint - _OptDiffusionFinder.template getNumLoads<BONE>() - _OptDiffusionFinder.template getNumLoads<POLYMER>();
                PfOpBones_Barycenter<MatOptConfigurator,BONE> ( _matOptConf, numConstraintBarycenter, _c[numConstraintBarycenter] ).apply( vCollabsedAndExtended, energy );
            }
        }
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
        
        if( numConstraint < _OptDiffusionFinder.template getNumLoads<BONE>() ){
                const int numConstraintBone = numConstraint;
                DiffusionConstraintMultipleLoad_DerivativeInM<MatOptConfigurator,BONE>( _OptDiffusionFinder ).apply( vCollabsedAndExtended, numConstraintBone, Deriv );
        }else{ 
            if( numConstraint < _OptDiffusionFinder.template getNumLoads<BONE>() + _OptDiffusionFinder.template getNumLoads<POLYMER>() ){
                 const int numConstraintPolymer = numConstraint - _OptDiffusionFinder.template getNumLoads<BONE>();
                 DiffusionConstraintMultipleLoad_DerivativeInM<MatOptConfigurator,POLYMER>( _OptDiffusionFinder ).apply( vCollabsedAndExtended, numConstraintPolymer, Deriv );
            } else {
                const int numConstraintBarycenter = numConstraint - _OptDiffusionFinder.template getNumLoads<BONE>() - _OptDiffusionFinder.template getNumLoads<POLYMER>();
                PfOpBones_BarycenterDerivative<MatOptConfigurator,BONE> ( _matOptConf, numConstraintBarycenter, _c[numConstraintBarycenter] ).apply( vCollabsedAndExtended, Deriv );
            }
        }
        
        _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
//     const PointType& getBarycenterPoint( ) const {return _c;};
      
};




}; //end namespace


#endif
