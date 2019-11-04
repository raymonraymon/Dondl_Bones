#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESCOMPLIANCEENERGIES_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADMATERIALOPTIMIZATIONBONESCOMPLIANCEENERGIES_H


#include "BonesMatOptDefines.h"
#include "BonesOptDeformSolver.h"


using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{    
    
    
//======================================================================================================================================
//================================= Interfaces  ===============================================================================
//======================================================================================================================================
    
    
//======================================================================================================================================
//================================= compliance = max_alpha( int f u^Bone, int f u^poly  ===============================================
//======================================================================================================================================

// Compl = E_PER + 2 E_MIXED + E_AFF 
 //   = (for optimal displacement) E_MIXED + E_AFF
template<typename MatOptConfiguratorType>
class ComplianceEnergyMultipleLoad_RegularizedMaxFunction {
  
  typedef typename MatOptConfiguratorType::RealType RealType;
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::VectorType VectorType;

protected:
  const OptimalDeformSolverMultipleLoad<MatOptConfiguratorType> & _OptDeformFinder;
  mutable RealType _complianceWeightBONE, _complianceWeightPOLYMER;
    mutable std::vector<RealType> _complianceBONE, _compliancePOLYMER;
  
public:
  ComplianceEnergyMultipleLoad_RegularizedMaxFunction( const OptimalDeformSolverMultipleLoad<MatOptConfiguratorType> & OptDeformFinder  ) : _OptDeformFinder ( OptDeformFinder ),
  _complianceBONE ( _OptDeformFinder.getNumLoads() ),
  _compliancePOLYMER ( _OptDeformFinder.getNumLoads() ) { }

  void apply( const VectorType & v, RealType& energy ) const { 
     _OptDeformFinder.updatePhasefield( v );

     for( int i=0; i<_OptDeformFinder.getNumLoads(); ++i ){
        _complianceBONE[i] = 0.0, _compliancePOLYMER[i] = 0.0;
       RealType _complianceBONE_Mixed = -1. * ( _OptDeformFinder.template getRHS<BONE>(i) ).dot( _OptDeformFinder.template getSolutionDisplacementAndMultiplier<BONE>(i) );
       RealType _compliancePOLYMER_Mixed = -1. * ( _OptDeformFinder.template getRHS<POLYMER>(i) ).dot( _OptDeformFinder.template getSolutionDisplacementAndMultiplier<POLYMER>(i) );
        RealType _complianceBONE_Affine = ( _OptDeformFinder.template getHessianLinElastAffinePart<BONE>() * _OptDeformFinder.template getAffineDisplacement<BONE>(i) ).dot( _OptDeformFinder.template getAffineDisplacement<BONE>(i));
        RealType _compliancePOLYMER_Affine = ( _OptDeformFinder.template getHessianLinElastAffinePart<POLYMER>() * _OptDeformFinder.template getAffineDisplacement<POLYMER>(i) ).dot( _OptDeformFinder.template getAffineDisplacement<POLYMER>(i));
          _complianceBONE[i] = _complianceBONE_Mixed + _complianceBONE_Affine;
          _compliancePOLYMER[i] = _compliancePOLYMER_Mixed + _compliancePOLYMER_Affine;
     }

     _complianceWeightBONE = _OptDeformFinder.getMatOptConfigurator().template getWeightFunctionLoad<BONE>().evaluate( _complianceBONE );
     _complianceWeightPOLYMER = _OptDeformFinder.getMatOptConfigurator().template getWeightFunctionLoad<POLYMER>().evaluate( _compliancePOLYMER );

     energy = _OptDeformFinder.getMatOptConfigurator().regularizedMaxFunction( _complianceWeightBONE, _complianceWeightPOLYMER );

  }
  
  RealType getLastComplianceWeightBone( ) const { return _complianceWeightBONE;}
  RealType getLastComplianceWeightPolymer( ) const { return _complianceWeightPOLYMER;}
  RealType getLastComplianceBone(const int i)  const { return _complianceBONE[i];}
  RealType getLastCompliancePolymer(const int i) const { return _compliancePOLYMER[i];}
  
};


//======================================================================================================================================
//================================= Derivative of Jphys  ===============================================================================
//======================================================================================================================================


template<typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType >
class ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM
: public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM<MatOptConfiguratorType,MaterialType> > {
protected: 

    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename MatOptConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
    typedef typename ConfiguratorType::DTContainer DataTypeContainer;

    const OptimalDeformSolverMultipleLoad<MatOptConfiguratorType> &_OptDeformFinder;
    const Material<RealType> &_MaterialBone, &_MaterialPolymer;
    RealType _mu, _lambda;
    const RealType _factorVoidMaterial;
    const int _dimDomain;
    mutable RealType _div_dispAffine;
    mutable DerivativeVectorValuedType _GradSym_dispAffine;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    mutable QuocDiscreteVectorFunctionDefault<ConfiguratorType> *_dispPeriodicPtr;

public:

    ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM ( const OptimalDeformSolverMultipleLoad<MatOptConfiguratorType> & OptDeformFinder ) 
    : QuocFENonlinOpIntegrator<ConfiguratorType, ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM<MatOptConfiguratorType,MaterialType> > (
    OptDeformFinder.getMatOptConfigurator()._conf ),
    _OptDeformFinder ( OptDeformFinder ),
    _MaterialBone ( OptDeformFinder.getMatOptConfigurator()._materialInfo._MaterialBone ),
    _MaterialPolymer ( OptDeformFinder.getMatOptConfigurator()._materialInfo._MaterialPolymer ),
    _factorVoidMaterial( OptDeformFinder.getMatOptConfigurator()._factorVoidMaterial ),
    _dimDomain( OptDeformFinder.getMatOptConfigurator()._conf.dimDomain ),
    _PfPtr ( NULL ), _dispPeriodicPtr ( NULL ) {
        switch( MaterialType ){
            case BONE :
                _mu = _MaterialBone.getMu(); _lambda = _MaterialBone.getLambda();
                break;
            case POLYMER :
                _mu = _MaterialPolymer.getMu(); _lambda = _MaterialPolymer.getLambda();
                break;
            default:
            throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
            break;
        }   
    }

    ~ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM() {
        delete _PfPtr;
        delete _dispPeriodicPtr;
    };

    void apply ( const VectorType &Pf, VectorType &Dest ) const {
        Dest.setZero();
        _OptDeformFinder.updatePhasefield( Pf ); // This computes opt disp u(Pf) and by reference this updated in adjoint problem
  
        std::vector<RealType> _complianceBONE ( _OptDeformFinder.getNumLoads() ),
                              _complianceBONE_Periodic ( _OptDeformFinder.getNumLoads() ),
                              _complianceBONE_Mixed ( _OptDeformFinder.getNumLoads() ),
                              _complianceBONE_Affine ( _OptDeformFinder.getNumLoads() ),
                              _compliancePOLYMER ( _OptDeformFinder.getNumLoads() ),
                              _compliancePOLYMER_Periodic ( _OptDeformFinder.getNumLoads() ),
                              _compliancePOLYMER_Mixed ( _OptDeformFinder.getNumLoads() ),
                              _compliancePOLYMER_Affine ( _OptDeformFinder.getNumLoads() );
        
        for( int i=0; i<_OptDeformFinder.getNumLoads(); ++i ){
            _complianceBONE[i] = 0.0; _compliancePOLYMER[i] = 0.0;
            _complianceBONE_Periodic[i] = 0.0; _complianceBONE_Mixed[i] = 0.0; _complianceBONE_Affine[i] = 0.0;
            _compliancePOLYMER_Periodic[i] = 0.0; _compliancePOLYMER_Mixed[i] = 0.0; _compliancePOLYMER_Affine[i] = 0.0;
            
                _complianceBONE_Periodic[i] = ( _OptDeformFinder.template getHessianLinElast<BONE>() * _OptDeformFinder.template getSolutionDisplacementAndMultiplier<BONE>(i) ).dot( _OptDeformFinder.template getSolutionDisplacementAndMultiplier<BONE>(i));
                _compliancePOLYMER_Periodic[i] = ( _OptDeformFinder.template getHessianLinElast<POLYMER>() * _OptDeformFinder.template getSolutionDisplacementAndMultiplier<POLYMER>(i) ).dot( _OptDeformFinder.template getSolutionDisplacementAndMultiplier<POLYMER>(i));
                _complianceBONE_Mixed[i] = ( _OptDeformFinder.template getRHS<BONE>(i) ).dot( _OptDeformFinder.template getSolutionDisplacementAndMultiplier<BONE>(i) );
                _compliancePOLYMER_Mixed[i] = ( _OptDeformFinder.template getRHS<POLYMER>(i) ).dot( _OptDeformFinder.template getSolutionDisplacementAndMultiplier<POLYMER>(i) );
                _complianceBONE_Affine[i] = ( _OptDeformFinder.template getHessianLinElastAffinePart<BONE>() * _OptDeformFinder.template getAffineDisplacement<BONE>(i) ).dot( _OptDeformFinder.template getAffineDisplacement<BONE>(i));
                _compliancePOLYMER_Affine[i] = ( _OptDeformFinder.template getHessianLinElastAffinePart<POLYMER>() * _OptDeformFinder.template getAffineDisplacement<POLYMER>(i) ).dot( _OptDeformFinder.template getAffineDisplacement<POLYMER>(i));
              _complianceBONE[i] = - 1. * _complianceBONE_Periodic[i] + 2. *_complianceBONE_Mixed[i] - _complianceBONE_Affine[i];
              _compliancePOLYMER[i] = - 1. * _compliancePOLYMER_Periodic[i] + 2. * _compliancePOLYMER_Mixed[i] - _compliancePOLYMER_Affine[i];
        }
        

        RealType complianceBone = _OptDeformFinder.getMatOptConfigurator().template getWeightFunctionLoad<BONE>().evaluate( _complianceBONE );
        RealType compliancePolymer = _OptDeformFinder.getMatOptConfigurator().template getWeightFunctionLoad<POLYMER>().evaluate( _compliancePOLYMER );
        
        RealType fac = 0.0;
        switch( MaterialType ){
            case BONE :
                fac = _OptDeformFinder.getMatOptConfigurator().regularizedMaxFunctionPartialDerivative1( complianceBone, compliancePolymer );
                break;
            case POLYMER :
                fac = _OptDeformFinder.getMatOptConfigurator().regularizedMaxFunctionPartialDerivative2( complianceBone, compliancePolymer );
                break;
            default:
            throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
            break;
        }  
             
        delete _PfPtr;
        _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _OptDeformFinder.getMatOptConfigurator()._conf, _OptDeformFinder.getPhaseFieldPeriodicallyExtended() );
        
        for( int i=0; i<_OptDeformFinder.getNumLoads(); ++i ){
        
            delete _dispPeriodicPtr;
            _dispPeriodicPtr = new QuocDiscreteVectorFunctionDefault<ConfiguratorType> ( _OptDeformFinder.getMatOptConfigurator()._conf, _OptDeformFinder.template getSolutionDisplacementPeriodicallyExtended<MaterialType>(i) );

            QuocDiscreteFunctionDefaultAffine<DataTypeContainer, ConfiguratorType::dimDomain> dispAffine ( _OptDeformFinder.template getAffineDisplacement<MaterialType>(i) );
            _div_dispAffine = dispAffine.getDiv();
            _GradSym_dispAffine = dispAffine.getSymGrad();
        
            RealType facWeightFunction = 0.0;
            switch( MaterialType ){
                case BONE :
                    facWeightFunction = _OptDeformFinder.getMatOptConfigurator().template getWeightFunctionLoad<BONE>().evaluateDerivative( _complianceBONE, i );
                    break;
                case POLYMER :
                    facWeightFunction = _OptDeformFinder.getMatOptConfigurator().template getWeightFunctionLoad<POLYMER>().evaluateDerivative( _compliancePOLYMER, i );
                    break;
                default:
                throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                break;
            }  
            
            QuocFENonlinOpIntegrator< ConfiguratorType, ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM<MatOptConfiguratorType,MaterialType> >::assembleAdd( Dest, facWeightFunction * fac );
        }
    }

    RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {

        RealType div_disp = _dispPeriodicPtr->evaluateDivergenceAtQuadPoint( El, QuadPoint );
        DerivativeVectorValuedType GradSym_disp;
        _dispPeriodicPtr->evaluateSymmetrizedGradientAtQuadPoint( El, QuadPoint, GradSym_disp );

        const RealType Dchi = _OptDeformFinder.getMatOptConfigurator().template approxCharFct_material_Derivative<MaterialType>( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
        const RealType materialfactor = Dchi * ( 1. - _factorVoidMaterial );
        
        RealType aux = 0.0;
        RealType auxPeriodic = 0.0, auxMixed = 0.0, auxAffine = 0.0;
        auxPeriodic = -1. * materialfactor * ( _lambda * div_disp * div_disp + 2. * _mu * aol::ddProd<RealType,DerivativeVectorValuedType>(GradSym_disp, GradSym_disp) );
        auxMixed = -1. *  materialfactor * ( _lambda * div_disp * _div_dispAffine + 2. * _mu * aol::ddProd<RealType,DerivativeVectorValuedType>(GradSym_disp, _GradSym_dispAffine) );
        auxAffine = -1. * materialfactor * ( _lambda * _div_dispAffine * _div_dispAffine + 2. * _mu * aol::ddProd<RealType,DerivativeVectorValuedType>( _GradSym_dispAffine, _GradSym_dispAffine) );
        aux = auxPeriodic + 2. * auxMixed + auxAffine;
        
        return aux; 
    }

};




template<typename MatOptConfiguratorType >
class ComplianceEnergyMultipleLoad_RegularizedMaxFunction_Derivative {
    
protected: 
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
    
  const OptimalDeformSolverMultipleLoad<MatOptConfiguratorType> &_OptDeformFinder;
  
public:
  ComplianceEnergyMultipleLoad_RegularizedMaxFunction_Derivative ( 
   const OptimalDeformSolverMultipleLoad<MatOptConfiguratorType> & OptDeformFinder ) :
        _OptDeformFinder ( OptDeformFinder ) {}
  
  void apply( const VectorType &Pf, VectorType &Dest ) const {
      
      VectorType tmpBoneInM ( Dest.size() ), tmpPolymerInM ( Dest.size() ); 
      ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM<MatOptConfiguratorType,BONE> ( _OptDeformFinder ).apply( Pf, tmpBoneInM );
      ComplianceEnergyMultipleLoad_RegularizedMaxFunction_DerivativeInM<MatOptConfiguratorType,POLYMER> ( _OptDeformFinder ).apply( Pf, tmpPolymerInM );
     
      Dest = tmpBoneInM + tmpPolymerInM;
      
      _OptDeformFinder.getQuocHandler().collabseVectorPeriodicallyAdditive( Dest );
  }
};

};

#endif //__MATERIALOPTIMIZATIONCOMPLIANCEENERGIES_H
