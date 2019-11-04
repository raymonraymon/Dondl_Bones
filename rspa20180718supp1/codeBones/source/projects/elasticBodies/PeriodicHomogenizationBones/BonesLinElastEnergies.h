#ifndef __PERIODICHOMOGENIZATIONBONESLINELASTENERGYBONES_H
#define __PERIODICHOMOGENIZATIONBONESLINELASTENERGYBONES_H

#include <quocDiscreteFunction.h>
#include <quocAffineSymGradIntegrator.h>
#include "BonesMatOptDefines.h"

using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{
    
//int 2mu eps(u)::eps(u) + lambda div(u) div(u)
template <typename MatOpType, MaterialTypeBonePolymer MaterialType>
class LinElastEnergy :
      public QuocIntegrator< typename MatOpType::ConfiguratorType, LinElastEnergy<MatOpType,MaterialType> > {
protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeMixedAffineSymGrad LocalMatrixTypeMixedAffineSymGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineSymGrad LocalMatrixTypeAffineSymGrad;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::GlobalAffineSymGradBaseFuncSet GlobalAffineSymGradBaseFuncSet;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  
  const MatOpType &_matOpConf;
  const ConfiguratorType &_config;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
  const QuocDiscreteVectorFunctionDefault<ConfiguratorType> _dispPeriodic;
  const QuocDiscreteFunctionDefaultAffine<DataTypeContainer, ConfiguratorType::dimDomain> _dispAffine;
  RealType _div_Affine;
  DerivativeVectorValuedType _GradSym_Affine;
  const Material<RealType> & _MaterialBone,  &_MaterialPolymer;
  RealType _mu, _lambda;
  const RealType _factorVoidMaterial;
  QuadRuleType _quadRule;
  
public:
  LinElastEnergy ( const MatOpType &matOpConf, const VectorType &material, const VectorType &dispPeriodic, const VectorType &dispAffine ) : 
   QuocIntegrator< ConfiguratorType, LinElastEnergy<MatOpType,MaterialType> > ( matOpConf._conf ),
   _matOpConf( matOpConf ),
   _config ( matOpConf._conf ), 
   _pf( _config, material ),
   _dispPeriodic( _config, dispPeriodic ), _dispAffine( dispAffine ),
   _MaterialBone ( matOpConf._materialInfo._MaterialBone ),_MaterialPolymer( matOpConf._materialInfo._MaterialPolymer ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) {
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
      
      _GradSym_Affine = _dispAffine.getSymGrad();
      _div_Affine = _dispAffine.getDiv();
   }
      
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, QuadPoint )  );
    const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
    
    RealType div_Periodic = _dispPeriodic.evaluateDivergenceAtQuadPoint( El, QuadPoint );
    DerivativeVectorValuedType GradSym_Periodic;
    _dispPeriodic.evaluateSymmetrizedGradientAtQuadPoint( El, QuadPoint, GradSym_Periodic );
    
    return materialfactor * 
    ( 2. * _mu * aol::ddProd<RealType,DerivativeVectorValuedType>( GradSym_Periodic + _GradSym_Affine, GradSym_Periodic + _GradSym_Affine  ) 
     + _lambda * aol::Sqr( div_Periodic + _div_Affine ) );
  }
};


template <typename MatOpType, MaterialTypeBonePolymer MaterialType>
class LinElastHessianApplied :
      public QuocIntegrator< typename MatOpType::ConfiguratorType, LinElastHessianApplied<MatOpType,MaterialType> > {
protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeMixedAffineSymGrad LocalMatrixTypeMixedAffineSymGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineSymGrad LocalMatrixTypeAffineSymGrad;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::GlobalAffineSymGradBaseFuncSet GlobalAffineSymGradBaseFuncSet;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  
  const MatOpType &_matOpConf;
  const ConfiguratorType &_config;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
  const QuocDiscreteVectorFunctionDefault<ConfiguratorType> _dispPeriodic1, _dispPeriodic2;
  const QuocDiscreteFunctionDefaultAffine<DataTypeContainer, ConfiguratorType::dimDomain> _dispAffine1, _dispAffine2;
  RealType _div_Affine1, _div_Affine2;
  DerivativeVectorValuedType _GradSym_Affine1, _GradSym_Affine2;
  const Material<RealType> & _MaterialBone,  &_MaterialPolymer;
  RealType _mu, _lambda;
  const RealType _factorVoidMaterial;
  QuadRuleType _quadRule;
  
public:
  LinElastHessianApplied ( const MatOpType &matOpConf, const VectorType &material, 
                           const VectorType &dispPeriodic1, const VectorType &dispAffine1,
                           const VectorType &dispPeriodic2, const VectorType &dispAffine2 ) : 
   QuocIntegrator< ConfiguratorType, LinElastHessianApplied<MatOpType,MaterialType> > ( matOpConf._conf ),
   _matOpConf( matOpConf ),
   _config ( matOpConf._conf ), 
   _pf( _config, material ),
   _dispPeriodic1( _config, dispPeriodic1 ), _dispAffine1( dispAffine1 ),
   _dispPeriodic2( _config, dispPeriodic2 ), _dispAffine2( dispAffine2 ),
   _MaterialBone ( matOpConf._materialInfo._MaterialBone ),_MaterialPolymer( matOpConf._materialInfo._MaterialPolymer ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) {
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
      
      _GradSym_Affine1 = _dispAffine1.getSymGrad(); _GradSym_Affine2 = _dispAffine2.getSymGrad();
      _div_Affine1 = _dispAffine1.getDiv(); _div_Affine2 = _dispAffine2.getDiv();
   }
      
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, QuadPoint )  );
    const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
    
    const RealType div_Periodic1 = _dispPeriodic1.evaluateDivergenceAtQuadPoint( El, QuadPoint );
    const RealType div_Periodic2 = _dispPeriodic2.evaluateDivergenceAtQuadPoint( El, QuadPoint );
    DerivativeVectorValuedType GradSym_Periodic1, GradSym_Periodic2;
    _dispPeriodic1.evaluateSymmetrizedGradientAtQuadPoint( El, QuadPoint, GradSym_Periodic1 );
    _dispPeriodic2.evaluateSymmetrizedGradientAtQuadPoint( El, QuadPoint, GradSym_Periodic2 );
    
    return materialfactor * 
    ( 2. * _mu * aol::ddProd<RealType,DerivativeVectorValuedType>( GradSym_Periodic1 + _GradSym_Affine1, GradSym_Periodic2 + _GradSym_Affine2  ) 
     + _lambda * ( div_Periodic1 + _div_Affine1 ) * ( div_Periodic2 + _div_Affine2 ) );
  }
};
    

template <typename MatOpType, MaterialTypeBonePolymer MaterialType>
class LinElastHessian :
      public QuocPlusAffineSymGradBlockMatrixValuedIntegratorBase< typename MatOpType::ConfiguratorType, LinElastHessian<MatOpType,MaterialType> > {
protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeMixedAffineSymGrad LocalMatrixTypeMixedAffineSymGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineSymGrad LocalMatrixTypeAffineSymGrad;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::GlobalAffineSymGradBaseFuncSet GlobalAffineSymGradBaseFuncSet;
  
  const MatOpType &_matOpConf;
  const ConfiguratorType &_config;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
  const Material<RealType> & _MaterialBone,  &_MaterialPolymer;
  RealType _mu, _lambda;
  const RealType _factorVoidMaterial;
  QuadRuleType _quadRule;
  
public:
  LinElastHessian ( const MatOpType &matOpConf, const VectorType &material ) : 
   QuocPlusAffineSymGradBlockMatrixValuedIntegratorBase< ConfiguratorType, LinElastHessian<MatOpType,MaterialType> > ( matOpConf._conf ),
   _matOpConf( matOpConf ),
   _config ( matOpConf._conf ), 
   _pf( _config, material ),
   _MaterialBone ( matOpConf._materialInfo._MaterialBone ),_MaterialPolymer( matOpConf._materialInfo._MaterialPolymer ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) {
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
   
   
  void prepareLocalMatrix( const typename ConfiguratorType::ElementType &El, LocalMatrixType (&localMatrix)[ConfiguratorType::dimDomain][ConfiguratorType::dimDomain]) const {
      
      for (int argComp = 0; argComp < ConfiguratorType::dimDomain; ++argComp)
          for (int destComp = 0; destComp < ConfiguratorType::dimDomain; ++destComp)
              localMatrix[argComp][destComp].setZero();
            
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
      const int numDofs = _config.getNumLocalDofs ( El );   
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {

          RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int i = 0; i < numDofs; ++i ) {
              const DomVecType& grad_b_i = bfs.evaluateGradient( i, quadPoint );
              for ( int j = 0; j < numDofs; ++j ) {
                  const DomVecType& grad_b_j = bfs.evaluateGradient( j, quadPoint );
                  for(int argComp=0; argComp < ConfiguratorType::dimDomain; ++argComp){
                      localMatrix[argComp][argComp](j,i) +=  materialfactor * _mu * grad_b_i.dot(grad_b_j) * bfs.getWeight ( quadPoint );
                      for(int destComp=0; destComp < ConfiguratorType::dimDomain; ++destComp ){
                        localMatrix[argComp][destComp](j,i) += materialfactor * _mu * grad_b_i[destComp] * grad_b_j[argComp] * bfs.getWeight ( quadPoint );
                        localMatrix[argComp][destComp](j,i) += materialfactor * _lambda * grad_b_i[argComp] * grad_b_j[destComp] * bfs.getWeight ( quadPoint );
                      }
                  }
              }
          }
      }
  }

  void prepareLocalMatrixMixedPart( const typename ConfiguratorType::ElementType &El, LocalMatrixTypeMixedAffineSymGrad (&localMatrix)[ConfiguratorType::dimDomain] ) const {
      for (int comp = 0; comp < ConfiguratorType::dimDomain; ++comp) localMatrix[comp].setZero();
      
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
      const int numDofs = _config.getNumLocalDofs ( El );    
     
      GlobalAffineSymGradBaseFuncSet globAffBfs;
      DomVecType globalCoord;
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
          _config.getGlobalCoords ( El, _quadRule.getRefCoord( quadPoint ), globalCoord );
          const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int i = 0; i < numDofs; ++i ) {
              const DomVecType& grad_b_i = bfs.evaluateGradient( i, quadPoint );
              for ( int affDof = 0; affDof < globAffBfs.numBaseFuncs; ++affDof ) {
                  DerivativeVectorValuedType symgrad; globAffBfs.evaluateSymGrad( affDof, globalCoord, symgrad );
                  RealType div = globAffBfs.evaluateDiv( affDof, globalCoord );
                  for(int argComp=0; argComp < ConfiguratorType::dimDomain; ++argComp ){
                    localMatrix[argComp](affDof,i) +=  materialfactor * _mu * grad_b_i.dot( symgrad.col(argComp) ) * bfs.getWeight ( quadPoint );
                    for(int destComp=0; destComp < ConfiguratorType::dimDomain; ++destComp ){
                      localMatrix[destComp](affDof,i) += materialfactor * _mu * grad_b_i[argComp] * symgrad(argComp,destComp) * bfs.getWeight ( quadPoint );
                      localMatrix[destComp](affDof,i) += materialfactor * _lambda * grad_b_i[destComp] * symgrad(argComp,argComp) * bfs.getWeight ( quadPoint );
                    }
                  }
              }
          }
      }
  }

   
  void prepareLocalMatrixAffinePart( const typename ConfiguratorType::ElementType &El, LocalMatrixTypeAffineSymGrad &localMatrix ) const {
      localMatrix.setZero();   
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El); 
      GlobalAffineSymGradBaseFuncSet globAffBfs;
      DomVecType globalCoord;
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
          _config.getGlobalCoords ( El, _quadRule.getRefCoord( quadPoint ), globalCoord );
          const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int affDofArg = 0; affDofArg < globAffBfs.numBaseFuncs; ++affDofArg ) {
              DerivativeVectorValuedType symgrad_arg; globAffBfs.evaluateSymGrad( affDofArg, globalCoord, symgrad_arg  );
              RealType div_arg = globAffBfs.evaluateDiv( affDofArg, globalCoord );
              for ( int affDofDest = 0; affDofDest < globAffBfs.numBaseFuncs; ++affDofDest ) {
                  DerivativeVectorValuedType symgrad_dest; globAffBfs.evaluateSymGrad( affDofDest, globalCoord, symgrad_dest  );
                  RealType div_dest = globAffBfs.evaluateDiv( affDofDest, globalCoord );
                      localMatrix(affDofDest,affDofArg) += materialfactor * 2. * _mu * aol::ddProd<RealType,DerivativeVectorValuedType>( symgrad_arg, symgrad_dest ) * bfs.getWeight ( quadPoint );
                      localMatrix(affDofDest,affDofArg) += materialfactor * _lambda * div_arg * div_dest * bfs.getWeight ( quadPoint );
              }
          }
      }
  }

};

}//end namespace

#endif