#ifndef __PERIODICHOMOGENIZATIONBONESLINDIFFUSIONENERGYBONES_H
#define __PERIODICHOMOGENIZATIONBONESLINDIFFUSIONENERGYBONES_H

#include <quocDiscreteFunction.h>
#include <quocAffineSymGradIntegrator.h>
#include "BonesMatOptDefines.h"

using namespace quocFE;

namespace shapeOptBonePolymerPeriodicHomogenization{
    
//int a  nabla(u) \cdot eps(u)
template <typename MatOpType, MaterialTypeBonePolymer MaterialType>
class LinDiffusionEnergy :
      public QuocIntegrator< typename MatOpType::ConfiguratorType, LinDiffusionEnergy<MatOpType,MaterialType> > {
protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeMixed LocalMatrixTypeMixed;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineGrad LocalMatrixTypeAffineGrad;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::GlobalAffineGradBaseFuncSet GlobalAffineGradBaseFuncSet;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  
  const MatOpType &_matOpConf;
  const ConfiguratorType &_config;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _diffusionPeriodic;
  DomVecType _Grad_Affine;
  const Material<RealType> & _MaterialBone,  &_MaterialPolymer;
  RealType _diffusionCoefficient;
  const RealType _factorVoidMaterial;
  QuadRuleType _quadRule;
  
public:
  LinDiffusionEnergy ( const MatOpType &matOpConf, const VectorType &material, const VectorType &diffusionPeriodic, const DomVecType &GradDiffusionAffine ) : 
   QuocIntegrator< ConfiguratorType, LinDiffusionEnergy<MatOpType,MaterialType> > ( matOpConf._conf ),
   _matOpConf( matOpConf ),
   _config ( matOpConf._conf ), 
   _pf( _config, material ),
   _diffusionPeriodic( _config, diffusionPeriodic ), _Grad_Affine( GradDiffusionAffine ),
   _MaterialBone ( matOpConf._materialInfo._MaterialBone ),_MaterialPolymer( matOpConf._materialInfo._MaterialPolymer ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) {
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
      
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, QuadPoint )  );
    const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
    DomVecType Grad_Periodic;
    _diffusionPeriodic.evaluateGradientAtQuadPoint( El, QuadPoint, Grad_Periodic );
    return 0.5 * materialfactor * _diffusionCoefficient * (Grad_Periodic + _Grad_Affine).squaredNorm();
  }
};

    

template <typename MatOpType, MaterialTypeBonePolymer MaterialType>
class LinDiffusionHessian :
      public QuocPlusAffineGradMatrixValuedIntegratorBase< typename MatOpType::ConfiguratorType, LinDiffusionHessian<MatOpType,MaterialType> > {
protected:
  typedef typename MatOpType::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeMixedAffineGrad LocalMatrixTypeMixedAffineGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineGrad LocalMatrixTypeAffineGrad;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::GlobalAffineGradBaseFuncSet GlobalAffineGradBaseFuncSet;
  
  const MatOpType &_matOpConf;
  const ConfiguratorType &_config;
  const QuocDiscreteFunctionDefault<ConfiguratorType> _pf;
  DomVecType _Grad_Affine;
  const Material<RealType> & _MaterialBone,  &_MaterialPolymer;
  RealType _diffusionCoefficient;
  const RealType _factorVoidMaterial;
  QuadRuleType _quadRule;
  
public:
  LinDiffusionHessian ( const MatOpType &matOpConf, const VectorType &material ) : 
   QuocPlusAffineGradMatrixValuedIntegratorBase< ConfiguratorType, LinDiffusionHessian<MatOpType,MaterialType> > ( matOpConf._conf ),
   _matOpConf( matOpConf ),
   _config ( matOpConf._conf ), 
   _pf( _config, material ),
   _MaterialBone ( matOpConf._materialInfo._MaterialBone ),_MaterialPolymer( matOpConf._materialInfo._MaterialPolymer ),
   _factorVoidMaterial( matOpConf._factorVoidMaterial ) {
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
   
  void prepareLocalMatrix( const typename ConfiguratorType::ElementType &El, LocalMatrixType &localMatrix ) const {
      
      localMatrix.setZero();
            
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
      const int numDofs = _config.getNumLocalDofs ( El );   
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {

          RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int i = 0; i < numDofs; ++i ) {
              const DomVecType& grad_b_i = bfs.evaluateGradient( i, quadPoint );
              for ( int j = 0; j < numDofs; ++j ) {
                  const DomVecType& grad_b_j = bfs.evaluateGradient( j, quadPoint );
                  localMatrix(j,i) +=  materialfactor * _diffusionCoefficient * grad_b_i.dot(grad_b_j) * bfs.getWeight ( quadPoint );
              }
          }
      }
  }

  void prepareLocalMatrixMixedPart( const typename ConfiguratorType::ElementType &El, LocalMatrixTypeMixedAffineGrad &localMatrix ) const {
      
      localMatrix.setZero();
      
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El);
      const int numDofs = _config.getNumLocalDofs ( El );    
     
      GlobalAffineGradBaseFuncSet globAffBfs;
      DomVecType globalCoord;
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
          _config.getGlobalCoords ( El, _quadRule.getRefCoord( quadPoint ), globalCoord );
          const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int i = 0; i < numDofs; ++i ) {
              const DomVecType& grad_b_i = bfs.evaluateGradient( i, quadPoint );
              for ( int affDof = 0; affDof < globAffBfs.numBaseFuncs; ++affDof ) {
                  DomVecType grad_aff; globAffBfs.evaluateGrad( affDof, globalCoord, grad_aff );
                  localMatrix(affDof,i) +=  materialfactor * _diffusionCoefficient * grad_b_i.dot( grad_aff ) * bfs.getWeight ( quadPoint );
              }
          }
      }
  }

   
  void prepareLocalMatrixAffinePart( const typename ConfiguratorType::ElementType &El, LocalMatrixTypeAffineGrad &localMatrix ) const {
      localMatrix.setZero();   
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet(El); 
      GlobalAffineGradBaseFuncSet globAffBfs;
      DomVecType globalCoord;
      
      for (int quadPoint = 0; quadPoint < _config.maxNumQuadPoints(); ++quadPoint) {
          _config.getGlobalCoords ( El, _quadRule.getRefCoord( quadPoint ), globalCoord );
          const RealType chi = _matOpConf.template approxCharFct_material<MaterialType> ( _pf.evaluateAtQuadPoint( El, quadPoint )  );
          const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
        
          for ( int affDofArg = 0; affDofArg < globAffBfs.numBaseFuncs; ++affDofArg ) {
              DomVecType grad_arg; globAffBfs.evaluateGrad( affDofArg, globalCoord, grad_arg  );
              for ( int affDofDest = 0; affDofDest < globAffBfs.numBaseFuncs; ++affDofDest ) {
                  DomVecType grad_dest; globAffBfs.evaluateGrad( affDofDest, globalCoord, grad_dest  );
                    localMatrix(affDofDest,affDofArg) += materialfactor * _diffusionCoefficient * grad_arg.dot( grad_dest ) * bfs.getWeight ( quadPoint );
              }
          }
      }
  }

};

}//end namespace

#endif
