#ifndef __QUOCPHBONESENERGIESWITHMATERIAL_H
#define __QUOCPHBPMESENERGIESWITHMATERIAL_H

#include <general.h>
#include <energyDefines.h>
#include <quocIntegrator.h>
#include <quocHandler.h>
#include "BonesMatOptDefines.h" //for approxCharFct, doubleWell

using namespace quocFE;
    

namespace shapeOptBonePolymerPeriodicHomogenization{    

//==========================================================================================================================
// volume \int_domain chi(m)
//==========================================================================================================================
template <typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class PfOp_Volume :
public QuocIntegrator < typename MatOptConfiguratorType::ConfiguratorType, PfOp_Volume<MatOptConfiguratorType,MaterialType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  
  const MatOptConfiguratorType & _matOptConf;
  mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
  
  public:
    PfOp_Volume ( const MatOptConfiguratorType & matOptConf ) : 
     QuocIntegrator<ConfiguratorType, PfOp_Volume<MatOptConfiguratorType,MaterialType> > ( matOptConf._conf ),
     _matOptConf ( matOptConf ), 
     _PfPtr ( NULL ) {}
     
  ~PfOp_Volume() {
     delete _PfPtr;
  };
  
  void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocIntegrator< ConfiguratorType, PfOp_Volume<MatOptConfiguratorType,MaterialType>  >::assembleAdd( energy );
  }

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     return _matOptConf.template approxCharFct_vol<MaterialType>( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


// D Area(v) (hat v) = \int chi'(v) (hat v)
template <typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class PfOp_VolumeDerivative : public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_VolumeDerivative<MatOptConfiguratorType,MaterialType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    
    const MatOptConfiguratorType & _matOptConf;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_VolumeDerivative ( const MatOptConfiguratorType & matOptConf ) 
  : QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_VolumeDerivative<MatOptConfiguratorType,MaterialType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _PfPtr ( NULL ) {}
     
  ~PfOp_VolumeDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFENonlinOpIntegrator< ConfiguratorType, PfOp_VolumeDerivative<MatOptConfiguratorType,MaterialType> >::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return _matOptConf.template approxCharFct_vol_Derivative<MaterialType> ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


// D^2 Area(v) (hat v) (dhat v) = \int chi''(v) (hat v)(dhat v)
template <typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class PfOp_VolumeSecondDerivative : public QuocFELinWeightedMassIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOp_VolumeSecondDerivative<MatOptConfiguratorType,MaterialType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    
    const MatOptConfiguratorType & _matOptConf;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOp_VolumeSecondDerivative ( const MatOptConfiguratorType & matOptConf ) 
  : QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<MatOptConfiguratorType,MaterialType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), 
    _PfPtr ( NULL ) {}
     
  ~PfOp_VolumeSecondDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOp_VolumeSecondDerivative<MatOptConfiguratorType,MaterialType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      return _matOptConf.template approxCharFct_vol_SecondDerivative<MaterialType> ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) );
  }
};


// int_S chi(v)
template <typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class VolumeConstraint : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  
    protected:
        const MatOptConfigurator &_matOptConf;

    public:
    VolumeConstraint ( const MatOptConfigurator &matOptConf ) : _matOptConf ( matOptConf ){ } 
      
    void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
         PfOp_Volume<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( v, energy );
    }
    
    void evaluateJacobian( const VectorType & v, VectorType& Deriv ) const {
         PfOp_VolumeDerivative<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( v, Deriv );
    }
    
    void evaluateHessian( const VectorType & v, SparseMatrixType& Hessian ) const {
         PfOp_VolumeSecondDerivative<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( v, Hessian );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
      
};


template <typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class VolumeConstraintPeriodicBC : public aol::NonlinearEnergyOp< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  
  protected:
        const MatOptConfigurator &_matOptConf;
        const QuocHandler<ConfiguratorType> &_quocHandler;

  public:
  VolumeConstraintPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler )  : _matOptConf ( matOptConf ), _quocHandler ( quocHandler ) { } 
      
    void evaluateEnergy( const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOp_Volume<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( vCollabsedAndExtended, energy );
    }
    
    void evaluateJacobian( const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOp_VolumeDerivative<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( vCollabsedAndExtended, Deriv );
         _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
    void evaluateTripletListHessianSym( const VectorType &Arg, std::vector<typename ConfiguratorType::TripletType> & tripletListHessian ) const {
        throw std::invalid_argument( aol::strprintf ( "HessianTriplet not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
      
};







 
//==========================================================================================================================
// barycenter_i(m) = (1/int_domain chi(m) ) \int_domain chi(m) x_i  -  c_i
//==========================================================================================================================
template <typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class PfOpBones_Barycenter :
public QuocIntegrator < typename MatOptConfiguratorType::ConfiguratorType, PfOpBones_Barycenter<MatOptConfiguratorType,MaterialType> >{
  
protected :
  typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
  
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::PointType PointType;
  
  const MatOptConfiguratorType & _matOptConf;
  const int _direction; 
  const RealType _c;
  mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
  
  public:
    PfOpBones_Barycenter ( const MatOptConfiguratorType & matOptConf, const int direction, const RealType c ) : 
     QuocIntegrator<ConfiguratorType, PfOpBones_Barycenter<MatOptConfiguratorType,MaterialType> > ( matOptConf._conf ),
     _matOptConf ( matOptConf ), _direction ( direction ), _c ( c ),
     _PfPtr ( NULL ) {}
     
  ~PfOpBones_Barycenter() {
     delete _PfPtr;
  };
  
  void apply( const VectorType & v, RealType& energy ) const {   
      energy = 0.;
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocIntegrator< ConfiguratorType, PfOpBones_Barycenter<MatOptConfiguratorType,MaterialType>  >::assembleAdd( energy );
  }

  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint) const{
     PointType GlobalCoord; _matOptConf._conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
//      return _matOptConf.template approxCharFct_vol<MaterialType>( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
     return _matOptConf.template approxCharFct_material<MaterialType>( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


template <typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class PfOpBones_BarycenterDerivative : public QuocFENonlinOpIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOpBones_BarycenterDerivative<MatOptConfiguratorType,MaterialType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::PointType PointType;
    
    const MatOptConfiguratorType & _matOptConf;
    const int _direction; 
    const RealType _c;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOpBones_BarycenterDerivative ( const MatOptConfiguratorType & matOptConf, const int direction, const RealType c  ) 
  : QuocFENonlinOpIntegrator< ConfiguratorType, PfOpBones_BarycenterDerivative<MatOptConfiguratorType,MaterialType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _direction ( direction ), _c ( c ), _PfPtr ( NULL ) {}
     
  ~PfOpBones_BarycenterDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, VectorType& Deriv ) const {   
      Deriv.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFENonlinOpIntegrator< ConfiguratorType, PfOpBones_BarycenterDerivative<MatOptConfiguratorType,MaterialType> >::assembleAdd( Deriv );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      PointType GlobalCoord; _matOptConf._conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
//       return _matOptConf.template approxCharFct_vol_Derivative<MaterialType> ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
      return _matOptConf.template approxCharFct_material_Derivative<MaterialType> ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


// D^2 Area(v) (hat v) (dhat v) = \int chi''(v) (hat v)(dhat v)
template <typename MatOptConfiguratorType, MaterialTypeBonePolymer MaterialType>
class PfOpBones_BarycenterSecondDerivative : public QuocFELinWeightedMassIntegrator< typename MatOptConfiguratorType::ConfiguratorType, PfOpBones_BarycenterSecondDerivative<MatOptConfiguratorType,MaterialType> > {
  
  protected:
    typedef typename MatOptConfiguratorType::ConfiguratorType ConfiguratorType;
    
    typedef typename ConfiguratorType::RealType RealType;
    typedef typename ConfiguratorType::VectorType VectorType;
    typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
    typedef typename ConfiguratorType::PointType PointType;
    
    const MatOptConfiguratorType & _matOptConf;
    const int _direction; 
    const RealType _c;
    mutable QuocDiscreteFunctionDefault<ConfiguratorType> *_PfPtr;
    
  public:
  PfOpBones_BarycenterSecondDerivative ( const MatOptConfiguratorType & matOptConf, const int direction, const RealType c  ) 
  : QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOpBones_BarycenterSecondDerivative<MatOptConfiguratorType,MaterialType>  > ( matOptConf._conf ),
    _matOptConf ( matOptConf ), _direction ( direction ), _c ( c ),  _PfPtr ( NULL ) {}
     
  ~PfOpBones_BarycenterSecondDerivative() {
     delete _PfPtr;
  };
  
   void apply( const VectorType & v, SparseMatrixType& Hessian ) const {   
      Hessian.setZero();
      delete _PfPtr;
     _PfPtr = new QuocDiscreteFunctionDefault<ConfiguratorType> ( _matOptConf._conf, v );
      QuocFELinWeightedMassIntegrator< ConfiguratorType, PfOpBones_BarycenterSecondDerivative<MatOptConfiguratorType,MaterialType> >::assemble( Hessian );
  }
  
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
      PointType GlobalCoord; _matOptConf._conf.getGlobalCoords ( El, QuadPoint, GlobalCoord );
//       return _matOptConf.template approxCharFct_vol_SecondDerivative<MaterialType> ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
      return _matOptConf.template approxCharFct_material_SecondDerivative<MaterialType> ( _PfPtr->evaluateAtQuadPoint( El, QuadPoint ) ) * ( GlobalCoord[_direction] - _c);
  }
};


// int_S chi(v)
template <typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class BaryCenterConstraintBones : public aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
    protected:
        const MatOptConfigurator &_matOptConf;
        const PointType &_c;

    public:
    BaryCenterConstraintBones ( const MatOptConfigurator &matOptConf, const PointType &c ) 
    : aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  > ( c.size() ),
      _matOptConf ( matOptConf ), _c ( c ) { } 
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
         PfOpBones_Barycenter<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( v, energy );
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
         PfOpBones_BarycenterDerivative<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( v, Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & v, SparseMatrixType& Hessian ) const {
         PfOpBones_BarycenterSecondDerivative<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( v, Hessian );
    }
      
};


template <typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class BaryCenterConstraintBonesPeriodicBC : public aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
protected:
  const MatOptConfigurator &_matOptConf;
  const QuocHandler<ConfiguratorType> & _quocHandler;
  const PointType &_c;

public:
  BaryCenterConstraintBonesPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler, const PointType &c ) 
  : aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  > ( c.size() ),
  _matOptConf ( matOptConf ), _quocHandler ( quocHandler ), _c ( c ) { }
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOpBones_Barycenter<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, energy );
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         PfOpBones_BarycenterDerivative<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, Deriv );
         _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & v, SparseMatrixType& Hessian ) const {
        throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
//         VectorType vCollabsedAndExtended( v );
//         _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
//         _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
//          PfOpBones_BarycenterSecondDerivative<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint, _c[numConstraint] ).apply( vCollabsedAndExtended, Hessian );
//          // TODO _quocHandler.collabseMatrixPeriodicallyAdditive( Hessian );
    }
    
    const PointType& getBarycenterPoint( ) const {return _c;};
      
};












template <typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class VolumeAndBarycenterConstraintPeriodicBC : public aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  >
{
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename MatOptConfigurator::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::SparseMatrixType SparseMatrixType;
  typedef typename ConfiguratorType::PointType PointType;
  
protected:
  const MatOptConfigurator &_matOptConf;
  const QuocHandler<ConfiguratorType> & _quocHandler;
  const PointType &_c;

  public:
  VolumeAndBarycenterConstraintPeriodicBC ( const MatOptConfigurator &matOptConf, const QuocHandler<ConfiguratorType> &quocHandler, const PointType &c )  :
  aol::NonlinearConstraintOps< typename MatOptConfigurator::ConfiguratorType::DTContainer  > ( 1 + c.size() ),
  _matOptConf ( matOptConf ), _quocHandler ( quocHandler ), _c ( c ) { } 
      
    void evaluateEnergy( const int numConstraint, const VectorType & v, RealType& energy ) const {
         energy = 0;
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
        
        if( numConstraint == 0 )  PfOp_Volume<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( vCollabsedAndExtended, energy );
        else PfOpBones_Barycenter<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint - 1, _c[numConstraint-1] ).apply( vCollabsedAndExtended, energy );
         
    }
    
    void evaluateJacobian( const int numConstraint, const VectorType & v, VectorType& Deriv ) const {
        VectorType vCollabsedAndExtended( v );
        _quocHandler.collabseVectorPeriodically( vCollabsedAndExtended );
        _quocHandler.extendVectorPeriodically( vCollabsedAndExtended );
         
        if( numConstraint == 0 ) PfOp_VolumeDerivative<MatOptConfigurator,MaterialType> ( _matOptConf ).apply( vCollabsedAndExtended, Deriv );
        else PfOpBones_BarycenterDerivative<MatOptConfigurator,MaterialType> ( _matOptConf, numConstraint - 1, _c[numConstraint-1] ).apply( vCollabsedAndExtended, Deriv );
        
        _quocHandler.collabseVectorPeriodicallyAdditive( Deriv );
    }
    
    void evaluateHessian( const int numConstraint, const VectorType & /*v*/, SparseMatrixType& /*Hessian*/ ) const {
         throw std::invalid_argument( aol::strprintf ( "Hessian is so far not implemented. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    
//     const PointType& getBarycenterPoint( ) const {return _c;};
      
};

}//end namespace


#endif