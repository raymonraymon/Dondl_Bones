#ifndef __BONESMATERIALOPTIMIZATIONDEFINES__H
#define __BONESMATERIALOPTIMIZATIONDEFINES__H

#include <materialOptimizationDefines.h>
#include <multipleLoad.h>

#define APPROXCHARFCT_BONES_FOURTHORDER

enum MaterialTypeBonePolymer {BONE, POLYMER};

template<typename RealType>
class materialOptBonesInfo{
    
private :
     MaterialDictionary<RealType> materialDict;
    
public:
    Material<RealType> _MaterialBone, _MaterialPolymer;
    
    template<typename ParameterParserType>
    materialOptBonesInfo ( const ParameterParserType &parser ) :
    materialDict ( ),
    _MaterialBone ( ), _MaterialPolymer (  ){ 
    
        if( parser.template get<int>( "Material.MaterialFromDictionary" ) == 1 ){
            _MaterialBone.set( materialDict.getMaterialFromDictionary( parser.template get<std::string>( "Material.MaterialBone" ) ) );
            _MaterialPolymer.set( materialDict.getMaterialFromDictionary( parser.template get<std::string>( "Material.MaterialPolymer" ) ) );
        }else{
            _MaterialBone.set(  Material<RealType> ( "Bone", 1.0, parser.template get<RealType> ( "Material.DiffusionCoefficient_Bone" ), parser.template get<RealType> ( "Material.ElastModulus_Bone" ), parser.template get<RealType> ( "Material.PoissonRatio_Bone" ) ) );
            _MaterialPolymer.set(  Material<RealType> ( "Polymer", 1.0, parser.template get<RealType> ( "Material.DiffusionCoefficient_Polymer" ), parser.template get<RealType> ( "Material.ElastModulus_Polymer" ), parser.template get<RealType> ( "Material.PoissonRatio_Polymer" ) ) );
        }
    }
    
    
    materialOptBonesInfo ( const Material<RealType> &MaterialBone, const Material<RealType> &MaterialPolymer ) :
    materialDict ( ), _MaterialBone ( MaterialBone ), _MaterialPolymer ( MaterialPolymer  ){ 
    }
};




template < typename _ConfiguratorType >
class QuocMaterialConfiguratorBones {
  public :  
      
   typedef _ConfiguratorType ConfiguratorType;
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    
   const ConfiguratorType &_conf;
   const RealType _regParameterMaxFunction;
   const RealType _factorDoubleWell;
  
   QuocMaterialConfiguratorBones ( const ConfiguratorType &conf, const RealType regParameterMaxFunction, const RealType factorDoubleWell = 9. / 16. ) : 
   _conf ( conf ), _regParameterMaxFunction ( regParameterMaxFunction ), _factorDoubleWell ( factorDoubleWell ) {}


#ifdef APPROXCHARFCT_BONES_FOURTHORDER
    
    template<MaterialTypeBonePolymer MaterialType>
    RealType approxCharFct_vol ( const RealType v ) const { 
      switch( MaterialType ){
        case BONE :
           return 0.5 * ( v + 1.0 );
            break;
        case POLYMER :
            return  0.5 * ( 1.0 - v );
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    RealType approxCharFct_vol_Derivative ( const RealType v ) const {
      switch( MaterialType ){
        case BONE :
            return 0.5;
            break;
        case POLYMER :
            return -0.5;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
    }
    
    RealType approxCharFct_vol ( const RealType v ) const { return aol::Sqr( aol::Sqr( v + 1.0 ) ) / 16.; }
    RealType approxCharFct_vol_Derivative ( const RealType v ) const { return 0.25 * aol::Cub( v + 1.0 );}
    RealType approxCharFct_vol_SecondDerivative ( const RealType v ) const { return 0.75 * aol::Sqr( v + 1.0 );}  
    
    template<MaterialTypeBonePolymer MaterialType>
    RealType approxCharFct_material ( const RealType v ) const { 
      switch( MaterialType ){
        case BONE :
           return aol::Sqr( aol::Sqr( v + 1.0 ) ) / 16.;
            break;
        case POLYMER :
            return aol::Sqr( aol::Sqr( v - 1.0 ) ) / 16.;
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
    }
    
    template<MaterialTypeBonePolymer MaterialType>
    RealType approxCharFct_material_Derivative ( const RealType v ) const {
      switch( MaterialType ){
        case BONE :
            return 0.25 * aol::Cub( v + 1.0 );
            break;
        case POLYMER :
            return 0.25 * aol::Cub( v - 1.0 );
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
    }
#endif


    RealType doubleWell ( const RealType v ) const { return _factorDoubleWell * aol::Sqr( v * v - 1.0 );}
    RealType doubleWellDerivative ( const RealType v ) const { return _factorDoubleWell * 4.0 * ( v * v - 1.0 ) * v;}
    RealType doubleWellSecondDerivative ( const RealType v ) const { return _factorDoubleWell * 4.0 * ( 3. * v * v - 1.0 );}
    
    
    
    RealType regularizedMaxFunction ( const RealType &x, const RealType &y ) const {
     return 0.5 * ( x + y + std::sqrt( aol::Sqr(x-y) + _regParameterMaxFunction ) );   
    }
    RealType regularizedMaxFunctionPartialDerivative1 ( const RealType &x, const RealType &y ) const {
     return 0.5 * ( 1. + (x - y) / std::sqrt( aol::Sqr(x-y) + _regParameterMaxFunction ) );   
    }
    RealType regularizedMaxFunctionPartialDerivative2 ( const RealType &x, const RealType &y ) const {
     //return 0.5 * ( 1. + (y - x) / std::sqrt( aol::Sqr(x-y) + _regParameterMaxFunction ) );  
     return this->regularizedMaxFunctionPartialDerivative1(y,x);
    }

};


template < typename _ConfiguratorType>
class QuocMaterialOptimizationConfiguratorBones : public QuocMaterialConfiguratorBones<_ConfiguratorType> {
  
public :
  
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename DataTypeContainer::ParameterParserType ParameterParserType;
  
  const materialOptBonesInfo<RealType> _materialInfo;
  const RealType _factorVoidMaterial;
  RealType _factorComplianceCost, _factorInterfaceCost;
  const RealType _epsInterfaceLength;
  const RealType _epsFactor;
  
protected:
    WeightFunctionMultipleLoad<RealType> _weightFctLoad_BONE, _weightFctLoad_POLYMER;
 
public:
  
  //used e.g. if one is only interested in optimal deformation
  QuocMaterialOptimizationConfiguratorBones ( const ParameterParserType &parser, const ConfiguratorType &conf ) :
    QuocMaterialConfiguratorBones<_ConfiguratorType> ( conf, 1.0 ),
    _materialInfo ( parser ),
    _factorVoidMaterial( parser.template get<double>( "MaterialOptimization.factorVoidMaterial" ) ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( 1. ), _epsFactor ( 1. ) {
        cout << endl << endl << "WARNING: this constructor of QuocMaterialOptimizationConfigurator should not be used for material optimization, since eps=1" << endl << endl;
    }

  QuocMaterialOptimizationConfiguratorBones ( const ParameterParserType &parser, const ConfiguratorType &conf, const RealType epsInterfaceLength, const RealType epsFactor ) :
    QuocMaterialConfiguratorBones<_ConfiguratorType> ( conf, parser.template get<double> ("MaterialOptimization.regParameterMaxFunction" ), parser.template get<double> ("MaterialOptimization.factorDoubleWell" ) ),
    _materialInfo ( parser ),
    _factorVoidMaterial( parser.template get<double>( "MaterialOptimization.factorVoidMaterial" ) ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( epsInterfaceLength ), _epsFactor ( epsFactor ) {}
    
  QuocMaterialOptimizationConfiguratorBones ( const ParameterParserType &parser, const ConfiguratorType &conf, const RealType epsInterfaceLength, const RealType epsFactor, 
                                              const std::vector<RealType> &weightFct_Bone_weightVec, const std::vector<RealType> &weightFct_Polymer_weightVec ) :
    QuocMaterialConfiguratorBones<_ConfiguratorType> ( conf, parser.template get<double> ("MaterialOptimization.regParameterMaxFunction" ), parser.template get<double> ("MaterialOptimization.factorDoubleWell" ) ),
    _materialInfo ( parser ),
    _factorVoidMaterial( parser.template get<double>( "MaterialOptimization.factorVoidMaterial" ) ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( epsInterfaceLength ), _epsFactor ( epsFactor ), 
    _weightFctLoad_BONE( parser.template get<RealType>("MaterialOptimization.weightFct_Bone_weight"), parser.template get<RealType>("MaterialOptimization.weightFct_Bone_p"), parser.template get<RealType>("MaterialOptimization.weightFct_Bone_q"), weightFct_Bone_weightVec ), 
    _weightFctLoad_POLYMER( parser.template get<RealType>("MaterialOptimization.weightFct_Polymer_weight"), parser.template get<RealType>("MaterialOptimization.weightFct_Polymer_p"), parser.template get<RealType>("MaterialOptimization.weightFct_Polymer_q"), weightFct_Polymer_weightVec ) {}

    
  template<MaterialTypeBonePolymer MaterialType>
    const WeightFunctionMultipleLoad<RealType> &getWeightFunctionLoad( ) const{ 
     switch( MaterialType ){
        case BONE :    return _weightFctLoad_BONE;    break;
        case POLYMER : return _weightFctLoad_POLYMER; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
      }
    }
    
};



#endif
