#ifndef __QUOCMATERIALOPTIMIZATIONDEFINES__H
#define __QUOCMATERIALOPTIMIZATIONDEFINES__H

// #define QUOC_APPROXCHARFCT_WIRTH
// #define APPROXCHARFCT_1DHOMOGENIZATION
#define QUOC_APPROXCHARFCT_FOURTHORDER

#include <materialOptimizationDefines.h>

template<typename RealType>
class QuocMaterialOptInfo{
    
private :
     MaterialDictionary<RealType> materialDict;
    
public:
    Material<RealType> _HardMaterial, _SoftMaterial;
    
    template<typename ParameterParserType>
    QuocMaterialOptInfo ( const ParameterParserType &parser ) :
    materialDict ( ),
    _HardMaterial ( ), _SoftMaterial ( ){ 
        if( parser.template get<int>( "Material.MaterialFromDictionary" ) == 1 ){
            //cout << "load materials from dictonary" << endl;
            _HardMaterial.set( materialDict.getMaterialFromDictionary( parser.template get<std::string>( "Material.MaterialHard" ) ) );
            _SoftMaterial.set( materialDict.getMaterialFromDictionary( parser.template get<std::string>( "Material.MaterialSoft" ) ) );
        }else{
            //cout << "read material constants from parser" << endl;
            _HardMaterial.set(  Material<RealType> ( "Hard", 1.0, parser.template get<RealType> ( "Material.ElastModulus_Hard" ), parser.template get<RealType> ( "Material.PoissonRatio_Hard" ) ) );
            _SoftMaterial.set(  Material<RealType> ( "Soft", 1.0, parser.template get<RealType> ( "Material.ElastModulus_Soft" ), parser.template get<RealType> ( "Material.PoissonRatio_Soft" ) ) );
        }
        //_HardMaterial.print(); _SoftMaterial.print();
    }
    
    
    QuocMaterialOptInfo ( const Material<RealType> &MaterialHard, const Material<RealType> &MaterialSoft ) :
    materialDict ( ), _HardMaterial ( MaterialHard ), _SoftMaterial ( MaterialSoft  ){ 
        _HardMaterial.print(); _SoftMaterial.print();
    }
};


template < typename _ConfiguratorType >
class QuocMaterialConfigurator {
  public :  
      
   typedef _ConfiguratorType ConfiguratorType;
   typedef typename ConfiguratorType::RealType RealType;
   typedef typename ConfiguratorType::DTContainer DataTypeContainer;
    
   const ConfiguratorType &_conf;
  
   QuocMaterialConfigurator ( const ConfiguratorType &conf ) : _conf ( conf ) {}

    
// #ifdef APPROXCHARFCT_1DHOMOGENIZATION
//     RealType approxCharFct_vol ( const RealType v ) const { return 0.5 * ( v + 1.0 ); }
//     RealType approxCharFct_vol_Derivative ( const RealType /*v*/ ) const { return 0.5;}
//     RealType approxCharFct_vol_SecondDerivative ( const RealType /*v*/ ) const { return 0.0;}
// #endif

#ifdef QUOC_APPROXCHARFCT_WIRTH
    RealType approxCharFct_vol ( const RealType v ) const { return 0.25 * aol::Sqr( v + 1.0 ); }
    RealType approxCharFct_vol_Derivative ( const RealType v ) const { return 0.5 * ( v + 1.0 );}
    RealType approxCharFct_vol_SecondDerivative ( const RealType /*v*/ ) const { return 0.5;}  
    
    RealType approxCharFct_material ( const RealType v ) const { return 0.25 * aol::Sqr( v + 1.0 ); }
    RealType approxCharFct_material_Derivative ( const RealType v ) const { return 0.5 * ( v + 1.0 );}
    RealType approxCharFct_material_SecondDerivative ( const RealType /*v*/ ) const { return 0.5;}
#endif

#ifdef QUOC_APPROXCHARFCT_FOURTHORDER
    RealType approxCharFct_vol ( const RealType v ) const { return aol::Sqr( aol::Sqr( v + 1.0 ) ) / 16.; }
    RealType approxCharFct_vol_Derivative ( const RealType v ) const { return 0.25 * aol::Cub( v + 1.0 );}
    RealType approxCharFct_vol_SecondDerivative ( const RealType v ) const { return 0.75 * aol::Sqr( v + 1.0 );}  
    
    RealType approxCharFct_material ( const RealType v ) const { return aol::Sqr( aol::Sqr( v + 1.0 ) ) / 16.; }
    RealType approxCharFct_material_Derivative ( const RealType v ) const { return 0.25 * aol::Cub( v + 1.0 );}
    RealType approxCharFct_material_SecondDerivative ( const RealType v ) const { return 0.75 * aol::Sqr( v + 1.0 );}  
#endif


    RealType doubleWell ( const RealType v ) const { return 9.0 / 16.0 * aol::Sqr( v * v - 1.0 );}
    RealType doubleWellDerivative ( const RealType v ) const { return 9.0 / 4.0 * ( v * v - 1.0 ) * v;}
    RealType doubleWellSecondDerivative ( const RealType v ) const { return 9.0 / 4.0 * ( 3. * v * v - 1.0 );}
    
};


template < typename _ConfiguratorType >
class QuocMaterialOptimizationConfigurator : public QuocMaterialConfigurator<_ConfiguratorType> {
  
public :
  
  typedef _ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename DataTypeContainer::ParameterParserType ParameterParserType;
  
  const QuocMaterialOptInfo<RealType> _materialInfo;
  const RealType _factorVoidMaterial;
  RealType _factorComplianceCost, _factorVolumeCost, _factorInterfaceCost;
  const RealType _epsInterfaceLength;

  //used e.g. if one is only interested in optimal deformation
  QuocMaterialOptimizationConfigurator ( const ParameterParserType &parser, const ConfiguratorType &conf ) :
    QuocMaterialConfigurator<_ConfiguratorType> ( conf ),
    _materialInfo ( parser ),
    _factorVoidMaterial( parser.template get<double>( "MaterialOptimization.factorVoidMaterial" ) ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorVolumeCost ( parser.template get<double>( "MaterialOptimization.factorVolumeCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( 1. ){
        cout << endl << endl << "WARNING: this constructor of QuocMaterialOptimizationConfigurator should not be used for material optimization, since eps=1" << endl << endl;
    }
    

  QuocMaterialOptimizationConfigurator ( const ParameterParserType &parser, const ConfiguratorType &conf, const RealType epsInterfaceLength ) :
    QuocMaterialConfigurator<_ConfiguratorType> ( conf ),
    _materialInfo ( parser ),
    _factorVoidMaterial( parser.template get<double>( "MaterialOptimization.factorVoidMaterial" ) ),
    _factorComplianceCost ( parser.template get<double>( "MaterialOptimization.factorComplianceCost" ) ),
    _factorVolumeCost ( parser.template get<double>( "MaterialOptimization.factorVolumeCost" ) ),
    _factorInterfaceCost ( parser.template get<double> ( "MaterialOptimization.factorInterfaceCost" ) ),
    _epsInterfaceLength ( epsInterfaceLength ) {}

};


#endif
