#ifndef __MATERIALOPTIMIZATIONDEFINES__H
#define __MATERIALOPTIMIZATIONDEFINES__H

#include <unordered_map>

template<typename RealType>
class Material {
protected:
  const std::string _name;
  RealType _density;       // in kg/m^3
  RealType _elastModulus;  // in GPa
  RealType _poissonRatio;  // ( nu )
  RealType _lambda;        // 1st lame constant in GPa = N/m^2
  RealType _mu;            // 2nd lame constant (shear modulus) in GPa = N/m^2

public:
  Material() :
    _name ( "DefaultMaterial" ),
    _density ( 0 ),_elastModulus ( 0 ), _poissonRatio ( 0 ), _lambda ( 0 ), _mu ( 0 ) {}

  Material ( std::string Name,
             RealType Density,
             RealType ElastModulus, RealType PoissonRatio ) :
    _name ( Name ),
    _density ( Density ),
    _elastModulus ( ElastModulus ),
    _poissonRatio ( PoissonRatio )
     { 
        this->setLambda( _elastModulus, _poissonRatio );
        this->setMu( _elastModulus, _poissonRatio );
     }
    
   Material ( const Material<RealType>& material ) :
    _name ( material.getName() ),
    _density ( material.getDensity() ),
    _elastModulus ( material.getElastModulus() ),
    _poissonRatio ( material.getPoissonRatio() ),
    _lambda ( material.getLambda() ),
    _mu ( material.getMu() ){ }
    
  ~Material() {}

  // ************** inspectors ************** //
  const std::string getName() const { return _name;}
  RealType getDensity() const {return _density;}
  RealType getElastModulus() const {return _elastModulus;}
  RealType getPoissonRatio() const {return _poissonRatio;}
  RealType getLambda() const {return _lambda;}
  RealType getMu() const {return _mu;}

  void print() const {
    std::cout << "-------- " << getName() << " --------" << std::endl
              << "density              : " << getDensity() << " kg/m^3" << std::endl
              << "elast modulus   (E)  : " << getElastModulus() << " GPa" << std::endl
              << "poisson ratio   (nu) : " << getPoissonRatio()  << std::endl
              << "1st Lame    (lambda) : " << getLambda() << " GPa" << std::endl
              << "2nd Lame        (mu) : " << getMu() << " GPa" << std::endl;
  }

  void setDensity ( RealType value_KgPerM3 ){ _density = value_KgPerM3;}
  void setElastModulus ( RealType value_GPa ) {_elastModulus = value_GPa;}
  void setPoissonRatio ( RealType value ){ _poissonRatio = value;}
  void setLambda( const RealType ElastModulus, const RealType PoissonRatio) { _lambda = ElastModulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );}
  void setMu( const RealType ElastModulus, const RealType PoissonRatio ) { _mu = ElastModulus / ( 2. * ( 1. + PoissonRatio ) ) ;}
  void setLambdaMu( const RealType lambda, const RealType mu ) {
    _lambda = lambda; _mu = mu;
  }
  
  void set( const Material<RealType>& material ) {
    setDensity( material._density ); setElastModulus( material._elastModulus ); setPoissonRatio( material._poissonRatio );
    setLambda( _elastModulus, _poissonRatio ); setMu( _elastModulus, _poissonRatio );
  }
};



template<typename RealType>
class MaterialDictionary {
protected :   
  std::unordered_map<std::string, Material<RealType> > _MaterialDictionary;
  
public :
    MaterialDictionary ( ) {
    //                                                    Name                  Density     E       nu   
    static const Material<RealType> DummyMaterial_Nu0   ( "DummyNu0"            , 1         ,   1     , 0.0    );
    static const Material<RealType> DummyMaterial_Nu025 ( "DummyNu025"          , 1         ,   1     , 0.25   );
    static const Material<RealType> DummyMaterial_E100  ( "DummyE100"           , 1         , 100     , 0.25   );
    static const Material<RealType> DummyMaterial_E10   ( "DummyE10"            , 1         ,  10     , 0.25  );
    static const Material<RealType> DummyMaterial_E5    ( "DummyE5"             , 1         ,   5     , 0.25  );
    static const Material<RealType> DummyMaterial_E5nu0125 ( "DummyE5nu0125"    , 1         ,   5     , 0.125  );
    static const Material<RealType> DummyMaterial_E5nu0375 ( "DummyE5nu0375"    , 1         ,   5     , 0.375  );
    static const Material<RealType> DummyMaterial_E1    ( "DummyE1"             , 1         ,   1     , 0.25  );
    static const Material<RealType> Air                 ( "Air"                 , 0.001225  ,   1.e-8 , 0.0    );
    static const Material<RealType> Aluminum            ( "Aluminum"            , 2.71e3    ,  70     , 0.33  );
    static const Material<RealType> Brass               ( "Brass"               , 8.5e3     , 103     , 0.33  );
    static const Material<RealType> Bone                ( "Bone"                , 1         , 10   ,    0.2   );
    static const Material<RealType> Carbon              ( "Carbon"              , 2.25e3    ,   6.9   , 0.0     );
    static const Material<RealType> Ceramic             ( "Ceramic"             , 2.5e3     , 350     , 0.0     );
    static const Material<RealType> Concrete            ( "Concrete"            , 2.35e3    ,  24     , 0.15  );
    static const Material<RealType> Copper              ( "Copper"              , 8.94e3    , 115     , 0.345 );
    static const Material<RealType> Glass               ( "Glass"               , 2.6e3     ,  65     , 0.235 );
    static const Material<RealType> Gold                ( "Gold"                , 19.32e3   ,  83     , 0.44  );
    static const Material<RealType> Iron                ( "Iron"                , 7.87e3    , 170     , 0.25  );
    static const Material<RealType> Magnesium           ( "Magnesium"           , 1.74e3    ,  41     , 0.35  );
    static const Material<RealType> Nickel              ( "Nickel"              , 8.89e3    , 210     , 0.31  );
    static const Material<RealType> Nitinol_Austenite   ( "Nitinol Austinite"   , 6.45      ,  83     , 0.33  );
    static const Material<RealType> Nitinol_Martensite  ( "Nitinol Martensite"  , 6.45      ,  24     , 0.33  );
    static const Material<RealType> Nylon               ( "Nylon"               , 1.1e3     ,   2.4   , 0.4   );
    static const Material<RealType> Platinum            ( "Platinum"            , 21.4e3    , 145     , 0.38  );
    static const Material<RealType> Polymer             ( "Polymer"             , 1.        , 0.25    , 0.3  );
    static const Material<RealType> Rubber              ( "Rubber"              , 1.15e3    ,   0.0023, 0.45 );
    static const Material<RealType> Silver              ( "Silver"              , 10.49e3   ,  76     , 0.37  );
    static const Material<RealType> Steel               ( "Steel"               , 7.85e3    , 200     , 0.285 );
    static const Material<RealType> Stone_Granite       ( "Stone Granite"       , 2.6e3     ,  55     , 0.25  );
    static const Material<RealType> Stone_Marble        ( "Stone Marble"        , 2.75e3    ,  75     , 0.25  );
    static const Material<RealType> Tin                 ( "Tin"                 , 7.3e3     ,  42     , 0.36  );
    static const Material<RealType> Titanium            ( "Titanium"            , 4.54e3    , 110     , 0.33  );
    static const Material<RealType> Zinc                ( "Zinc"                , 7.14e3    , 108     , 0.25  );
    
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyNu0"      , DummyMaterial_Nu0 ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyNu025"    , DummyMaterial_Nu025 ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyE100"     , DummyMaterial_E100 ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyE10"      , DummyMaterial_E10 ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyE5"       , DummyMaterial_E5 ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyE5nu0125"       , DummyMaterial_E5nu0125 ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyE5nu0375"       , DummyMaterial_E5nu0375 ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "DummyE1"       , DummyMaterial_E1 ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Air"                , Air ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Aluminum"           , Aluminum ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Brass"              , Brass ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Bone"               , Bone ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Carbon"             , Carbon ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Ceramic"            , Ceramic ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Concrete"           , Concrete ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Copper"             , Copper ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Glass"              , Glass ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Gold"               , Gold ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Iron"               , Iron ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Magnesium"          , Magnesium ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Nickel"             , Nickel ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Nitinol Austenite"  , Nitinol_Austenite ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Nitinol Martensite" , Nitinol_Martensite ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Nylon"              , Nylon ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Platinum"           , Platinum  ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Polymer"            , Polymer ) ); 
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Rubber"             , Rubber ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Silver"             , Silver  ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Steel"              , Steel ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Stone Granite"      , Stone_Granite ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Stone Marble"       , Stone_Marble ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Tin"                , Tin  ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Titanium"           , Titanium  ) );
    _MaterialDictionary.insert( std::pair< std::string, Material<RealType> >( "Zinc"               , Zinc ) );
        
    }
    
    Material<RealType> getMaterialFromDictionary( const std::string name ) const{
        typename std::unordered_map< std::string, Material<RealType> >::const_iterator iter;
        iter = _MaterialDictionary.find( name );
        if( iter == _MaterialDictionary.end() ) throw std::invalid_argument( aol::strprintf ( "Unknown Material. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
        return iter->second;
    }
    
};


#endif
