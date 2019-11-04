#ifndef __QUOCHANDLER_H
#define __QUOCHANDLER_H

#include <loadAndSave.h>
#include <meshWithData.h>

#include "quocDefines.h"
#include "quocDiscreteFunction.h"

namespace quocFE {


template< typename ConfiguratorType >
void generateDirichletBoundaryMaskUponBoundaryType ( const QuocBoundaryType boundaryType, const typename ConfiguratorType::InitType &mesh,
                                                     typename ConfiguratorType::MaskType & mask ) {
  switch ( boundaryType ){
      
      case NOBOUNDARY : {
      } break;
      
      case LEFT : {
            mask = mesh._boundaryLeft;
      } break;
  
      default :
        throw std::invalid_argument( aol::strprintf ( "Wrong boundary condition. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        break;
    }
}



template< typename ConfiguratorType >
void generatePeriodicBoundaryMaskUponBoundaryType ( const QuocBoundaryType boundaryType, const typename ConfiguratorType::InitType &mesh,
                                                    typename ConfiguratorType::MaskType & mask, std::vector<int> & periodicIndices ) {

    periodicIndices = mesh._periodicIdentificationIndices;
  
  switch ( boundaryType ){
      
      case NOBOUNDARY : {
      } break;
      
      case ALL : {
            mask = mesh._boundaryPeriodic;
      } break;
  
      default :
        throw std::invalid_argument( aol::strprintf ( "Wrong boundary condition. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
        break;
    }
    
}


template< typename ConfiguratorType >
class QuocHandler{
  
public:
  
  typedef typename ConfiguratorType::RealType       RealType;
  typedef typename ConfiguratorType::InitType       MeshType;
  typedef typename ConfiguratorType::MaskType       MaskType;
  typedef typename ConfiguratorType::PointType      PointType;
  typedef typename ConfiguratorType::VectorType     VectorType;
  typedef typename ConfiguratorType::DTContainer    DataTypeContainer;
  typedef typename DataTypeContainer::ParameterParserType ParameterParserType;
  
  const ParameterParserType &_parser;
  const ConfiguratorType &_conf;
  const MeshType &_mesh;
  const int _numVertices, _numGlobalDofs;
  MaskType _DirichletMask;
  MaskType _PeriodicMask;
  std::vector<int> _PeriodicIndices;
  
public:
  
  QuocHandler( const ParameterParserType &Parser, const ConfiguratorType &conf ) : 
  _parser ( Parser),
  _conf ( conf ),
  _mesh( conf.getInitializer() ),
  _numVertices( conf.getNumGlobalDofs() ), _numGlobalDofs ( conf.getNumGlobalDofs() ){
    generateDirichletBoundaryMask( _DirichletMask );
    generatePeriodicBoundaryMask( _PeriodicMask, _PeriodicIndices );
  }
  
  void generateDirichletBoundaryMask ( MaskType & mask ) const{
    mask.resize( _numGlobalDofs, false );
    generateDirichletBoundaryMaskUponBoundaryType<ConfiguratorType> ( static_cast<QuocBoundaryType>( _parser.template get<int>( "InputMesh.DirichletBoundaryType" ) ), _mesh, mask );
  }
  
  void generatePeriodicBoundaryMask ( MaskType & mask, std::vector<int> &periodicIndices ) const{
    mask.resize( _numGlobalDofs, false );
    periodicIndices.resize( _numGlobalDofs );
    generatePeriodicBoundaryMaskUponBoundaryType<ConfiguratorType> ( static_cast<QuocBoundaryType>( _parser.template get<int>( "InputMesh.PeriodicBoundaryType" ) ), _mesh, mask, periodicIndices );
  }
  
  //! \note the following operations only delete the periodic value (i.e. set it to zero) 
  void collabseVectorPeriodically ( VectorType & vec ) const{
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
        if( this->getPeriodicMask()[nodeIdx] ){
            vec[nodeIdx] = 0.0;
        }
    }
  }
  
  void extendVectorPeriodically ( VectorType & vec ) const{
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ) vec[nodeIdx] = vec[this->getPeriodicIndices()[nodeIdx]];
    }
  }
  
  //! \note the following operations is additive: it adds the periodic values to the corresponding entry 
  void collabseVectorPeriodicallyAdditive ( VectorType & vec ) const{
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
        if( this->getPeriodicMask()[nodeIdx] ){
            vec[this->getPeriodicIndices()[nodeIdx]] += vec[nodeIdx];
            vec[nodeIdx] = 0.0;
        }
    }
  }
  
  void collabseMultiVectorPeriodically ( VectorType & dispPeriodic ) const{
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp) 
              dispPeriodic[nodeIdx + comp * _mesh.getNumVertices()] = 0.0;
      }
    }
  }
  
  void collabseMultiVectorPeriodicallyAdditive ( VectorType & dispPeriodic ) const{
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp){
              dispPeriodic[this->getPeriodicIndices()[nodeIdx] + comp * _mesh.getNumVertices()] += dispPeriodic[nodeIdx + comp * _mesh.getNumVertices()];
              dispPeriodic[nodeIdx + comp * _mesh.getNumVertices()] = 0.0;
          }
      }
    }
  }
  
  void extendMultiVectorPeriodically ( const VectorType & dispPeriodic, VectorType & dispPeriodicallyExtended ) const{
    dispPeriodicallyExtended = dispPeriodic;
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp) 
              dispPeriodicallyExtended[nodeIdx + comp * _mesh.getNumVertices()] = dispPeriodic[this->getPeriodicIndices()[nodeIdx] + comp * _mesh.getNumVertices()];
      }
    }
  }
  
  void extendMultiVectorPeriodically ( VectorType & dispPeriodicallyExtended ) const{
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp) 
              dispPeriodicallyExtended[nodeIdx + comp * _mesh.getNumVertices()] = dispPeriodicallyExtended[this->getPeriodicIndices()[nodeIdx] + comp * _mesh.getNumVertices()];
      }
    }
  }
  
  const MaskType & getDirichletMask ( ) const { return _DirichletMask;}
  const MaskType & getPeriodicMask ( ) const { return _PeriodicMask;}
  const std::vector<int> & getPeriodicIndices ( ) const { return _PeriodicIndices;}
  
 //==========================================================================================================================
 //==================================   Material     ========================================================================
 //==========================================================================================================================
  void constructConstantMaterial( VectorType &material, RealType materialConstant ) const{
    for( int i=0; i<_numVertices; ++i ) material[i] = materialConstant;
  }
  
  void constructLayeredMaterial( VectorType &material, const RealType startLayer, const RealType endLayer, const int direction = 0, 
                                const bool LayerHard = true, const RealType mhard = 1.0, const RealType msoft = -1.0 ) const{
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( (coords[direction] <= endLayer) && (coords[direction] >= startLayer) )  material[nodeIdx] = mhard;
        else material[nodeIdx] = msoft;
    }
    if( LayerHard == false ) material *= -1.0;        
  }
  
  void constructHoles( VectorType &material ) const{ 
     
    RealType lx = 0.0, ly = 0.0;
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        RealType tmp = cos( 6. * pi * coords[0] / lx ) * cos( 4. * pi * coords[1] ) + 0.6;
        tmp -= std::max<RealType>( 200.*(0.01 - coords[0] * coords[0] - (coords[1] - 0.5*ly) * (coords[1] - 0.5*ly) ) , 0. );
        tmp -= std::max<RealType>( 100. * ( coords[0] + coords[1] - lx - ly + 0.1), 0. );
        tmp -= std::max<RealType>( 100. * ( coords[0] - coords[1] - lx + 0.1 ), 0. );
       
        if( tmp > 1.0 ) tmp = 1.0;
        if( tmp < -1.0 ) tmp = -1.0;
        material[nodeIdx] = tmp;
    }
  }
  
  
  void constructSchwarzPSurface( VectorType &material ) const{ 
    RealType lx = 0.0, ly = 0.0, lz = 0.0;
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
        if( coords[2] > lz ) lz = coords[2];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        RealType tmp = cos( 2. * pi * coords[0] / lx ) + cos( 2. * pi *  coords[1] / ly ) + cos( 2. * pi *  coords[2] / lz );
        if( tmp > 0.0 ) material[nodeIdx] = 1.0;
        if( tmp < 0.0 ) material[nodeIdx] = -1.0;
        if( tmp == 0.0 ) material[nodeIdx] = 0.0;
    }
  }
  
  void constructGyroid( VectorType &material ) const{ 
    RealType lx = 0.0, ly = 0.0, lz = 0.0;
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        if( coords[0] > lx ) lx = coords[0];
        if( coords[1] > ly ) ly = coords[1];
        if( coords[2] > lz ) lz = coords[2];
    }
    const RealType pi = 4 * atan ( 1.0 );
    for ( int nodeIdx=0; nodeIdx < _numVertices; ++nodeIdx ) {
        const PointType& coords ( _mesh.getVertex(nodeIdx) );
        RealType tmp = sin( 2. * pi * coords[0] / lx ) * cos( 2. * pi * coords[1] / ly ) 
                     + sin( 2. * pi * coords[1] / ly ) * cos( 2. * pi * coords[2] / lz ) 
                     + sin( 2. * pi * coords[2] / lz ) * cos( 2. * pi * coords[0] / lx );
        if( tmp > 0.0 ) material[nodeIdx] = 1.0;
        if( tmp < 0.0 ) material[nodeIdx] = -1.0;
        if( tmp == 0.0 ) material[nodeIdx] = 0.0;
    }
  }
  
  void constructRandomMaterial( VectorType &material ) const{
      material = VectorType::Random( material.size() );
  }
  
  void switchMaterialType( VectorType &material ) const{
      switch( _parser.template get<int>( "Material.initMaterialType" ) ){
        case -1:
          aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s", _parser.template get<string> ( "Material.materialFile" ).c_str () ) ); 
          material *= -1.;
          break;
        case 0:
          cout << endl << endl << "load vector from file " << _parser.template get<string> ( "Material.materialFile" ).c_str() << endl << endl;
          aol::loadVectorFromFile<VectorType>( material, aol::strprintf ( "%s", _parser.template get<string> ( "Material.materialFile" ).c_str() ) ); 
          break;
        case 1:
          constructConstantMaterial( material, _parser.template get<double>( "Material.materialConstant" ) );
          break;
        case 2:
          constructLayeredMaterial( material, _parser.template get<double>( "Material.materialStartLayer" ), _parser.template get<double>( "Material.materialEndLayer" ), _parser.template get<double>( "Material.materialDirectionLayer" ), _parser.template get<bool>( "Material.materialLayerHard" ) );
          break;
        case 10:
            constructHoles( material );
            break;
        case 20:
            constructSchwarzPSurface( material );
            break;
        case 21:
            constructSchwarzPSurface( material );
            material *= -1.;
            break;
        case 22:
            constructGyroid( material );
            break;
        case 23:
            constructGyroid( material );
            material *= -1.;
            break;
        case 1000 :
            constructRandomMaterial( material );
            break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
  }
  
  void switchMaterialTypeForFixedVolume( const int designType, const RealType volHardMaterial, VectorType &material, string &designTypeName ) const{
      switch(designType){
        case 0 :
           designTypeName = "SchwarzPSurface";
           constructSchwarzPSurface( material );
           break;
        case 1 :
           designTypeName = "SchwarzPSurfaceSwitch";
           constructSchwarzPSurface( material );
           material *= -1.;
           break;
        case 2 :
           designTypeName = "Gyroid";
           constructGyroid( material );
           break;
        case 3 :
           designTypeName = "GyroidSwitch";
           constructGyroid( material );
           material *= -1.;
           break;
        default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
          break;
      }
  }
  
  
 //==========================================================================================================================
 //==================================   Plot     ============================================================================
 //==========================================================================================================================
 
 void plotUndeformedWithMaterial ( const VectorType & material, const string outfile_base_name ) const{
    MeshWithData<MeshType> meshSaver ( _mesh );
    meshSaver.addScalarData ( material, "material", VERTEX_DATA );
    meshSaver.saveAsVTKSTRUCTUREDPOINTS( aol::strprintf ( "%s/%s.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
 }
 
 void plot ( const VectorType & disp, const string outfile_base_name ) const{
    MeshType mesh ( _mesh );
    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        PointType coords = _mesh.getVertex( i );
        for( int comp = 0; comp < _conf.dimDomain; ++comp ) coords[comp] += disp[i + comp * _numGlobalDofs];
        mesh.setVertex( i, coords );
    }

    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
  }

  void plotWithMaterial ( const VectorType & disp, const VectorType & material, const string outfile_base_name ) const{
    MeshType mesh ( _mesh );
    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        PointType coords = _mesh.getVertex( i );
        for( int comp = 0; comp < _conf.dimDomain; ++comp ) coords[comp] += disp[i + comp * _numGlobalDofs];
        mesh.setVertex( i, coords );
    }
    
    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.addScalarData ( material, "material", VERTEX_DATA );
    meshSaver.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
  }
  
  
  void plotPeriodicDisp ( const VectorType & dispPeriodic, const string outfile_base_name ) const{
    
    VectorType disp ( dispPeriodic );
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp) disp[nodeIdx + comp * _mesh.getNumVertices()] = dispPeriodic[this->getPeriodicIndices()[nodeIdx] + comp * _mesh.getNumVertices()];
      }
    }
      
    MeshType mesh ( _mesh );
    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        PointType coords = _mesh.getVertex( i );
        for( int comp = 0; comp < _conf.dimDomain; ++comp ) coords[comp] += disp[i + comp * _numGlobalDofs];
        mesh.setVertex( i, coords );
    }
    
    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
  }
  
  
  void plotPeriodicDispWithMaterial ( const VectorType & dispPeriodic, const VectorType & material, const string outfile_base_name ) const{
    
    VectorType disp ( dispPeriodic );
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp) disp[nodeIdx + comp * _mesh.getNumVertices()] = dispPeriodic[this->getPeriodicIndices()[nodeIdx] + comp * _mesh.getNumVertices()];
      }
    }
      
    MeshType mesh ( _mesh );
    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        PointType coords = _mesh.getVertex( i );
        for( int comp = 0; comp < _conf.dimDomain; ++comp ) coords[comp] += disp[i + comp * _numGlobalDofs];
        mesh.setVertex( i, coords );
    }
    
    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.addScalarData ( material, "material", VERTEX_DATA );
    meshSaver.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
  }
  
  void plotPeriodicPlusAffineSymGradDispWithMaterial ( const VectorType & dispPeriodic, const VectorType & dispAffine, const VectorType & material, const string outfile_base_name ) const{
    
    VectorType disp ( dispPeriodic );
    for ( int nodeIdx=0; nodeIdx<_mesh.getNumVertices(); nodeIdx++ ) {
      if( this->getPeriodicMask()[nodeIdx] ){
          for( int comp=0; comp<_conf.dimDomain; ++comp) disp[nodeIdx + comp * _mesh.getNumVertices()] = dispPeriodic[this->getPeriodicIndices()[nodeIdx] + comp * _mesh.getNumVertices()];
      }
    }
      
    typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
    DerivativeVectorValuedType symMat; 
    QuocDiscreteFunctionDefaultAffine<DataTypeContainer, ConfiguratorType::dimDomain> dispAffineDFD ( dispAffine );
    symMat = dispAffineDFD.getSymGrad();
    
    MeshType mesh ( _mesh );
    MeshType meshPeriodicDeformed ( _mesh );
    MeshType meshAffineDeformed ( _mesh );
    for( int i = 0; i < mesh.getNumVertices(); ++i ){
        PointType coords = _mesh.getVertex( i );
        PointType offsetAffine = symMat * coords;
        PointType coordsAffine = coords + offsetAffine;
        meshAffineDeformed.setVertex( i, coordsAffine );
        PointType coordsPeriodic = coords;
        for( int comp = 0; comp < _conf.dimDomain; ++comp ) coordsPeriodic[comp] += disp[i + comp * _numGlobalDofs];
        meshPeriodicDeformed.setVertex(i, coordsPeriodic );
        mesh.setVertex( i, coordsPeriodic + offsetAffine );
    }
    
    MeshWithData<MeshType> meshSaverPeriodic ( meshPeriodicDeformed );
    meshSaverPeriodic.addScalarData ( material, "material", VERTEX_DATA );
    meshSaverPeriodic.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s_Periodic.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
    
    MeshWithData<MeshType> meshSaverAffine ( meshAffineDeformed );
    meshSaverAffine.addScalarData ( material, "material", VERTEX_DATA );
    meshSaverAffine.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s_Affine.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
    
    MeshWithData<MeshType> meshSaver ( mesh );
    meshSaver.addScalarData ( material, "material", VERTEX_DATA );
    meshSaver.saveAsVTKUNSTRUCTUREDGRID( aol::strprintf ( "%s/%s_Total.vtk", _parser.template get<std::string> ( "saving.saveDirectory" ).c_str (), outfile_base_name.c_str() ) );
  }
  
};

} //end namespace

#endif //__QUOCHANDLER_H
