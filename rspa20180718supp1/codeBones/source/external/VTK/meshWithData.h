#ifndef __MESHWITHDATA_H
#define __MESHWITHDATA_H

#include <general.h>

//! data can either belong to vertices or to faces
enum DataSupp { VERTEX_DATA, FACE_DATA };
//! vector-valued data can be saved as vectors, normals or texture coordinates (the file format also supports color scalars and lookup tables, which are not use here).
enum VectorSpec { VECTORS, NORMALS };



//! Use this class like:  MeshWithData<> ( mesh ) -> .addData ( result, "color", VERTEX_DATA ) ->  .saveAsVTK ( "result.vtk" );
//! This class is a container for a mesh plus data vectors on vertices and faces.  We will not store any of the data vectors here, but only keep pointers to them.
template <class MeshType>
class MeshWithData {
public:
  typedef typename MeshType::RealType RealType;
  typedef typename MeshType::IndicesOfElementType IndicesOfElementType;
  typedef typename MeshType::PointType PointType;
  typedef typename MeshType::VectorType VectorType;

  MeshWithData ( const MeshType & mesh ) : _mesh ( mesh ), _precision( 8 ) {}
  MeshWithData ( const MeshType & mesh, int precision ) : _mesh ( mesh ), _precision( precision ) {}

  MeshWithData & addScalarData ( const VectorType & data, string dataDescr, DataSupp supp ) {
    
    ScalarData entry = { dataDescr, &data };

    switch ( supp ) {
    case VERTEX_DATA:
      if ( data.size() != static_cast<unsigned>( _mesh.getNumVertices () ) ) throw std::invalid_argument( aol::strprintf ( "Wrong size. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      _scalarVertexData.push_back ( entry );
      break;

    case FACE_DATA:
      if ( data.size() != static_cast<unsigned>( _mesh.getNumElements () ) ) throw std::invalid_argument( aol::strprintf ( "Wrong size. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
      _scalarFaceData.push_back ( entry );
      break;
      
    default:
      throw std::invalid_argument( aol::strprintf ( "Unknown DataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }

  MeshWithData & addVectorData ( const VectorType & data, const int numComponents, string dataDescr, DataSupp supp, VectorSpec vSpec = VECTORS ) {

    VectorData entry = { dataDescr, vSpec, &data, numComponents };

    switch ( supp ) {
    case VERTEX_DATA:
      _vectorVertexData.push_back ( entry );
      break;
    case FACE_DATA:
      _vectorFaceData.push_back ( entry );
      break;
    default:
      throw std::invalid_argument( aol::strprintf ( "Unknown DataSupp. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    }
    return *this;
  }

  //! set precision in saving methods
  void setPrecisionTo( int prec ){ _precision = prec; }

  
  void saveAsVTKPOLYDATA ( string filename, const PointType &offsetPoints ) const {

    const bool polydata = true;
      
    const int numVertices = _mesh.getNumVertices();
       
    std::ofstream out ( filename.c_str() );  
    out << "# vtk DataFile Version 3.0" << endl
        << "written by method MeshWithData::saveAsVTK" << endl
        << "ASCII" << endl
        << "DATASET POLYDATA" << endl;
        
    // vertex coordinates
    PointType tmpPoint; const int dimPoint = tmpPoint.size();
    out << "POINTS " << numVertices << " float" << endl;
    for ( int nodeIter = 0; nodeIter < _mesh.getNumVertices(); ++nodeIter ) {
      const PointType& coords ( _mesh.getVertex(nodeIter) );
      for ( short i = 0; i < dimPoint; ++i ) out << ( i == 0 ? "" : " " ) << coords[i] - offsetPoints[i];
      if( dimPoint == 1 ) out << " 0 0";
      if( dimPoint == 2 ) out << " 0";
      out << endl;
    }

    IndicesOfElementType tmpElementIndices; const int dimPolygon = tmpElementIndices.size();
    out << "POLYGONS " << _mesh.getNumElements() << " " << (1 + dimPolygon ) * _mesh.getNumElements() << endl;
    // vertex indices of element
    for ( int elementIter = 0; elementIter < _mesh.getNumElements(); ++elementIter ) {
      out << dimPolygon << " ";
      for ( short i = 0; i < dimPolygon; ++i )
        out << ( i == 0 ? "" : " " ) << _mesh.getElementNodeIdx( elementIter, i );
      out << endl;
    }
    
//     out << "CELL_TYPES " << _mesh.getNumElements() << endl;
//     // VTK CELL TYPES
//     for ( int elementIter = 0; elementIter < _mesh.getNumElements(); ++elementIter ) {
//       out << _mesh._VTKCELLTYPE << endl;
//     }

    this->writePointData( out );
    this->writeCellData( out, polydata );
  }
  
  void saveAsVTKPOLYDATA ( const string filename ) const {
      PointType zeroOffset; zeroOffset.setZero();
      this->saveAsVTKPOLYDATA( filename,  zeroOffset );
  }
  
  
  
  
  void saveAsVTKUNSTRUCTUREDGRID ( string filename, const PointType &offsetPoints ) const {

    const int numVertices = _mesh.getNumVertices();
       
    std::ofstream out ( filename.c_str() );  
    out << "# vtk DataFile Version 3.0" << endl
        << "written by method MeshWithData::saveAsVTK" << endl
        << "ASCII" << endl
        << "DATASET UNSTRUCTURED_GRID" << endl;
        
    // vertex coordinates
    PointType tmpPoint; const int dimPoint = tmpPoint.size();
    out << "POINTS " << numVertices << " float" << endl;
    for ( int nodeIter = 0; nodeIter < _mesh.getNumVertices(); ++nodeIter ) {
      const PointType& coords ( _mesh.getVertex(nodeIter) );
      for ( short i = 0; i < dimPoint; ++i ) out << ( i == 0 ? "" : " " ) << coords[i] - offsetPoints[i];
      if( dimPoint == 1 ) out << " 0 0";
      if( dimPoint == 2 ) out << " 0";
      out << endl;
    }

    IndicesOfElementType tmpElementIndices; const int dimPolygon = tmpElementIndices.size();
//     out << "POLYGONS " << _mesh.getNumElements() << " " << (1 + dimPolygon ) * _mesh.getNumElements() << endl;
    out << "CELLS " << _mesh.getNumElements() << " " << (1 + dimPolygon ) * _mesh.getNumElements() << endl;
    // vertex indices of element
    for ( int elementIter = 0; elementIter < _mesh.getNumElements(); ++elementIter ) {
      out << dimPolygon << " ";
      for ( short i = 0; i < dimPolygon; ++i )
        out << ( i == 0 ? "" : " " ) << _mesh.getElementNodeIdx( elementIter, i );
      out << endl;
    }
    
    out << "CELL_TYPES " << _mesh.getNumElements() << endl;
    // VTK CELL TYPES
    for ( int elementIter = 0; elementIter < _mesh.getNumElements(); ++elementIter ) {
      out << _mesh._VTKCELLTYPE << endl;
    }

    this->writePointData( out );
    this->writeCellData( out );
  }
  
  void saveAsVTKUNSTRUCTUREDGRID ( const string filename ) const {
      PointType zeroOffset; zeroOffset.setZero();
      this->saveAsVTKUNSTRUCTUREDGRID( filename,  zeroOffset );
  }
  
  
  
  
  void saveAsVTKSTRUCTUREDPOINTS ( const string filename, const PointType &offsetPoints ) const {

    std::ofstream out ( filename.c_str() );  
    out << "# vtk DataFile Version 3.0" << endl
        << "written by method MeshWithData::saveAsVTKSTRUCTUREDPOINTS" << endl
        << "ASCII" << endl
        << "DATASET STRUCTURED_POINTS" << endl;
        
      
    out << "DIMENSIONS";
    for( int i=0; i<3; ++i ) out << " " << _mesh.getNumDofs(i); 
    out << endl;
        
    out << "SPACING";
    for( int i=0; i<3; ++i ) out << " " << _mesh.getMeshSize(i);
    out << endl;
    
    out << "ORIGIN";
    PointType tmpPoint; const int dimPoint = tmpPoint.size();
    for(int i=0; i<dimPoint; ++i) out << " " << -1. * offsetPoints[i];
    if( dimPoint == 2 ) out << " 0";
    out << endl;
    
    this->writePointData( out );
  }
  
  void saveAsVTKSTRUCTUREDPOINTS ( const string filename ) const {
      PointType zeroOffset; zeroOffset.setZero();
      this->saveAsVTKSTRUCTUREDPOINTS( filename,  zeroOffset );
  }
  

protected:

  struct ScalarData {
    string             _descr;
    const VectorType * _data;
  };

  struct VectorData {
    string             _descr;
    VectorSpec         _spec;
    const VectorType * _data;
    const int          _numComponents;
  };

  std::vector<ScalarData> _scalarVertexData;
  std::vector<VectorData> _vectorVertexData;
  std::vector<ScalarData> _scalarFaceData;
  std::vector<VectorData> _vectorFaceData;

  const MeshType _mesh;
  int _precision;      
  
  
  void writePointData( std::ofstream &out ) const {
   
    const int numVertices = _mesh.getNumVertices();  
      
    if ( _scalarVertexData.size() > 0 || _vectorVertexData.size() > 0 ) out << "POINT_DATA " << numVertices << endl;
    
    // scalar data on vertices
    for (int i = 0; i < _scalarVertexData.size(); ++i) {
      out << "SCALARS " <<  _scalarVertexData[i]._descr << " float" << endl;
      out << "LOOKUP_TABLE default" << endl;
      for ( unsigned vx = 0; vx < (*_scalarVertexData[i]._data).size(); ++vx ) out << (*_scalarVertexData[i]._data)[vx] << endl;
    }
    
    // vector data on vertices
    for (int i = 0; i < _vectorVertexData.size(); ++i) {
      string spec;
      if ( _vectorVertexData[i]._spec == VECTORS ) spec = "VECTORS";
      if ( _vectorVertexData[i]._spec == NORMALS ) spec = "NORMALS";
      out << spec << " " << _vectorVertexData[i]._descr << " float" << endl;
      for ( int vx = 0; vx < numVertices; ++vx ) {
        for ( int comp = 0; comp < _vectorVertexData[i]._numComponents; ++comp )
          out << ( comp == 0 ? "" : " " ) << (*_vectorVertexData[i]._data)(vx + comp * numVertices); //TODO index mapper
        out << endl;
      }
    }
  }
  
  
  //as in writePointData: if ( _scalarFaceData.size() > 0 || _vectorFaceData.size() > 0 ) out << "CELL_DATA" << _mesh.getNumElements() << endl;
  void writeCellData( std::ofstream &out, const bool polydata = false ) const {
   
    // scalar data on faces
    for( int i=0; i < _scalarFaceData.size(); ++i ) {
        out << "CELL_DATA "; 
        out << _mesh.getNumElements() << endl;
        out << "SCALARS " << _scalarFaceData[i]._descr << " float" << endl;
        out << "LOOKUP_TABLE default" << endl;
        for ( int vx = 0; vx < (*_scalarFaceData[i]._data).size(); ++vx ) out << (*_scalarFaceData[i]._data)[vx] << endl;
    }
    
    // vector data on faces
    for (int i = 0; i < _vectorFaceData.size(); ++i) {
      string spec;
      if ( _vectorFaceData[i]._spec == VECTORS ) spec = "VECTORS";
      if ( _vectorFaceData[i]._spec == NORMALS ) spec = "NORMALS";
      out << spec << " " << _vectorFaceData[i]._descr << " float" << endl;
      //TODO
    }
      
  }
  
  
};















//! Use this class like:  PointCloud<> ( mesh ) -> .addData ( result, "color", VERTEX_DATA ) ->  .saveAsVTK ( "result.vtk" );
//! This class is a container for a mesh plus data vectors on vertices and faces.  We will not store any of the data vectors here, but only keep pointers to them.
template <typename DataTypeContainer>
class PointCloud {
public:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef typename DataTypeContainer::TangentVecType TangentVecType;

  PointCloud ( ) : _precision( 8 ) {}
  
  //! set precision in saving methods
  void setPrecisionTo( int prec ){ _precision = prec; }

  
  void saveAsLegacyVTK ( const string filename, std::vector<PointType> &pointsVec ) const {
       
    std::ofstream out ( filename.c_str() );  
    out << "# vtk DataFile Version 3.0" << endl
        << "written by method PointCloud::saveAsLegacyVTK" << endl
        << "ASCII" << endl
        << "DATASET UNSTRUCTURED_GRID" << endl;
        

    // vertex coordinates
    PointType tmpPoint; const int dimPoint = tmpPoint.size();
    out << "POINTS " << pointsVec.size() << " float" << endl;
    for ( int nodeIter = 0; nodeIter < pointsVec.size(); ++nodeIter ) {
      const PointType coords ( pointsVec[nodeIter] );
      for ( short i = 0; i < dimPoint; ++i )
        out << ( i == 0 ? "" : " " ) << coords[i];
      if( dimPoint == 1 ) out << " 0 0";
      if( dimPoint == 2 ) out << " 0";
      out << endl;
    }

    const int dimPolygon = 1;
    out << "CELLS " << pointsVec.size() << " " << (2 ) * pointsVec.size() << endl;
    // vertex indices of element
    for ( int elementIter = 0; elementIter < pointsVec.size(); ++elementIter ) {
      out << dimPolygon << " " << elementIter << endl;
    }
    
    // VTK CELL TYPES (1 is Point)
    out << "CELL_TYPES " << pointsVec.size() << endl;
    for ( int elementIter = 0; elementIter < pointsVec.size(); ++elementIter ) out << 1 << endl;

  }
  
  void saveAsLegacyVTK ( const string filename, std::vector<PointType> &pointsVec, std::vector<RealType> &pointDataVec ) const {
       
    std::ofstream out ( filename.c_str() );  
    out << "# vtk DataFile Version 3.0" << endl
        << "written by method PointCloud::saveAsLegacyVTK" << endl
        << "ASCII" << endl
        << "DATASET UNSTRUCTURED_GRID" << endl;
        

    // vertex coordinates
    PointType tmpPoint; const int dimPoint = tmpPoint.size();
    out << "POINTS " << pointsVec.size() << " float" << endl;
    for ( int nodeIter = 0; nodeIter < pointsVec.size(); ++nodeIter ) {
      const PointType coords ( pointsVec[nodeIter] );
      for ( short i = 0; i < dimPoint; ++i )
        out << ( i == 0 ? "" : " " ) << coords[i];
      if( dimPoint == 1 ) out << " 0 0";
      if( dimPoint == 2 ) out << " 0";
      out << endl;
    }

    const int dimPolygon = 1;
    out << "CELLS " << pointsVec.size() << " " << (2 ) * pointsVec.size() << endl;
    // vertex indices of element
    for ( int elementIter = 0; elementIter < pointsVec.size(); ++elementIter ) {
      out << dimPolygon << " " << elementIter << endl;
    }
    
    out << "CELL_TYPES " << pointsVec.size() << endl;
    // VTK CELL TYPES
    for ( int elementIter = 0; elementIter < pointsVec.size(); ++elementIter ) {
      out << 1 << endl;
    }
    
    // scalar data on vertices //todo: nameOfScalarVec
    out << "POINT_DATA " << pointsVec.size() << endl;
    out << "SCALARS " << "nameOfScalarVec" << " float" << endl;
    out << "LOOKUP_TABLE default" << endl;
    for ( int nodeIter = 0; nodeIter < pointDataVec.size(); ++nodeIter ) {
        out << pointDataVec[nodeIter] << endl;
    }
  }
  
  
  void saveAsLegacyVTK ( const string filename, std::vector<PointType> &pointsVec, std::vector<TangentVecType> &pointDataVec ) const {
       
    std::ofstream out ( filename.c_str() );  
    out << "# vtk DataFile Version 3.0" << endl
        << "written by method PointCloud::saveAsLegacyVTK" << endl
        << "ASCII" << endl
        << "DATASET UNSTRUCTURED_GRID" << endl;
        

    // vertex coordinates
    PointType tmpPoint; const int dimPoint = tmpPoint.size();
    out << "POINTS " << pointsVec.size() << " float" << endl;
    for ( int nodeIter = 0; nodeIter < pointsVec.size(); ++nodeIter ) {
      const PointType coords ( pointsVec[nodeIter] );
      for ( short i = 0; i < dimPoint; ++i )
        out << ( i == 0 ? "" : " " ) << coords[i];
      if( dimPoint == 1 ) out << " 0 0";
      if( dimPoint == 2 ) out << " 0";
      out << endl;
    }

    const int dimPolygon = 1;
    out << "CELLS " << pointsVec.size() << " " << (2 ) * pointsVec.size() << endl;
    // vertex indices of element
    for ( int elementIter = 0; elementIter < pointsVec.size(); ++elementIter ) {
      out << dimPolygon << " " << elementIter << endl;
    }
    
    out << "CELL_TYPES " << pointsVec.size() << endl;
    // VTK CELL TYPES
    for ( int elementIter = 0; elementIter < pointsVec.size(); ++elementIter ) {
      out << 1 << endl;
    }

    // scalar data on vertices
    out << "POINT_DATA " << pointsVec.size() << endl;
    out << "VECTORS " << "nameOfVec" << " float" << endl;
    for ( int nodeIter = 0; nodeIter < pointDataVec.size(); ++nodeIter ) {
        for ( int comp = 0; comp < 3; ++comp )
          out << ( comp == 0 ? "" : " " ) << pointDataVec[nodeIter](comp);
        out << endl;
     }
  }
  
  protected:

  int _precision;      
};

  
 

#endif
