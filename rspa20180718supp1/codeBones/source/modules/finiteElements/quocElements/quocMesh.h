#ifndef __QUOCMESH_H
#define __QUOCMESH_H

#include <general.h>
#include "quocDefines.h"

namespace quocFE {

    
    
    
template<typename DataTypeContainer>
class QuocElement1D {

  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::DomVecType      DomVecType;
  typedef typename DataTypeContainer::PointType       PointType;
  typedef typename DataTypeContainer::Indices1DType   Indices1DType;
  typedef std::vector<PointType>                      VertexIterator;

  public :   
    // global indices of element an nodes
    int _globIdx;
    Indices1DType _globNodeIdx;
    PointType _nodes[2];
    
public:
  
  QuocElement1D () : _globIdx(-1) {}
  
  QuocElement1D( const int globalIdx, const Indices1DType globalNodeIndex, const VertexIterator &nodes  ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex ){
        for ( int i = 0; i < 2; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
    }
  
  ~QuocElement1D(){}
    
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;} //TODO old: globIdx, getIndex
  const Indices1DType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}

};



//Line (0) to (lx)
// DOFs in x: Nx
// meshsize: (hx) = (lx/(Nx-1))
template< typename DataTypeContainer, typename QuocElementType = QuocElement1D<DataTypeContainer> >
class QuocMesh1D {
public:
  typedef QuocElementType ElementType;
  typedef typename DataTypeContainer::RealType      RealType;
  typedef typename DataTypeContainer::PointType     PointType;
  typedef typename DataTypeContainer::Indices1DType Indices1DType;
  typedef Indices1DType                             IndicesOfElementType;
  typedef typename DataTypeContainer::VectorType    VectorType;
  typedef typename DataTypeContainer::MaskType      MaskType;
  typedef typename DataTypeContainer::IntVecType    IntVecType;
  typedef typename DataTypeContainer::DomVecType    DomVecType;
  
protected:
  const int _Nx;
  IntVecType _NumDofVec;
  const RealType _lx;
  PointType _LengthVec;
  const RealType _hx;
  std::vector< PointType > _vertexIterator;
  std::vector< ElementType > _elementIterator;
public :
  MaskType _boundaryLeft, _boundaryRight;
  MaskType _boundaryPeriodic;
  std::vector<int> _periodicIdentificationIndices;
  
public:
  //! Create empty QuocMesh
  QuocMesh1D ( ) : _vertexIterator(), _elementIterator() { }
  
  QuocMesh1D ( const int Nx, const RealType lx ) : 
    _Nx (Nx), _lx( lx ), _hx( _lx / static_cast<RealType>(_Nx - 1 ) )
  {
        _NumDofVec[0] = _Nx;
        _LengthVec = _lx;
      
        const int numGlobalDofs = _Nx;
        _boundaryLeft.resize( numGlobalDofs, false ); _boundaryRight.resize( numGlobalDofs, false );
        _boundaryPeriodic.resize( numGlobalDofs, false );
                
        for( int xIdx=0; xIdx<_Nx; ++xIdx ){
              PointType vertex ( _lx * xIdx / static_cast<RealType>( _Nx - 1 ) );
              int nodeIdx = this->pushBackVertex( vertex );
              if( xIdx == 0 ) _boundaryLeft[nodeIdx] = true;
              if( xIdx == _Nx - 1 ) _boundaryRight[nodeIdx] = true;
              if( xIdx == _Nx - 1 ) _boundaryPeriodic[nodeIdx] = true;
        }
          
      for( int nodeIdx = 0; nodeIdx < numGlobalDofs; ++nodeIdx ){
          if( (!_boundaryRight[nodeIdx]) ){
              //! \note ordering of elements has to coincide with ordering of base functions
              const int numNewElement = this->pushBackElement( Indices1DType(nodeIdx, nodeIdx + 1 ) );  
          }
      }
      
      //Periodic indices to identify facing presenting node
      _periodicIdentificationIndices.resize( numGlobalDofs );
      for( int xIdx=0; xIdx<_Nx; ++xIdx ){
            int periodicXIdx = xIdx;
            if( xIdx == _Nx - 1 ) periodicXIdx = 0;
            _periodicIdentificationIndices[this->getGlobalNodeIndex(xIdx)] = this->getGlobalNodeIndex(periodicXIdx);
      }
      
  }
    
  virtual ~QuocMesh1D ( ) {} 
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumElements ( ) const { return ( static_cast<int> ( _elementIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new element and return global index
  int pushBackElement ( const Indices1DType nodeIdx ) {
    int globalIdx = getNumElements();
    _elementIterator.push_back ( QuocElementType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }

  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  int getGlobalNodeIndex ( const int xIdx ) const {  return xIdx;}
  int getGlobalElementIndex ( const int xFac ) const { return xFac;} 
  
  const QuocElementType& getElement ( const int num ) const {return _elementIterator[num];}
  QuocElementType & getElement ( const int num ) {return _elementIterator[num];}
  void setElement ( const int num, const QuocElementType Arg ) { _elementIterator[num] = Arg;}

  int getElementNodeIdx  ( const int num, const int localNode ) const { return _elementIterator[num].getGlobalNodeIdx(localNode);}
  void setElementNodeIdx ( const int num, const int localNode, const int newIdx ) {_elementIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  
  const RealType getMeshSize( const int direction ) const{ 
    if( direction == 0 ) return _hx;
    if( direction == 1 ) return 0.;
    if( direction == 2 ) return 0.;
    throw std::invalid_argument( aol::strprintf ( "Wrong direction in getMeshSize in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  const RealType getInterfaceWith( ) const{ return _hx; }
  const RealType getWidth( const int direction ) const {
   switch( direction ){
        case 0 :    return _lx;    break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  const int getNumDofs( const int direction ) const {
   switch( direction ){
        case 0 : return _Nx; break;
        case 1 : return 1; break; //for compatibility with vtk
        case 2 : return 1; break; //for compatibility with vtk
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel in QuocMesh1D get num dofs. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  void getNumDofVec( IntVecType &vec ) const { vec = _NumDofVec;}
  void getLenghtVec( PointType &vec ) const { vec = _LengthVec;}
  
};
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 
    
    
template<typename DataTypeContainer>
class QuocElement2D {

  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::DomVecType      DomVecType;
  typedef typename DataTypeContainer::PointType       PointType;
  typedef typename DataTypeContainer::Indices2DType   Indices2DType;
  typedef std::vector<PointType>                      VertexIterator;

  public :   
    // global indices of element an nodes
    int _globIdx;
    Indices2DType _globNodeIdx;
    PointType _nodes[4];
    
public:
  
  QuocElement2D () : _globIdx(-1) {}
  
  QuocElement2D( const int globalIdx, const Indices2DType globalNodeIndex, const VertexIterator &nodes  ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex ){
        for ( int i = 0; i < 4; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
    }
  
  ~QuocElement2D(){}
    
  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const Indices2DType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}

};




template<typename DataTypeContainer>
class QuocBoundaryElement2D {

  typedef typename DataTypeContainer::RealType              RealType;
  typedef typename DataTypeContainer::DomVecType            DomVecType;
  typedef typename DataTypeContainer::DomVecTypeBoundary    DomVecTypeBoundary;
  typedef typename DataTypeContainer::PointType             PointType;
  typedef typename DataTypeContainer::Indices2DType         Indices2DType;
  typedef std::vector<PointType>                            VertexIterator;
  typedef QuocElement2D<DataTypeContainer>                  ElementType;

  protected :   
    // global indices of element an nodes
    int _globBoundaryIdx;
    const ElementType _element;
    const QuocBoundaryType _boundaryType;
    DomVecType _normal;
    PointType _boundaryNodes[2];
    PointType _boundaryNodesRefCoord[2];
    const RealType _h;
    int _nodeIndicesOf2DElement[2];
    
public:
  
  QuocBoundaryElement2D () : _globBoundaryIdx(-1) {}
  
  QuocBoundaryElement2D( const int globalBoundaryIdx, const ElementType & element, const QuocBoundaryType &boundaryType, const RealType h  ) : 
     _globBoundaryIdx(globalBoundaryIdx), _element (element), _boundaryType ( boundaryType ), _h(h) {
         
         _normal.setZero();
         switch( boundaryType ){
             case LEFT:{
               _boundaryNodesRefCoord[0] << 0.0, 0.0; _boundaryNodesRefCoord[1] << 0.0, 1.0;
               _boundaryNodes[0] = element.getNode(0);_boundaryNodes[1] = element.getNode(2);
               _normal[0] = -1.;
               _nodeIndicesOf2DElement[0] = 0; _nodeIndicesOf2DElement[1] = 2;
             }break;
                 
             case RIGHT:{
               _boundaryNodesRefCoord[0] << 1.0, 0.0; _boundaryNodesRefCoord[1] << 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(1);_boundaryNodes[1] = element.getNode(3);
               _normal[0] = 1.;
               _nodeIndicesOf2DElement[0] = 1; _nodeIndicesOf2DElement[1] = 3;
             }break;
                 
             case BOTTOM:{
               _boundaryNodesRefCoord[0] << 0.0, 0.0; _boundaryNodesRefCoord[1] << 1.0, 0.0;
               _boundaryNodes[0] = element.getNode(0); _boundaryNodes[1] = element.getNode(1);
               _normal[1] = -1.;
               _nodeIndicesOf2DElement[0] = 0; _nodeIndicesOf2DElement[1] = 1;
             }break;
                 
             case TOP:{
               _boundaryNodesRefCoord[0] << 0.0, 1.0; _boundaryNodesRefCoord[1] << 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(2); _boundaryNodes[1] = element.getNode(3);
               _normal[1] = 1.;
               _nodeIndicesOf2DElement[0] = 2; _nodeIndicesOf2DElement[1] = 3;
             }break;
                 
             default:
                 break;
         }
    }
  
  ~QuocBoundaryElement2D(){}

  const ElementType & getElement() const {return _element;}
  int getNodeIndexOfElement(const int index ) const{ return _nodeIndicesOf2DElement[index];}
  const DomVecType& getNormal() const {return _normal;}
  const PointType& getBoundaryNode ( int i ) const { return _boundaryNodes[i];}
  PointType& getBoundaryNode ( int i ) { return _boundaryNodes[i];}
  const PointType& getBoundaryNodeRefCoord ( int i ) const { return _boundaryNodesRefCoord[i];}
  PointType& getBoundaryNodeRefCoord ( int i ) { return _boundaryNodesRefCoord[i];}
  const RealType getH() const {return _h;}
  const QuocBoundaryType getBoundaryType() const {return _boundaryType;}
  void getRefCoord( const DomVecTypeBoundary &refCoordBoundary, DomVecType & refCoord ) const {
    refCoord = ( 1. - refCoordBoundary ) * this->getBoundaryNodeRefCoord(0) + refCoordBoundary * this->getBoundaryNodeRefCoord(1);   
  }

};




//Rectangle (0,0) to (lx,ly)
// DOFs in x: Nx, y:Ny
// meshsize: (hx,hy) = (lx/(Nx-1), ly/(Ny-1))
template< typename DataTypeContainer, typename QuocElementType = QuocElement2D<DataTypeContainer>, typename QuocBoundaryElementType = QuocBoundaryElement2D<DataTypeContainer> >
class QuocMesh2D {
public:
  typedef QuocElementType ElementType;
  typedef QuocBoundaryElementType                   BoundaryElementType;
  typedef typename DataTypeContainer::RealType      RealType;
  typedef typename DataTypeContainer::PointType     PointType;
  typedef typename DataTypeContainer::Indices2DType Indices2DType;
  typedef Indices2DType                             IndicesOfElementType;
  typedef typename DataTypeContainer::VectorType    VectorType;
  typedef typename DataTypeContainer::MaskType      MaskType;
  typedef typename DataTypeContainer::IntVecType    IntVecType;
  
  static const int _VTKCELLTYPE = 8; //VTK_PIXEL
  
protected:
  const int _Nx, _Ny;
  IntVecType _NumDofVec;
  const RealType _lx, _ly;
  PointType _LengthVec;
  const RealType _hx, _hy;
  std::vector< PointType > _vertexIterator;
  std::vector< ElementType > _elementIterator;
  std::vector< BoundaryElementType > _boundaryElementIterator;
public :
  MaskType _boundaryTop, _boundaryBottom, _boundaryLeft, _boundaryRight;
  MaskType _boundaryPeriodic;
  std::vector<int> _periodicIdentificationIndices;
  
public:
  //! Create empty QuocMesh
  QuocMesh2D ( ) : _vertexIterator(), _elementIterator(), _boundaryElementIterator()  { }
  
  QuocMesh2D ( const int Nx, const int Ny, const RealType lx, const RealType ly ) : 
    _Nx (Nx), _Ny(Ny), 
    _lx( lx ), _ly( ly ),
    _hx( _lx / static_cast<RealType>(_Nx - 1 ) ), _hy( _ly / static_cast<RealType>(_Ny - 1 ) )
  {
        _NumDofVec[0] = _Nx; _NumDofVec[1] = _Ny;
        _LengthVec[0] = _lx; _LengthVec[1] = _ly;
      
        const int numGlobalDofs = _Nx * _Ny;
        _boundaryTop.resize( numGlobalDofs, false ); _boundaryBottom.resize( numGlobalDofs, false );
        _boundaryLeft.resize( numGlobalDofs, false ); _boundaryRight.resize( numGlobalDofs, false );
        _boundaryPeriodic.resize( numGlobalDofs, false );
        
        for( int yIdx=0; yIdx<_Ny; ++yIdx )
            for( int xIdx=0; xIdx<_Nx; ++xIdx ){
              PointType vertex ( _lx * xIdx / static_cast<RealType>( _Nx - 1 ),  _ly * yIdx / static_cast<RealType>( _Ny - 1 ) );
              int nodeIdx = this->pushBackVertex( vertex );
              if( xIdx == 0 ) _boundaryLeft[nodeIdx] = true;
              if( xIdx == _Nx - 1 ) _boundaryRight[nodeIdx] = true;
              if( yIdx == 0 ) _boundaryBottom[nodeIdx] = true;
              if( yIdx == _Ny - 1 ) _boundaryTop[nodeIdx] = true;
              if( ( xIdx == _Nx - 1 ) || ( yIdx == _Ny - 1 ) ) _boundaryPeriodic[nodeIdx] = true;
          }
          
      for( int nodeIdx = 0; nodeIdx < numGlobalDofs; ++nodeIdx ){
          if( (!_boundaryRight[nodeIdx]) && (! _boundaryTop[nodeIdx]) ){
              //! \note ordering of elements has to coincide with ordering of base functions
              const int numNewElement = this->pushBackElement( Indices2DType(nodeIdx, nodeIdx + 1, nodeIdx + _Nx, nodeIdx + _Nx + 1 ) );  
              if( this->getVertex(nodeIdx)[0] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::LEFT, _hy );
              if( this->getVertex(nodeIdx + 1)[0] == _lx ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::RIGHT, _hy );
              if( this->getVertex(nodeIdx)[1] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::BOTTOM, _hx );
              if( this->getVertex(nodeIdx + _Nx)[1] == _ly ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::TOP, _hx );
          }
      }
      
      //Periodic indices to identify facing presenting node
      _periodicIdentificationIndices.resize( numGlobalDofs );
      for( int yIdx=0; yIdx<_Ny; ++yIdx )
          for( int xIdx=0; xIdx<_Nx; ++xIdx ){
                int periodicXIdx = xIdx, periodicYIdx = yIdx;
                if( xIdx == _Nx - 1 ) periodicXIdx = 0;
                if( yIdx == _Ny - 1 ) periodicYIdx = 0;
                _periodicIdentificationIndices[this->getGlobalNodeIndex(xIdx,yIdx)] = this->getGlobalNodeIndex(periodicXIdx,periodicYIdx);
          }
      
  }
  
  QuocMesh2D ( const IntVecType &NumDofVec, const PointType &LenghtVec ) : QuocMesh2D ( NumDofVec[0], NumDofVec[1], LenghtVec[0], LenghtVec[1] ) {} 
    
  virtual ~QuocMesh2D ( ) {} 
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumElements ( ) const { return ( static_cast<int> ( _elementIterator.size() ) );}
  int getNumBoundaryElements ( ) const { return ( static_cast<int> ( _boundaryElementIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new element and return global index
  int pushBackElement ( const Indices2DType nodeIdx ) {
    int globalIdx = getNumElements();
    _elementIterator.push_back ( QuocElementType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }
  
  int pushBackBoundaryElement ( const int elementIndex, QuocBoundaryType boundaryType, const RealType h ) {
    int globalBoundaryIdx = getNumBoundaryElements();
    _boundaryElementIterator.push_back ( BoundaryElementType( globalBoundaryIdx, this->getElement(elementIndex), boundaryType, h ) );
    return globalBoundaryIdx;
  }

  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  int getGlobalNodeIndex ( const int xIdx, const int yIdx ) const {  return ( xIdx + yIdx *_Nx );}
  int getGlobalElementIndex ( const int xFac, const int yFac ) const { return ( xFac + yFac * (_Nx - 1 ) );} 
  
  const QuocElementType& getElement ( const int num ) const {return _elementIterator[num];}
  QuocElementType & getElement ( const int num ) {return _elementIterator[num];}
  void setElement ( const int num, const QuocElementType Arg ) { _elementIterator[num] = Arg;}

  const BoundaryElementType& getBoundaryElement ( const int num ) const {return _boundaryElementIterator[num];}
  BoundaryElementType & getBoundaryElement ( const int num ) {return _boundaryElementIterator[num];}
  
  
  int getElementNodeIdx  ( const int num, const int localNode ) const { return _elementIterator[num].getGlobalNodeIdx(localNode);}
  void setElementNodeIdx ( const int num, const int localNode, const int newIdx ) {_elementIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  
  const RealType getMeshSize( const int direction ) const{ 
    if( direction == 0 ) return _hx;
    if( direction == 1 ) return _hy;
    if( direction == 2 ) return 0.;
    throw std::invalid_argument( aol::strprintf ( "Wrong direction in getMeshSize in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  const RealType getInterfaceWith( ) const{ return std::sqrt(_hx * _hy ); }
  const RealType getWidth( const int direction ) const {
   switch( direction ){
        case 0 :    return _lx;    break;
        case 1 : return _ly; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  const int getNumDofs( const int direction ) const {
   switch( direction ){
        case 0 : return _Nx; break;
        case 1 : return _Ny; break;
        case 2 : return 1; break; //for compatibility with vtk
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel in QuocMesh2D get num dofs. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  void getNumDofVec( IntVecType &vec ) const { vec = _NumDofVec;}
  void getLenghtVec( PointType &vec ) const { vec = _LengthVec;}
  
};


 














template<typename DataTypeContainer>
class QuocElement3D {

  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::DomVecType      DomVecType;
  typedef typename DataTypeContainer::PointType       PointType;
  typedef typename DataTypeContainer::Indices3DType   Indices3DType;
  typedef typename DataTypeContainer::MaskType        MaskType;
  typedef std::vector<PointType>                      VertexIterator;

  public :   
    // global indices of element an nodes
    int _globIdx;
    Indices3DType _globNodeIdx;
    PointType _nodes[8];
    
public:
  
  QuocElement3D () : _globIdx(-1) {}
  
  QuocElement3D( const int globalIdx, const Indices3DType globalNodeIndex, const VertexIterator &nodes  ) : 
     _globIdx(globalIdx),
     _globNodeIdx ( globalNodeIndex )
    {
        for ( int i = 0; i < 8; ++i ) _nodes[i] = nodes[ _globNodeIdx[i] ];
    }
  
  ~QuocElement3D(){}

  // get and set functions
  int getGlobalElementIdx(  ) const { return _globIdx;}
  const Indices3DType & getGlobalNodeIdx( ) const{ return _globNodeIdx;}
  int getGlobalNodeIdx(int localIndex) const{ return _globNodeIdx[localIndex];}
  void setGlobalNodeIdx(int localIndex, int globalIndex) { _globNodeIdx[localIndex] = globalIndex;}
  
  const PointType& getNode ( int i ) const { return _nodes[i];}
  PointType& getNode ( int i ) { return _nodes[i];}
  void setNode ( int i, const PointType& node ) {_nodes[i] = node;}
};


template<typename DataTypeContainer>
class QuocBoundaryElement3D {

  typedef typename DataTypeContainer::RealType        RealType;
  typedef typename DataTypeContainer::DomVecType      DomVecType;
  typedef typename DataTypeContainer::DomVecTypeBoundary DomVecTypeBoundary;
  typedef typename DataTypeContainer::PointType     PointType;
  typedef typename DataTypeContainer::Indices3DType   Indices3DType;
  typedef std::vector<PointType>                    VertexIterator;
  typedef QuocElement3D<DataTypeContainer>            ElementType;

  protected :   
    int _globBoundaryIdx;
    const ElementType _element;
    const QuocBoundaryType _boundaryType;
    DomVecType _normal;
    PointType _boundaryNodes[4];
    PointType _boundaryNodesRefCoord[4];
    const RealType _h;
    int _nodeIndicesOf3DElement[4];
    
public:
  
  QuocBoundaryElement3D () : _globBoundaryIdx(-1) {}
  
  QuocBoundaryElement3D( const int globalBoundaryIdx, const ElementType & element, const QuocBoundaryType &boundaryType, const RealType h  ) : 
     _globBoundaryIdx(globalBoundaryIdx), _element (element), _boundaryType ( boundaryType ), _h(h) {
         
         _normal.setZero();
         switch( boundaryType ){
             case LEFT:{
               _nodeIndicesOf3DElement[0] = 0; _nodeIndicesOf3DElement[1] = 2; _nodeIndicesOf3DElement[2] = 4; _nodeIndicesOf3DElement[3] = 6;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 0.0; 
               _boundaryNodesRefCoord[1] << 0.0, 1.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 0.0, 1.0; 
               _boundaryNodesRefCoord[3] << 0.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(0);
               _boundaryNodes[1] = element.getNode(2);
               _boundaryNodes[2] = element.getNode(4);
               _boundaryNodes[3] = element.getNode(6);
               _normal[0] = -1.;
             }break;
                 
             case RIGHT:{
               _nodeIndicesOf3DElement[0] = 1; _nodeIndicesOf3DElement[1] = 3; _nodeIndicesOf3DElement[2] = 5; _nodeIndicesOf3DElement[3] = 7;
               _boundaryNodesRefCoord[0] << 1.0, 0.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 1.0, 0.0;
               _boundaryNodesRefCoord[2] << 1.0, 0.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(1);
               _boundaryNodes[1] = element.getNode(3);
               _boundaryNodes[2] = element.getNode(5);
               _boundaryNodes[3] = element.getNode(7);
               _normal[0] = 1.;
               
             }break;
                 
             case BOTTOM:{
               _nodeIndicesOf3DElement[0] = 0; _nodeIndicesOf3DElement[1] = 1; _nodeIndicesOf3DElement[2] = 2; _nodeIndicesOf3DElement[3] = 3;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 0.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 1.0, 0.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 0.0;
               _boundaryNodes[0] = element.getNode(0);
               _boundaryNodes[1] = element.getNode(1);
               _boundaryNodes[2] = element.getNode(2);
               _boundaryNodes[3] = element.getNode(3);
               _normal[2] = -1.;
             }break;
                 
             case TOP:{
               _nodeIndicesOf3DElement[0] = 4; _nodeIndicesOf3DElement[1] = 5; _nodeIndicesOf3DElement[2] = 6; _nodeIndicesOf3DElement[3] = 7;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 1.0;
               _boundaryNodesRefCoord[1] << 1.0, 0.0, 1.0;
               _boundaryNodesRefCoord[2] << 0.0, 1.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(4);
               _boundaryNodes[1] = element.getNode(5);
               _boundaryNodes[2] = element.getNode(6);
               _boundaryNodes[3] = element.getNode(7);
               _normal[2] = 1.;
             }break;
             
             case FRONT:{
               _nodeIndicesOf3DElement[0] = 0; _nodeIndicesOf3DElement[1] = 1; _nodeIndicesOf3DElement[2] = 4; _nodeIndicesOf3DElement[3] = 5;
               _boundaryNodesRefCoord[0] << 0.0, 0.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 0.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 0.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 0.0, 1.0;
               _boundaryNodes[0] = element.getNode(0);
               _boundaryNodes[1] = element.getNode(1);
               _boundaryNodes[2] = element.getNode(4);
               _boundaryNodes[3] = element.getNode(5);
               _normal[1] = -1.;
             }break;
                 
             case BACK:{
               _nodeIndicesOf3DElement[0] = 2; _nodeIndicesOf3DElement[1] = 3; _nodeIndicesOf3DElement[2] = 6; _nodeIndicesOf3DElement[3] = 7;
               _boundaryNodesRefCoord[0] << 0.0, 1.0, 0.0;
               _boundaryNodesRefCoord[1] << 1.0, 1.0, 0.0;
               _boundaryNodesRefCoord[2] << 0.0, 1.0, 1.0; 
               _boundaryNodesRefCoord[3] << 1.0, 1.0, 1.0;
               _boundaryNodes[0] = element.getNode(2);
               _boundaryNodes[1] = element.getNode(3);
               _boundaryNodes[2] = element.getNode(6);
               _boundaryNodes[3] = element.getNode(7);
               _normal[1] = 1.;
             }break;
                 
             default:
                 break;
         }
    }
  
  ~QuocBoundaryElement3D(){}
  
  const ElementType & getElement() const {return _element;}
  int getNodeIndexOfElement(const int index ) const{ return _nodeIndicesOf3DElement[index];}
  const DomVecType& getNormal() const {return _normal;}
  const PointType& getBoundaryNode ( int i ) const { return _boundaryNodes[i];}
  PointType& getBoundaryNode ( int i ) { return _boundaryNodes[i];}
  const PointType& getBoundaryNodeRefCoord ( int i ) const { return _boundaryNodesRefCoord[i];}
  PointType& getBoundaryNodeRefCoord ( int i ) { return _boundaryNodesRefCoord[i];}
  const RealType getH() const {return _h;}
  const QuocBoundaryType getBoundaryType() const {return _boundaryType;}
  void getRefCoord( const DomVecTypeBoundary &refCoordBoundary, DomVecType & refCoord ) const {
    refCoord = (1. - refCoordBoundary[1]) * ( ( 1. - refCoordBoundary[0] ) * this->getBoundaryNodeRefCoord(0) + refCoordBoundary[0] * this->getBoundaryNodeRefCoord(1) )
                   + refCoordBoundary[1]  * ( ( 1. - refCoordBoundary[0] ) * this->getBoundaryNodeRefCoord(2) + refCoordBoundary[0] * this->getBoundaryNodeRefCoord(3) ); 
  }

};








//Rectangle (0,0,0) to (lx,ly,ly)
// DOFs in x: Nx, y:Ny, z:Nz
// meshsize: (hx,hy,hz) = (lx/(Nx-1), ly/(Ny-1), lz/(Nz-1))
template< typename DataTypeContainer, typename QuocElementType = QuocElement3D<DataTypeContainer>, typename QuocBoundaryElementType = QuocBoundaryElement3D<DataTypeContainer> >
class QuocMesh3D {
public:
  typedef QuocElementType ElementType;
  typedef QuocBoundaryElementType BoundaryElementType;
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef typename DataTypeContainer::Indices3DType Indices3DType;
  typedef Indices3DType IndicesOfElementType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::MaskType MaskType;
  typedef typename DataTypeContainer::IntVecType IntVecType;
  
  static const int _VTKCELLTYPE = 11; //VTK_VOXEL
  
protected:
  const int _Nx, _Ny, _Nz;
  IntVecType _NumDofVec;
  const RealType _lx, _ly, _lz;
  PointType _LengthVec;
  const RealType _hx, _hy, _hz;
  std::vector< PointType > _vertexIterator;
  std::vector< ElementType > _elementIterator;
  std::vector< BoundaryElementType > _boundaryElementIterator;
public:
  MaskType _boundaryFront, _boundaryBack, _boundaryTop, _boundaryBottom, _boundaryLeft, _boundaryRight;
  MaskType _boundaryPeriodic;
  std::vector<int> _periodicIdentificationIndices;
  
public:
  //! Create empty QuocMesh
  QuocMesh3D ( ) : _vertexIterator(), _elementIterator()  { }
  
  QuocMesh3D ( const int Nx, const int Ny, const int Nz, const RealType lx, const RealType ly, const RealType lz ) : 
    _Nx (Nx), _Ny(Ny), _Nz(Nz),
    _lx( lx ), _ly( ly ), _lz( lz ),
    _hx( _lx / static_cast<RealType>(_Nx - 1 ) ), _hy( _ly / static_cast<RealType>(_Ny - 1 ) ), _hz( _lz / static_cast<RealType>(_Nz - 1 ) )
    {
        
        _NumDofVec[0] = _Nx; _NumDofVec[1] = _Ny; _NumDofVec[2] = _Nz;
        _LengthVec[0] = _lx; _LengthVec[1] = _ly; _LengthVec[2] = _lz;
        
        const int numGlobalDofs = _Nx * _Ny * _Nz;
        _boundaryFront.resize( numGlobalDofs, false ); _boundaryBack.resize( numGlobalDofs, false );
        _boundaryTop.resize( numGlobalDofs, false ); _boundaryBottom.resize( numGlobalDofs, false );
        _boundaryLeft.resize( numGlobalDofs, false ); _boundaryRight.resize( numGlobalDofs, false );
        _boundaryPeriodic.resize( numGlobalDofs, false );
        
        for( int zIdx=0; zIdx<_Nz; ++zIdx )
            for( int yIdx=0; yIdx<_Ny; ++yIdx )
                for( int xIdx=0; xIdx<_Nx; ++xIdx ){
                    PointType vertex ( _lx * xIdx / static_cast<RealType>( _Nx - 1 ), 
                                         _ly * yIdx / static_cast<RealType>( _Ny - 1 ),
                                         _lz * zIdx / static_cast<RealType>( _Nz - 1 ) 
                                       );
                    int nodeIdx = this->pushBackVertex( vertex );
                    if( xIdx == 0 ) _boundaryLeft[nodeIdx] = true;
                    if( xIdx == _Nx - 1 ) _boundaryRight[nodeIdx] = true;
                    if( yIdx == 0 ) _boundaryFront[nodeIdx] = true;
                    if( yIdx == _Ny - 1 ) _boundaryBack[nodeIdx] = true;
                    if( zIdx == 0 ) _boundaryBottom[nodeIdx] = true;
                    if( zIdx == _Nz - 1 ) _boundaryTop[nodeIdx] = true;
                    
                    if( ( xIdx == _Nx - 1 ) || ( yIdx == _Ny - 1 ) || ( zIdx == _Nz - 1 ) ) _boundaryPeriodic[nodeIdx] = true;
          }
          
      for( int nodeIdx = 0; nodeIdx < numGlobalDofs; ++nodeIdx ){
          if( (!_boundaryRight[nodeIdx]) && (!_boundaryBack[nodeIdx]) && (!_boundaryTop[nodeIdx]) ){
              Indices3DType indices;
              indices << nodeIdx, nodeIdx + 1, nodeIdx + _Nx, nodeIdx + _Nx + 1,
                         nodeIdx + _Nx * _Ny, nodeIdx + _Nx * _Ny + 1, nodeIdx + _Nx * _Ny + _Nx, nodeIdx + _Nx * _Ny + _Nx + 1;
              const int numNewElement = this->pushBackElement( indices );   
              
              if( this->getVertex(nodeIdx)[0] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::LEFT, _hy * _hz );
              if( this->getVertex(nodeIdx + 1)[0] == _lx ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::RIGHT, _hy * _hz );
              if( this->getVertex(nodeIdx)[1] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::FRONT, _hx * _hz );
              if( this->getVertex(nodeIdx + _Nx)[1] == _ly ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::BACK, _hx * _hz );
              if( this->getVertex(nodeIdx)[2] == 0 ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::BOTTOM, _hx * _hy );
              if( this->getVertex(nodeIdx + _Nx * _Ny)[2] == _lz ) this->pushBackBoundaryElement( numNewElement, QuocBoundaryType::TOP, _hx * _hy );
          }
      }

      //Periodic indices
      _periodicIdentificationIndices.resize( numGlobalDofs, -1 );
      for( int zIdx=0; zIdx<_Nz; ++zIdx )
          for( int yIdx=0; yIdx<_Ny; ++yIdx )
             for( int xIdx=0; xIdx<_Nx; ++xIdx ){
                int periodicXIdx = xIdx, periodicYIdx = yIdx, periodicZIdx = zIdx;
                if( xIdx == _Nx - 1 ) periodicXIdx = 0;
                if( yIdx == _Ny - 1 ) periodicYIdx = 0;
                if( zIdx == _Nz - 1) periodicZIdx = 0;
                _periodicIdentificationIndices[this->getGlobalNodeIndex(xIdx,yIdx,zIdx)] = this->getGlobalNodeIndex(periodicXIdx,periodicYIdx,periodicZIdx);
             }


  }
  
  
  QuocMesh3D ( const IntVecType &NumDofVec, const PointType &LenghtVec ) : QuocMesh3D ( NumDofVec[0], NumDofVec[1], NumDofVec[2], LenghtVec[0], LenghtVec[1], LenghtVec[2] ) {} 
    
  virtual ~QuocMesh3D (  ) {} 
    
public:
  int getNumVertices ( ) const { return ( _vertexIterator.size() );}
  int getNumElements ( ) const { return ( static_cast<int> ( _elementIterator.size() ) );}
  int getNumBoundaryElements ( ) const { return ( static_cast<int> ( _boundaryElementIterator.size() ) );}
  
  //! insert new vertex and return global index
  int pushBackVertex ( const PointType newVertex ) {
    _vertexIterator.push_back ( newVertex );
    return getNumVertices() - 1;
  }
  //! insert new element and return global index
  int pushBackElement ( const Indices3DType nodeIdx ) {
    int globalIdx = getNumElements();
    _elementIterator.push_back ( ElementType( globalIdx, nodeIdx, _vertexIterator ) );
    return globalIdx;
  }
  
  int pushBackBoundaryElement ( const int elementIndex, QuocBoundaryType boundaryType, const RealType h ) {
    int globalBoundaryIdx = getNumBoundaryElements();
    _boundaryElementIterator.push_back ( BoundaryElementType( globalBoundaryIdx, this->getElement(elementIndex), boundaryType, h) );
    return globalBoundaryIdx;
  }

  const PointType& getVertex ( const int num ) const { return _vertexIterator[num];}
  void setVertex ( const int num, const PointType Arg ) { _vertexIterator[num] = Arg;}

  int getGlobalNodeIndex ( const int xIdx, const int yIdx, const int zIdx ) const { return ( xIdx + yIdx * _Nx + zIdx * _Nx * _Ny );}
  int getGlobalElementIndex ( const int xFac, const int yFac, const int zFac ) const { return ( xFac + yFac * (_Nx - 1 ) + zFac * (_Nx - 1 ) * ( _Ny - 1 ) );}  

  const ElementType& getElement ( const int num ) const {return _elementIterator[num];}
  ElementType & getElement ( const int num ) {return _elementIterator[num];}
  void setElement ( const int num, const ElementType Arg ) { _elementIterator[num] = Arg;}

  const QuocBoundaryElementType& getBoundaryElement ( const int num ) const {return _boundaryElementIterator[num];}
  QuocBoundaryElementType & getBoundaryElement ( const int num ) {return _boundaryElementIterator[num];}
  
  int getElementNodeIdx  ( const int num, const int localNode ) const { return _elementIterator[num].getGlobalNodeIdx(localNode);}
  void setElementNodeIdx ( const int num, const int localNode, const int newIdx ) {_elementIterator[num].setGlobalNodeIdx(localNode, newIdx );}
  
  const RealType getMeshSize( const int direction ) const{ 
    if( direction == 0 ) return _hx;
    if( direction == 1 ) return _hy;
    if( direction == 2 ) return _hz;
    throw std::invalid_argument( aol::strprintf ( "Wrong direction in getMeshSize in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
  }
  
  const RealType getInterfaceWith( ) const{ return std::cbrt(_hx * _hy * _hz);}
  const RealType getWidth( const int direction ) const {
     switch( direction ){
        case 0 : return _lx;    break;
        case 1 : return _ly; break;
        case 2 : return _lz; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  const int getNumDofs( const int direction ) const {
   switch( direction ){
        case 0 : return _Nx; break;
        case 1 : return _Ny; break;
        case 2 : return _Nz; break;
        default: throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); break;
    }
  }
  
  void getNumDofVec( IntVecType &vec ) const { vec = _NumDofVec;}
  void getLenghtVec( PointType &vec ) const { vec = _LengthVec;}
};

}//end namespace

#endif
