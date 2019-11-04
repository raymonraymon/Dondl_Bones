#ifndef __QUOCCONFIGURATORS_H
#define __QUOCCONFIGURATORS_H

#include <general.h>
#include <quocBasefunctionSets.h>
#include <quocMesh.h>
#include <quocQuadrature.h>

namespace quocFE {
    
    
template <typename DataTypeContainer,
          typename MeshType = QuocMesh1D<DataTypeContainer>,
          typename QuadType = SimpsonQuadrature1DVec<typename DataTypeContainer::RealType, typename DataTypeContainer::DomVecType> >
class QuocConfigurator1D {
public:
    static const int dimDomain = 1;
    static const int maxNumLocalDofs = 2;

  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
  typedef MeshType                                                      InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  
  typedef QuocCachedBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  //! globally affine functions with symmetric gradient
  static const int numAffineSymGradDofs = 1;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >         LocalMatrixTypeMixed;
  typedef QuocGlobalAffineSymGradBaseFunctionSet2D<DataTypeContainer> GlobalAffineSymGradBaseFuncSet;

protected:  
  const MeshType &_mesh;
  const QuadType _quad;
  const BaseFuncSetType _baseFuncSet;
  const RealType _volElement;
  
public:
  
  QuocConfigurator1D ( const InitType &Mesh ) : _mesh ( Mesh ), _baseFuncSet( _mesh.getMeshSize(0) ), _volElement( _mesh.getMeshSize(0) ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &El, int localIndex ) const {return El.getGlobalNodeIdx( localIndex );}

  const RealType getVolOfElement() const {return _volElement;}
  
  void getGlobalCoords ( const ElementType &El, const DomVecType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord[0]; 
      Coord *= _mesh.getMeshSize(0);
      Coord += El.getNode(0); //Coord Left Down
  }
  
  void getGlobalCoords ( const ElementType &El, const int QuadPoint, PointType &Coord ) const {
      const DomVecType & LocalCoord = _quad.getRefCoord ( QuadPoint ); 
      this->getGlobalCoords( El, LocalCoord, Coord );
  }
  
  void getLocalCoords ( const PointType &GlobalCoord, int &elementNumber, PointType &LocalCoord ) const {
      PointType BasePoint;
      std::vector<int> BasePointFac ( 1 );
      for ( int comp = 0; comp < 1; ++comp ) {
          const RealType sc = GlobalCoord[comp] / _mesh.getMeshSize(comp);
          BasePointFac[comp] = static_cast<int> ( sc );
          BasePoint[comp] = BasePointFac[comp] * _mesh.getMeshSize(comp);
          LocalCoord[comp] = sc - BasePointFac[comp];
          //special case: "right" boundary
          if( std::abs( _mesh.getWidth( comp ) - GlobalCoord[comp] ) < 2.e-16 ){
              BasePoint[comp] = _mesh.getWidth(comp) - _mesh.getMeshSize(comp);
              BasePointFac[comp]  = _mesh.getNumDofs( comp ) - 2.;
              LocalCoord[comp] = 1.;
          }
      }
      elementNumber = _mesh.getGlobalElementIndex( BasePointFac[0] );
  }

};
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
template <typename DataTypeContainer,
          typename MeshType = QuocMesh2D<DataTypeContainer>,
          typename QuadType = SimpsonQuadrature2D<typename DataTypeContainer::RealType, typename DataTypeContainer::DomVecType>,
          typename BdrQuadType = SimpsonQuadrature1D<typename DataTypeContainer::RealType> >
class QuocConfigurator2D {
public:
    static const int dimDomain = 2;
    static const int maxNumLocalDofs = 4;
    static const int maxNumLocalBoundaryDofs = 2;

  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
  typedef BdrQuadType                                                   BoundaryQuadType;
  typedef MeshType                                                      InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename MeshType::BoundaryElementType                        BoundaryElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::DomVecTypeBoundary                DomVecTypeBoundary;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  
  typedef QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  //! globally affine functions with symmetric gradient
  static const int numAffineSymGradDofs = 3;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >         LocalMatrixTypeMixed;
  typedef QuocGlobalAffineSymGradBaseFunctionSet2D<DataTypeContainer> GlobalAffineSymGradBaseFuncSet;

protected:  
  const MeshType &_mesh;
  const QuadType _quad;
  const BaseFuncSetType _baseFuncSet;
  const RealType _volElement;
  
public:
  
  QuocConfigurator2D ( const InitType &Mesh ) : _mesh ( Mesh ), _baseFuncSet( _mesh.getMeshSize(0), _mesh.getMeshSize(1) ), _volElement( _mesh.getMeshSize(0) * _mesh.getMeshSize(1) ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs (  ) const { return maxNumLocalDofs;}
  inline int getNumLocalBoundaryDofs (  ) const { return maxNumLocalBoundaryDofs;}

  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  int maxNumBoundaryQuadPoints() const {return BoundaryQuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &El, int localIndex ) const {return El.getGlobalNodeIdx( localIndex );}

  const RealType getVolOfElement() const {return _volElement;}
  
  void getGlobalCoords ( const ElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1);
      Coord += El.getNode(0); //Coord Left Down
  }
  
  void getGlobalCoords ( const ElementType &El, const int QuadPoint, PointType &Coord ) const {
      const PointType & LocalCoord = _quad.getRefCoord ( QuadPoint ); 
      this->getGlobalCoords( El, LocalCoord, Coord );
  }
  
  void getLocalCoords ( const PointType &GlobalCoord, int &elementNumber, PointType &LocalCoord ) const {
      PointType BasePoint;
      std::vector<int> BasePointFac ( 2 );
      for ( int comp = 0; comp < 2; ++comp ) {
          const RealType sc = GlobalCoord[comp] / _mesh.getMeshSize(comp);
          BasePointFac[comp] = static_cast<int> ( sc );
          BasePoint[comp] = BasePointFac[comp] * _mesh.getMeshSize(comp);
          LocalCoord[comp] = sc - BasePointFac[comp];
          //special case: "right" boundary
          if( std::abs( _mesh.getWidth( comp ) - GlobalCoord[comp] ) < 2.e-16 ){
              BasePoint[comp] = _mesh.getWidth(comp) - _mesh.getMeshSize(comp);
              BasePointFac[comp]  = _mesh.getNumDofs( comp ) - 2.;
              LocalCoord[comp] = 1.;
          }
      }
      
      elementNumber = _mesh.getGlobalElementIndex( BasePointFac[0], BasePointFac[1] );
      
  }
  
  void getGlobalCoordsForBoundaryElement ( const BoundaryElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1);
      Coord += El.getElement().getNode(0); //Coord Left Down
  }

};






template <typename DataTypeContainer, 
          typename MeshType = QuocMesh3D<DataTypeContainer>,
          typename QuadType = SimpsonQuadrature3D<typename DataTypeContainer::RealType, typename DataTypeContainer::DomVecType>,
          typename BdrQuadType = SimpsonQuadrature2D<typename DataTypeContainer::RealType, typename DataTypeContainer::DomVecTypeBoundary> >
class QuocConfigurator3D {
public:
  static const int dimDomain = 3;
  static const int maxNumLocalDofs = 8;
  static const int maxNumLocalBoundaryDofs = 4;

  typedef DataTypeContainer                                             DTContainer;
  typedef QuadType                                                      QuadRuleType;
  typedef BdrQuadType                                                   BoundaryQuadType;
  typedef MeshType                                                      InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef typename MeshType::ElementType                                ElementType;
  typedef typename MeshType::BoundaryElementType                        BoundaryElementType;
  typedef typename DataTypeContainer::RealType                          RealType;
  typedef typename DataTypeContainer::DomVecType                        DomVecType;
  typedef typename DataTypeContainer::DomVecTypeBoundary                DomVecTypeBoundary;
  typedef typename DataTypeContainer::PointType                         PointType;
  typedef typename DataTypeContainer::Matrix22                          Matrix22;
  typedef typename DataTypeContainer::Matrix32                          Matrix32;
  typedef typename DataTypeContainer::Matrix33                          Matrix33;
  typedef typename DataTypeContainer::DerivativeVectorValuedType        DerivativeVectorValuedType;
  typedef typename DataTypeContainer::VectorType                        VectorType;
  typedef typename DataTypeContainer::FullMatrixType                    FullMatrixType;
  typedef typename DataTypeContainer::TripletType                       TripletType;
  typedef typename DataTypeContainer::SparseMatrixType                  SparseMatrixType;
  typedef typename DataTypeContainer::MaskType                          MaskType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, 1 >                  LocalVectorType;
  typedef Eigen::Matrix<RealType, maxNumLocalDofs, maxNumLocalDofs >    LocalMatrixType;
  
  typedef QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> BaseFuncSetType;
  
  //! globally affine functions with symmetric gradient
  static const int numAffineSymGradDofs = 6;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, numAffineSymGradDofs >    LocalMatrixTypeAffineSymGrad;
  typedef Eigen::Matrix<RealType, numAffineSymGradDofs, maxNumLocalDofs >    LocalMatrixTypeMixed;
  typedef QuocGlobalAffineSymGradBaseFunctionSet3D<DataTypeContainer> GlobalAffineSymGradBaseFuncSet;

protected:  
  const MeshType &_mesh;
  const QuadType _quad;
  mutable BaseFuncSetType _baseFuncSet;
  const RealType _volElement;
  
public:
  
  QuocConfigurator3D ( const InitType &Mesh ) : _mesh ( Mesh ), _baseFuncSet( _mesh.getMeshSize(0), _mesh.getMeshSize(1), _mesh.getMeshSize(2) ),
  _volElement( _mesh.getMeshSize(0) * _mesh.getMeshSize(1) * _mesh.getMeshSize(2) ) {}

  const InitType& getInitializer( ) const { return this->_mesh; }

  inline int getNumLocalDofs ( const ElementType & ) const { return maxNumLocalDofs;}
  inline int getNumLocalDofs ( ) const { return maxNumLocalDofs;}
  inline int getNumLocalBoundaryDofs (  ) const { return maxNumLocalBoundaryDofs;}
  
  int getNumGlobalDofs( ) const {return this->_mesh.getNumVertices();}
  int maxNumQuadPoints( ) const { return QuadType::numQuadPoints;}
  int maxNumBoundaryQuadPoints() const {return BoundaryQuadType::numQuadPoints;}
  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &/*T*/ ) const {return _baseFuncSet;}
  const BaseFuncSetType& getBaseFunctionSet ( ) const {return _baseFuncSet;}
  
  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {return T.getGlobalNodeIdx( localIndex );}

  const RealType getVolOfElement() const {return _volElement;}
  
  void getGlobalCoords ( const ElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1); Coord[2] *= _mesh.getMeshSize(2);
      Coord += El.getNode(0); //Coord Left Down
  }
  
  void getGlobalCoords ( const ElementType &El, const int QuadPoint, PointType &Coord ) const {
      const PointType & LocalCoord = _quad.getRefCoord ( QuadPoint ); 
      this->getGlobalCoords( El, LocalCoord, Coord );
  }
  
  void getLocalCoords ( const PointType &GlobalCoord, int &elementNumber, PointType &LocalCoord ) const {
      PointType BasePoint;
      std::vector<int> BasePointFac ( 3 );
      
      for ( int comp = 0; comp < 3; ++comp ) {
          
         bool leftBoundary = false, rightBoundary = false, outside = false;
         if( GlobalCoord[comp] < 2.e-16 ) { leftBoundary = true; outside = true; }
         if( GlobalCoord[comp] > (_mesh.getWidth( comp ) - 2.e-16) ) { rightBoundary = true; outside = true;}
         
         if( outside ){
            if( leftBoundary ){
                BasePoint[comp] = 0.;
                BasePointFac[comp]  = 0.;
                LocalCoord[comp] = 0.;   
            }
            if( rightBoundary ){
                BasePoint[comp] = _mesh.getWidth(comp) - _mesh.getMeshSize(comp);
                BasePointFac[comp]  = _mesh.getNumDofs( comp ) - 2.;
                LocalCoord[comp] = 1.;   
            }
          }else{
              const RealType relativeCoord = GlobalCoord[comp] / _mesh.getMeshSize(comp);
              BasePointFac[comp] = static_cast<int> ( relativeCoord );
              BasePoint[comp] = BasePointFac[comp] * _mesh.getMeshSize(comp);
              LocalCoord[comp] = relativeCoord - BasePointFac[comp];
          }
      }
      
      elementNumber = _mesh.getGlobalElementIndex( BasePointFac[0], BasePointFac[1], BasePointFac[2] );
      
  }
  
  void getGlobalCoordsForBoundaryElement ( const BoundaryElementType &El, const PointType &LocalCoord, PointType &Coord ) const {
      Coord = LocalCoord; 
      Coord[0] *= _mesh.getMeshSize(0); Coord[1] *= _mesh.getMeshSize(1); Coord[2] *= _mesh.getMeshSize(2);
      Coord += El.getElement().getNode(0); //Coord Left Down
  }
};

} //end namespace

#endif //__QUOCMESHCONFIGURATORS_H
