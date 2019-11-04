#ifndef __QUOCINTEGRATORS_H
#define __QUOCINTEGRATORS_H

// #include <tensor.h>
#include <quocDiscreteFunction.h>
  
namespace quocFE {
  
//!===========================================================================================================================
//! FE OPERATOR INTERFACES (for any quadrature type specified in the template argument)
//!===========================================================================================================================


//!===========================================================================================================================
//! Scalar-Valued Intefaces
//!===========================================================================================================================
  
//! Integrator to compute \f$\int_\Omega f(...) dx\f$, where \f$\f$ is the argument of the operator.
template <typename ConfiguratorType, typename Imp>
class QuocIntegrator{
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
protected:
  const ConfiguratorType &_config;
public:

  QuocIntegrator ( const ConfiguratorType & Config ) : _config( Config ) {}
  
  virtual ~QuocIntegrator( ) {}

  void assembleAdd ( RealType &Dest ) const {

    RealType res = 0.;

    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      RealType a = 0.;
      for ( int q = 0; q < numQuadPoints; ++q ) {
        a += this->asImp().evaluateIntegrand ( El, q ) * bfs.getWeight ( q );
      }
      res += a;
    }
    Dest += _config.getVolOfElement() * res;
  }

  //! interface function, has to be provided in derived classes.
  RealType evaluateIntegrand ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().evaluateIntegrand ( El, QuadPoint );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


//!===========================================================================================================================
//! Vector-Valued Intefaces
//!===========================================================================================================================


//! Integrator for \f$ (\int_\Omega w_i(x) \da )_{i} \f$
template <typename ConfiguratorType>
class QuocFEMassIntegrator {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType & _config;
  
public:
  QuocFEMassIntegrator( const ConfiguratorType & Conf ) : _config(Conf){ }

  void assembleAdd ( VectorType &Dest ) const {
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );
      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );
      RealType aux;
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.;   
        for ( int q = 0; q < numQuadPoints; ++q ) aux += bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
        Dest[ _config.localToGlobal ( El, dof ) ] += _config.getVolOfElement() * aux;
      }
    }
  }

};

//! Integrator for \f$(\int_\Omega s(x)  w_i(x) \da )_i\f$, of some scalar valued function \f$s\f$.
template <typename ConfiguratorType, typename Imp>
class QuocFENonlinOpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType & _config;
  
public:
  QuocFENonlinOpIntegrator( const ConfiguratorType & Conf ) : _config( Conf){ }

  virtual ~QuocFENonlinOpIntegrator( ) {}

  void assembleAdd ( VectorType &Dest, const RealType factor = 1.0 ) const {
    
    RealType *nl_cache = new RealType[ _config.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        nl_cache[q] = this->asImp().getNonlinearity ( El, q );

      RealType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.;   
        for ( int q = 0; q < numQuadPoints; ++q ) aux += nl_cache[q] * bfs.evaluate ( dof, q ) * bfs.getWeight ( q );
        Dest[ _config.localToGlobal ( El, dof ) ] += factor * _config.getVolOfElement() * aux;
      }
    }
    delete[] nl_cache;
  }

  //! interface function, has to be provided in derived classes.
  RealType getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! Integrator for \f$ (\int_\Omega A(x) \nabla w_i(x) \da )_{ij} \f$, of some matrix valued function \f$ A\f$.
template <typename ConfiguratorType, typename Imp>
class QuocFENonlinDiffOpIntegrator  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::ElementType ElementType;
  const ConfiguratorType & _config;
  
public:
  QuocFENonlinDiffOpIntegrator ( const ConfiguratorType & Conf ) : _config( Conf ){ }

  virtual ~QuocFENonlinDiffOpIntegrator( ) {}

  void assembleAdd ( VectorType &Dest ) const {

    DomVecType *nl_cache = new DomVecType[ _config.maxNumQuadPoints() ];

    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity (  El, q, nl_cache[q] );

      RealType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux = 0.0;   
        
        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.getWeight ( q ) * ( nl_cache[q].dot(bfs.evaluateGradient( dof, q ) ) );

        Dest[ _config.localToGlobal ( El, dof ) ] += _config.getVolOfElement() * aux;
      }
    }
    delete[] nl_cache;
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint,
                         DomVecType &NL ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint,  NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//!===========================================================================================================================
//! Multi-Vector-Valued Intefaces
//!===========================================================================================================================

//! Integrator for \f$ (\int_\Omega (v(x) \cdot w_i(x)) \da )_{i} \f$, of some vector valued function \f$v\f$. \f$v\f$  is supposed to be of dimension 3.
template <typename ConfiguratorType, typename Imp, int dimDomain = ConfiguratorType::dimDomain>
class QuocFENonlinVectorOpIntegrator {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType & _config;
  
public:
  QuocFENonlinVectorOpIntegrator ( const ConfiguratorType & Conf ) : _config( Conf){ }

  virtual ~QuocFENonlinVectorOpIntegrator( ) {}

  void assembleAdd ( typename ConfiguratorType::VectorType &Dest ) const {
    
    PointType *nl_cache = new PointType[ _config.maxNumQuadPoints() ];
    const int numGlobalDofs = _config.getNumGlobalDofs();
    
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( El, q, nl_cache[q] );

      PointType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux.setZero();    

        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.evaluate ( dof, q ) * bfs.getWeight ( q ) * nl_cache[q];

        for ( int comp = 0; comp < dimDomain; ++comp )
          Dest[ _config.localToGlobal ( El, dof ) + comp * numGlobalDofs ] += _config.getVolOfElement() * aux[comp];
      }
    }
    delete[] nl_cache;
  }

  void assembleDirichlet ( typename ConfiguratorType::VectorType &Dest, const MaskType& boundaryMask ) const {
      Dest.setZero();
      this->assembleAdd( Dest );
      const int numGlobalDofs = _config.getNumGlobalDofs();
      for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<dimDomain; ++comp ) Dest[i + comp * numGlobalDofs] = 0.0;
          }
      } 
  }
  
  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, PointType &NL ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};


//! Integrator for \f$ (\int_\Omega \tr( f(x)^T w_i(x)) \da )_{i} \f$, of some matrix valued function \f$f\f$.
template <typename ConfiguratorType, typename Imp, int dimDomain = ConfiguratorType::dimDomain>
class QuocFENonlinVectorDiffOpIntegrator {
    
  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType & _config;
  
public:
  QuocFENonlinVectorDiffOpIntegrator ( const ConfiguratorType & Conf ) : _config( Conf){ }

  virtual ~QuocFENonlinVectorDiffOpIntegrator( ) {}

  
  void assembleAdd ( typename ConfiguratorType::VectorType &Dest ) const {
    
    DerivativeVectorValuedType *nl_cache = new DerivativeVectorValuedType[ _config.maxNumQuadPoints() ];
    const int numGlobalDofs = _config.getNumGlobalDofs();
    
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q ) this->asImp().getNonlinearity ( El, q, nl_cache[q] );

      PointType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux.setZero();    

        for ( int q = 0; q < numQuadPoints; ++q ) aux +=  bfs.getWeight ( q ) * nl_cache[q] * bfs.evaluateGradient ( dof, q );

        for ( int comp = 0; comp < dimDomain; ++comp )
          Dest[ _config.localToGlobal ( El, dof ) + comp * numGlobalDofs ] += _config.getVolOfElement() * aux[comp];
      }
    }
    delete[] nl_cache;
  }
  
  void assembleDirichlet ( typename ConfiguratorType::VectorType &Dest, const MaskType& boundaryMask ) const {
      Dest.setZero();
      this->assembleAdd( Dest );
      const int numGlobalDofs = _config.getNumGlobalDofs();
      for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<dimDomain; ++comp ) Dest[i + comp * numGlobalDofs] = 0.0;
          }
      } 
  }

  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, DerivativeVectorValuedType &NL ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};










template <typename ConfiguratorType, typename Imp, int dimDomain = ConfiguratorType::dimDomain>
class QuocFENonlinVectorOpIntegratorWithRefCoord  {

  protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  const ConfiguratorType & _config;
  
public:
  QuocFENonlinVectorOpIntegratorWithRefCoord ( const ConfiguratorType & Conf ) : _config( Conf){ }

  virtual ~QuocFENonlinVectorOpIntegratorWithRefCoord( ) {}

  void assembleAdd ( typename ConfiguratorType::VectorType &Dest ) const {
    
    PointType *nl_cache = new PointType[ _config.maxNumQuadPoints() ];
    const int numGlobalDofs = _config.getNumGlobalDofs();
    
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
      const int numQuadPoints = bfs.numQuadPoints( );

      for ( int q = 0; q < numQuadPoints; ++q )
        this->asImp().getNonlinearity ( El, q, bfs.getRefCoord(q), nl_cache[q] );

      PointType aux;
  
      for ( int dof = 0; dof < numLocalDofs; ++dof ) {
        aux.setZero();    

        for ( int q = 0; q < numQuadPoints; ++q )
          aux += bfs.evaluate ( dof, q ) * bfs.getWeight ( q ) * nl_cache[q];

        for ( int comp = 0; comp < dimDomain; ++comp )
          Dest[ _config.localToGlobal ( El, dof ) + comp * numGlobalDofs ] += _config.getVolOfElement() * aux[comp];
      }
    }
    delete[] nl_cache;
  }

  void assembleDirichlet ( typename ConfiguratorType::VectorType &Dest, const MaskType& boundaryMask ) const {
      Dest.setZero();
      this->assembleAdd( Dest );
      const int numGlobalDofs = _config.getNumGlobalDofs();
      for( int i = 0; i < numGlobalDofs; ++i ){
          if ( boundaryMask[i] ){
            for( int comp=0; comp<3; ++comp ) Dest[i + comp * numGlobalDofs] = 0.0;
          }
      } 
  }
  
  //! interface function, has to be provided in derived classes.
  void getNonlinearity ( const typename ConfiguratorType::ElementType & El, int QuadPoint, PointType LocalCoords, PointType &NL ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    this->asImp().getNonlinearity (  El, QuadPoint, LocalCoords, NL );
  }

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }

};



//!===========================================================================================================================
//! Matrix-Valued Intefaces
//!===========================================================================================================================


//! General Interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp >
class QuocMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;

  explicit QuocMatrixValuedIntegratorBase ( const ConfiguratorType &conf ) : _config ( conf ) {}

public:
  void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
    tripletList.reserve(_config.getInitializer().getNumElements() * aol::Sqr( _config.getNumLocalDofs() ) );
    LocalMatrixType localMatrix;
    int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
    for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
      const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
      this->asImp().prepareLocalMatrix ( El, localMatrix );
      const int numLocalDofs = _config.getNumLocalDofs ( El );

      for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _config.localToGlobal ( El, i );

      for ( int i = 0; i < numLocalDofs; ++i ) {
        int glob_i = globalDofs[ i ];
        for ( int j = 0; j < numLocalDofs; ++j ) {
          int glob_j = globalDofs[ j ];
          tripletList.push_back( TripletType( glob_i, glob_j, _config.getVolOfElement() * Factor * localMatrix(i,j) ) );
        }
      }
    }
      
  }
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
      
    tripletListMasked.reserve(_config.getInitializer().getNumElements() * aol::Sqr( _config.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row()]) || (boundaryMask[tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
  }

  template <typename SparseMatrixType>
  void assemble ( SparseMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets( tripletList.cbegin(), tripletList.cend() );
  }
  
  template <typename SparseMatrixType>
  void assembleDirichlet ( SparseMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
       
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve(_config.getInitializer().getNumElements() * aol::Sqr( _config.getNumLocalDofs() ) );

    for( unsigned iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row()]) || (boundaryMask[tripletList[iter].col()]) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] )
         tripletListMasked.push_back( TripletType( i, i, 1.0 ) );
    }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
};



//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  w(x) \phi_i \phi_j dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class QuocFELinWeightedMassIntegrator :
      public QuocMatrixValuedIntegratorBase<  ConfiguratorType, QuocFELinWeightedMassIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  const ConfiguratorType &_config;
  
public:
  QuocFELinWeightedMassIntegrator ( const ConfiguratorType & Config ) : 
   QuocMatrixValuedIntegratorBase<  ConfiguratorType, QuocFELinWeightedMassIntegrator<ConfiguratorType, Imp> > ( Config ),
  _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, 
                                  LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _config.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          RealType b_i = bfs.evaluate( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          RealType b_j = bfs.evaluate( j, q );
          LocalMatrix(j,i) += nonlinearity * b_i * b_j * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};



//! The corresponding matrix assembly yields \f$ \left(\int_\Omega  A(x) \phi_i' \phi_j' dx\right)_{ij} \f$ for FE basis functions \f$ \phi_i,\phi_j \f$.
template <typename ConfiguratorType, typename Imp >
class QuocFELinWeightedStiffIntegrator :
      public QuocMatrixValuedIntegratorBase<  ConfiguratorType, QuocFELinWeightedStiffIntegrator<ConfiguratorType, Imp> > {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  const ConfiguratorType &_config;
  
public:
  QuocFELinWeightedStiffIntegrator ( const ConfiguratorType & Config ) : 
  QuocMatrixValuedIntegratorBase<  ConfiguratorType, QuocFELinWeightedStiffIntegrator<ConfiguratorType, Imp> > ( Config ), _config ( Config ) {}


  //! this function has to be provided in the implementation (derived class) of the interface
  inline RealType getNonlinearity ( const typename ConfiguratorType::ElementType &El, int QuadPoint ) const {
    return this->asImp().getNonlinearity ( El, QuadPoint );
  }

  //! this function computes the numerical quadrature of the bilinear form and saves the values locally.
  inline void prepareLocalMatrix ( const typename ConfiguratorType::ElementType &El, LocalMatrixType &LocalMatrix ) const {
    
    const int numDofs = _config.getNumLocalDofs ( El );	
    
    for ( int i = 0; i < numDofs; ++i )
      for ( int j = 0; j < numDofs; ++j )
        LocalMatrix(i,j) = 0.;

    RealType nonlinearity;

    const typename ConfiguratorType::BaseFuncSetType &bfs = _config.getBaseFunctionSet ( El );
    const int numQuadPoints = bfs.numQuadPoints( );

    for ( int q = 0; q < numQuadPoints; ++q ) {
      nonlinearity = getNonlinearity ( El, q );
      for ( int i = 0; i < numDofs; ++i ) {
          DomVecType grad_i = bfs.evaluateGradient( i, q );
        for ( int j = 0; j < numDofs; ++j ) {
          DomVecType grad_j = bfs.evaluateGradient( j, q );
          LocalMatrix(j,i) += nonlinearity * grad_i.dot(grad_j) * bfs.getWeight ( q );
        }
      }
    }

  }

protected:
  inline Imp &asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

};




//!===========================================================================================================================
//! BlockMatrix-Valued Intefaces
//!===========================================================================================================================

//! General Interface for matrix valued integrators
template < typename ConfiguratorType, typename Imp, int dimDomain = ConfiguratorType::dimDomain>
class QuocBlockMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::TripletType TripletType;

  explicit QuocBlockMatrixValuedIntegratorBase ( const ConfiguratorType &conf ): _config ( conf ) {}

protected:
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( dimDomain * dimDomain * aol::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer ().getNumElements ());
        LocalMatrixType localMatrix[dimDomain][dimDomain];
        const int numGlobalDofs = _config.getNumGlobalDofs();
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrix ( El, localMatrix );

            const int numLocalDofs = _config.getNumLocalDofs ( El );

            for ( int i = 0; i < numLocalDofs; ++i )
                globalDofs[ i ] = _config.localToGlobal ( El, i );
            
            for ( int argComp = 0; argComp < dimDomain; ++argComp )
                for ( int destComp = 0; destComp < dimDomain; ++destComp )
                 for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int j = 0; j < numLocalDofs; ++j ) {
                      int glob_j = globalDofs[ j ];
                      tripletList.push_back( TripletType( glob_i + destComp * numGlobalDofs, glob_j + argComp * numGlobalDofs, _config.getVolOfElement() * Factor * localMatrix[argComp][destComp]( i, j ) ) );
                    }
                }
        }
    }
  
public:

  template <typename BlockMatrixType>
  void assemble ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );
    Dest.setFromTriplets( tripletList.begin(), tripletList.end() ); 
  }
  
  template <typename BlockMatrixType>
  void assembleDirichlet ( BlockMatrixType &Dest, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _config.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( dimDomain * dimDomain * aol::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer().getNumElements() );
    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < dimDomain; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
    /*
     * Dest has additional size of dimDomain for constraints on \int u_i = 0
     * periodic Mask contains values 0 (no periodic node) and 1 (periodic node)
     * periodicIndicesMask contains for every periodic node the corresponding node on the boundary
     */
  template <typename BlockMatrixType>
  void assemblePeriodic ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _config.getNumGlobalDofs();
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( dimDomain * dimDomain * aol::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer().getNumElements ());

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
        const int colIndex = tripletList[iter].col(); 
        const int colNodeIdx = colIndex % numGlobalDofs;
        const int rowIndex = tripletList[iter].row(); 
        const int rowNodeIdx = rowIndex % numGlobalDofs;
      if( (periodicMask[rowNodeIdx] ) || (periodicMask[colNodeIdx] ) ){
       //Periodic Boundary node! 
        const int colZ = colIndex / numGlobalDofs;
        const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
        const int rowZ = rowIndex / numGlobalDofs;
        const int rowIndexPeriodic = periodicIndicesMask[rowNodeIdx] + rowZ * numGlobalDofs;
        if( (periodicMask[rowNodeIdx] ) && (periodicMask[colNodeIdx] ) ) tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndexPeriodic, tripletList[iter].value() ) );
        if( (!periodicMask[rowNodeIdx]) && (periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndex, colIndexPeriodic, tripletList[iter].value() ) );
        if( (periodicMask[rowNodeIdx]) && (!periodicMask[colNodeIdx]) )  tripletListMasked.push_back( TripletType( rowIndexPeriodic, colIndex, tripletList[iter].value() ) );
      }else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    //diagonal
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( periodicMask[i] ){
        for ( int Comp = 0; Comp < dimDomain; ++Comp ) tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
    
    //int u_i = 0
    typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
    QuocFEMassIntegrator<ConfiguratorType> ( _config ).assembleAdd( constraintVec );
    //colabse periodically
    for ( int nodeIdx=0; nodeIdx<numGlobalDofs; nodeIdx++ ) {
        if( periodicMask[nodeIdx] ){
            constraintVec[periodicIndicesMask[nodeIdx]] += constraintVec[nodeIdx];
            constraintVec[nodeIdx] = 0.0;
        }
    }
    //insert into matrix
    for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
        for( int comp=0; comp<dimDomain; ++comp){
            tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  dimDomain * numGlobalDofs + comp,    constraintVec[nodeIdx] ) );
            tripletListMasked.push_back( TripletType( dimDomain * numGlobalDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
        }
    }
    
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  void assembleTripletListDirichlet ( std::vector<TripletType> & tripletListMasked, const MaskType& boundaryMask, const RealType Factor = 1.0 ) const {
    std::vector<TripletType> tripletList;
    assembleTripletList ( tripletList, Factor );

    // Boundary Mask
    const int numGlobalDofs = _config.getNumGlobalDofs();
    tripletListMasked.reserve( dimDomain * dimDomain * aol::Sqr( _config.getNumLocalDofs() ) *_config.getInitializer ().getNumElements ());

    for( unsigned int iter=0; iter < tripletList.size(); ++iter ){
      if( (boundaryMask[tripletList[iter].row() % numGlobalDofs] ) || (boundaryMask[tripletList[iter].col() % numGlobalDofs] ) ){
       //Boundary node!        
      } else {
        tripletListMasked.push_back( tripletList[iter] );
      }
    }
    
    for ( int i = 0; i < _config.getNumGlobalDofs(); ++i ){
      if ( boundaryMask[i] ){
        for ( int Comp = 0; Comp < dimDomain; ++Comp )
            tripletListMasked.push_back( TripletType( i + Comp * numGlobalDofs, i + Comp * numGlobalDofs, 1.0 ) );
      }
    }
  }
  

protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
};




//!===========================================================================================================================
//! Boundary Integration Intefaces
//!===========================================================================================================================




//! Interface class for computing \f$ \int_{\partial \Omega} f \cdot \psi_i \, da \f$. 
//The function \f$ f \f$ has to be implemented
// \psi_i are the (vector-valued!) basis functions
template <typename ConfiguratorType, typename Imp>
class QuocBoundaryIntegrationInterface {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BoundaryElementType BoundaryElementType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::DomVecTypeBoundary DomVecTypeBoundary;
  typedef typename ConfiguratorType::BoundaryQuadType BoundaryQuadType;

  const ConfiguratorType &_config;
  
public:

  QuocBoundaryIntegrationInterface ( const ConfiguratorType &config ) : _config ( config ) {}

  void assembleAdd ( VectorType &dest ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs();
    const int numLocalDofs = _config.getNumLocalBoundaryDofs();
   
    for ( int boundaryElementIdx = 0; boundaryElementIdx < _config.getInitializer().getNumBoundaryElements(); ++boundaryElementIdx){
      const BoundaryElementType& bdrEl ( _config.getInitializer().getBoundaryElement( boundaryElementIdx ) );
      BoundaryQuadType quad;
      DomVecType refCoord;
      for ( int q = 0; q < _config.maxNumBoundaryQuadPoints(); q++ ) {
        DomVecTypeBoundary refCoordBoundary = quad.getRefCoord ( q );
        bdrEl.getRefCoord( refCoordBoundary, refCoord );
        DomVecType aux; asImp().getNonlinearity ( bdrEl, refCoord, aux );

        for ( int locNodeIndex = 0; locNodeIndex < numLocalDofs; locNodeIndex++ ) {
          int b = bdrEl.getNodeIndexOfElement( locNodeIndex );
          for ( int comp = 0; comp < ConfiguratorType::dimDomain; ++comp ){
            dest[ _config.localToGlobal ( bdrEl.getElement(), b ) + comp * numGlobalDofs ] += quad.getWeight ( q ) * bdrEl.getH() * _config.getBaseFunctionSet( ).evaluate ( b, refCoord ) * aux[comp];
          }
        }
      }
    } 
    
  }
  
  //! this method computes \f$ f \f$ and has to be implemented in the derived class
 void getNonlinearity ( const BoundaryElementType &bdrEl, const DomVecType refCoord, DomVecType &aux ) const {
    throw std::invalid_argument( aol::strprintf ( "Called the interface function. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
    return asImp().getNonlinearity ( bdrEl, refCoord, aux );
  }

protected:
  /** barton-nackman **/
  Imp& asImp( ) { return static_cast<Imp&> ( *this ); }
  const Imp& asImp( ) const { return static_cast<const Imp&> ( *this ); }

};


template <typename ConfiguratorType>
class IntegrateDispOverBoundary
  : public QuocBoundaryIntegrationInterface < ConfiguratorType, IntegrateDispOverBoundary<ConfiguratorType> > {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BoundaryElementType BoundaryElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  const std::vector<DomVecType> &_forceVec; // LEFT, _BOTTOM, (in 3d: _FRONT)
  
public:
  IntegrateDispOverBoundary( const ConfiguratorType &conf, const std::vector<DomVecType> &f  ) 
  : QuocBoundaryIntegrationInterface< ConfiguratorType, IntegrateDispOverBoundary< ConfiguratorType> > ( conf ), _forceVec ( f ) {}

  void getNonlinearity( const BoundaryElementType &bdrEl, const DomVecType & refCoord, DomVecType &aux ) const {
      aux.setZero();
      switch( bdrEl.getBoundaryType() ){
             case LEFT:  { aux = _forceVec[0];}break;
             case BOTTOM:{ aux = _forceVec[1];}break;
             case RIGHT: { aux = _forceVec[2];}break;
             case TOP:   { aux = _forceVec[3];}break;
             case FRONT: { aux = _forceVec[4];}break;
             case BACK:  { aux = _forceVec[5];}break;
             default:  break;
         }
  }
};


}//end namespace


#endif
