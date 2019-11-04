#ifndef __QUOCAFFINESYMGRADINTEGRATORS_H
#define __QUOCAFFINESYMGRADINTEGRATORS_H

#include "quocQuadrature.h"
#include "quocIntegrator.h"

namespace quocFE {

//! Interface class for computing \f$ \int_{\partial \Omega} f \cdot \psi_i \, da \f$. 
//The function \f$ f \f$ has to be implemented
// \psi_i are the (vector-valued!) affine basis functions:
// in 2d psi \in ( (x,0), (0,y), (y,x) )
// in 3d psi \in ( (x,0,0), (0,y,0), (0,0,z), (y,x,0), (z,0,x), (0,z,y) )
template <typename ConfiguratorType, typename Imp>
class QuocAffineBoundaryIntegrationInterface {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BoundaryElementType BoundaryElementType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::DomVecTypeBoundary DomVecTypeBoundary;
  typedef typename ConfiguratorType::BoundaryQuadType BoundaryQuadType;
  typedef typename ConfiguratorType::GlobalAffineSymGradBaseFuncSet GlobalAffineSymGradBaseFuncSet;
  
  const ConfiguratorType &_config;
  
public:

  QuocAffineBoundaryIntegrationInterface ( const ConfiguratorType &config ) : _config ( config ) {}

  void assembleAdd ( VectorType &dest ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs();
    const int numLocalDofs = _config.getNumLocalBoundaryDofs();
   
    GlobalAffineSymGradBaseFuncSet globAffBfs;
    
    for ( int boundaryElementIdx = 0; boundaryElementIdx < _config.getInitializer().getNumBoundaryElements(); ++boundaryElementIdx){
      const BoundaryElementType& bdrEl ( _config.getInitializer().getBoundaryElement( boundaryElementIdx ) );
      BoundaryQuadType quad;
      DomVecType refCoord, globalCoord;
      for ( int q = 0; q < _config.maxNumBoundaryQuadPoints(); q++ ) {
        DomVecTypeBoundary refCoordBoundary = quad.getRefCoord ( q );
        bdrEl.getRefCoord( refCoordBoundary, refCoord );
        DomVecType aux; asImp().getNonlinearity ( bdrEl, refCoord, aux );

        _config.getGlobalCoordsForBoundaryElement ( bdrEl, refCoord, globalCoord );

        for( int affineIndex = 0; affineIndex < globAffBfs.numBaseFuncs; affineIndex++ ){
          DomVecType affinePart; globAffBfs.evaluate ( affineIndex, globalCoord, affinePart );
          dest[affineIndex] += quad.getWeight ( q ) * bdrEl.getH() * aux.dot( affinePart  );
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
class IntegrateAffineDispOverBoundary
  : public QuocAffineBoundaryIntegrationInterface < ConfiguratorType, IntegrateAffineDispOverBoundary<ConfiguratorType> > {

  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::BoundaryElementType BoundaryElementType;
  typedef typename ConfiguratorType::DomVecType DomVecType;

  const std::vector<DomVecType> &_forceVec; // LEFT, _BOTTOM, (in 3d: _FRONT)
  
public:
  IntegrateAffineDispOverBoundary( const ConfiguratorType &conf, const std::vector<DomVecType> &f  ) 
  : QuocAffineBoundaryIntegrationInterface< ConfiguratorType, IntegrateAffineDispOverBoundary< ConfiguratorType> > ( conf ), _forceVec ( f ) {}

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





template < typename ConfiguratorType, typename Imp, int dimDomain = ConfiguratorType::dimDomain>
class QuocPlusAffineSymGradBlockMatrixValuedIntegratorBase {
public:
  typedef typename ConfiguratorType::RealType   RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::ElementType ElementType;
  typedef typename ConfiguratorType::MaskType MaskType;
  typedef typename ConfiguratorType::TripletType TripletType;
  typedef typename ConfiguratorType::GlobalAffineSymGradBaseFuncSet GlobalAffineSymGradBaseFuncSet;
  typedef typename ConfiguratorType::LocalMatrixType LocalMatrixType;
  typedef typename ConfiguratorType::LocalMatrixTypeAffineSymGrad LocalMatrixTypeAffineSymGrad;
  typedef typename ConfiguratorType::LocalMatrixTypeMixed LocalMatrixTypeMixed;

  explicit QuocPlusAffineSymGradBlockMatrixValuedIntegratorBase ( const ConfiguratorType &conf ): _config ( conf ) {}

public:
    
    void assembleTripletList ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        tripletList.reserve( aol::Sqr( dimDomain * _config.getNumLocalDofs() ) * _config.getInitializer ().getNumElements ());
        LocalMatrixType localMatrix[dimDomain][dimDomain];
        const int numGlobalDofs = _config.getNumGlobalDofs();
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrix ( El, localMatrix );
            const int numLocalDofs = _config.getNumLocalDofs ( El );

            for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _config.localToGlobal ( El, i );
            
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
    
  //assembles only dimAffine x dimQuocVector Matrix (without offsets)
   void assembleTripletListMixedPart ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
       const int numGlobalDofs = _config.getNumGlobalDofs();
       const int numDofsOtherDisp = _config.getNumGlobalDofs() * _config.dimDomain;
        GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
        tripletList.reserve( dimDomain * _config.getNumLocalDofs() * numAffineDofs  *_config.getInitializer ().getNumElements ());
        int globalDofs[ ConfiguratorType::maxNumLocalDofs ];
        LocalMatrixTypeMixed localMatrix[dimDomain];
        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
            // assemble the local matrix for the current element
            this->asImp().prepareLocalMatrixMixedPart ( El, localMatrix );
            const int numLocalDofs = _config.getNumLocalDofs ( El );
            for ( int i = 0; i < numLocalDofs; ++i ) globalDofs[ i ] = _config.localToGlobal ( El, i );
            for ( int comp = 0; comp < dimDomain; ++comp )
                 for ( int i = 0; i < numLocalDofs; ++i ) {
                    int glob_i = globalDofs[ i ];
                    for ( int affIndex = 0; affIndex < globAffBfs.numBaseFuncs; ++affIndex ) {
                      tripletList.push_back( TripletType( affIndex, glob_i + comp * numGlobalDofs, _config.getVolOfElement() * Factor * localMatrix[comp]( affIndex, i ) ) );
                    }
                }
        }
    }
    
    //without offset
    void assembleTripletListAffinePart ( std::vector<TripletType> & tripletList, const RealType Factor ) const {
        const int numDofsOtherDisp = _config.getNumGlobalDofs() * _config.dimDomain;
        GlobalAffineSymGradBaseFuncSet globAffBfs;
        LocalMatrixTypeAffineSymGrad localMatrix;
        for ( int elementIdx = 0; elementIdx < _config.getInitializer().getNumElements(); ++elementIdx){
            const ElementType& El ( _config.getInitializer().getElement( elementIdx ) );
            this->asImp().prepareLocalMatrixAffinePart ( El, localMatrix );
            for( int affIndexArg = 0; affIndexArg < globAffBfs.numBaseFuncs; ++affIndexArg )
                for( int affIndexDest = 0; affIndexDest < globAffBfs.numBaseFuncs; ++affIndexDest ){
                   tripletList.push_back( TripletType( affIndexDest, affIndexArg, _config.getVolOfElement() * Factor * localMatrix( affIndexArg, affIndexDest ) ) );
                }
        }
    }
  
public:
  
    
    
  template <typename BlockMatrixType>
  void assemblePeriodic ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs(); const int numPeriodicDispDofs = dimDomain * numGlobalDofs;
    const int numLocalDofs = _config.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( aol::Sqr( numLocalDofs * dimDomain ) + 2 * numLocalDofs * dimDomain * numAffineDofs + aol::Sqr( numAffineDofs ) ) * _config.getInitializer().getNumElements() + dimDomain * numPeriodicDispDofs );
   
 
    //periodic part
        std::vector<TripletType> tripletList;
        assembleTripletList ( tripletList, Factor );

        // Boundary Mask
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
    
    
    //int u_i^per = 0
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
                tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numPeriodicDispDofs + comp,   constraintVec[nodeIdx] ) );
                tripletListMasked.push_back( TripletType( numPeriodicDispDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
            }
        }
    
//     cout << "size = " << tripletListMasked.size() << endl;
//     cout << "max StorageIndex = " << std::numeric_limits<int>::max() << endl;
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
    
    
    
    
      /*
     * Dest has additional size of dimDomain for constraints on \int u_i = 0
     * periodic Mask contains values 0 (no periodic node) and 1 (periodic node)
     * periodicIndicesMask contains for every periodic node the corresponding node on the boundary
     */
  template <typename BlockMatrixType>
  void assemblePeriodicPlusAffine ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs(); const int numPeriodicDispDofs = dimDomain * numGlobalDofs;
    const int numLocalDofs = _config.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( aol::Sqr( numLocalDofs * dimDomain ) + 2 * numLocalDofs * dimDomain * numAffineDofs + aol::Sqr( numAffineDofs ) ) * _config.getInitializer().getNumElements() + dimDomain * numPeriodicDispDofs );
 
    //periodic part
        std::vector<TripletType> tripletList;
        assembleTripletList ( tripletList, Factor );

        // Boundary Mask
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
    
    
    //Mixed Part
       std::vector<TripletType> tripletListMixed;
       assembleTripletListMixedPart( tripletListMixed, Factor );
       // Boundary Mask
       for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
            const int colIndex = tripletListMixed[iter].col(); 
            const int colNodeIdx = colIndex % numGlobalDofs;
            const int rowIndex = tripletListMixed[iter].row();
            if( periodicMask[colNodeIdx] ){
            //Periodic Boundary node! 
                const int colZ = colIndex / numGlobalDofs;
                const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
                tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndexPeriodic, tripletListMixed[iter].value() ) );
                tripletListMasked.push_back( TripletType( colIndexPeriodic, rowIndex + numPeriodicDispDofs, tripletListMixed[iter].value() ) );
            }else {
                tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndex, tripletListMixed[iter].value() ) );
                tripletListMasked.push_back( TripletType( colIndex, rowIndex + numPeriodicDispDofs, tripletListMixed[iter].value() ) );
            }
       }
    
    //Affine Part
       std::vector<TripletType> tripletListAffine;
       assembleTripletListAffinePart( tripletListAffine, Factor );
       for( unsigned int iter=0; iter < tripletListAffine.size(); ++iter ){
            const int colIndex = tripletListAffine[iter].col(); 
            const int rowIndex = tripletListAffine[iter].row();
            tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndex + numPeriodicDispDofs, tripletListAffine[iter].value() ) );
       }
    
    
    //int u_i^per = 0
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
                tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numPeriodicDispDofs + numAffineDofs  + comp,   constraintVec[nodeIdx] ) );
                tripletListMasked.push_back( TripletType( numPeriodicDispDofs + numAffineDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
            }
        }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  
  
  template <typename BlockMatrixType>
  void assemblePeriodicPlusAffineSeparated ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs(); const int numPeriodicDispDofs = dimDomain * numGlobalDofs;
    const int numLocalDofs = _config.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( aol::Sqr( numLocalDofs * dimDomain ) + 2 * numLocalDofs * dimDomain * numAffineDofs + aol::Sqr( numAffineDofs ) ) * _config.getInitializer().getNumElements() + dimDomain * numPeriodicDispDofs );
 
    //periodic part
        std::vector<TripletType> tripletList;
        assembleTripletList ( tripletList, Factor );

        // Boundary Mask
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
    
    
    //Mixed Part
       std::vector<TripletType> tripletListMixed;
       assembleTripletListMixedPart( tripletListMixed, Factor );
       // Boundary Mask
       for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
            const int colIndex = tripletListMixed[iter].col(); 
            const int colNodeIdx = colIndex % numGlobalDofs;
            const int rowIndex = tripletListMixed[iter].row();
            if( periodicMask[colNodeIdx] ){
            //Periodic Boundary node! 
                const int colZ = colIndex / numGlobalDofs;
                const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
                //tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndexPeriodic, tripletListMixed[iter].value() ) );
                tripletListMasked.push_back( TripletType( colIndexPeriodic, rowIndex + numPeriodicDispDofs, tripletListMixed[iter].value() ) );
            }else {
                //tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndex, tripletListMixed[iter].value() ) );
                tripletListMasked.push_back( TripletType( colIndex, rowIndex + numPeriodicDispDofs, tripletListMixed[iter].value() ) );
            }
       }
    
    //Affine Part
       std::vector<TripletType> tripletListAffine;
       assembleTripletListAffinePart( tripletListAffine, Factor );
       for( unsigned int iter=0; iter < tripletListAffine.size(); ++iter ){
            const int colIndex = tripletListAffine[iter].col(); 
            const int rowIndex = tripletListAffine[iter].row();
            tripletListMasked.push_back( TripletType( rowIndex + numPeriodicDispDofs, colIndex + numPeriodicDispDofs, tripletListAffine[iter].value() ) );
       }
    
    
    //int u_i^per = 0
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
                tripletListMasked.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numPeriodicDispDofs + numAffineDofs  + comp,   constraintVec[nodeIdx] ) );
                tripletListMasked.push_back( TripletType( numPeriodicDispDofs + numAffineDofs + comp,  nodeIdx + comp * numGlobalDofs,    constraintVec[nodeIdx] ) );
            }
        }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  template <typename BlockMatrixType>
  void assembleMixedPeriodic ( BlockMatrixType &Dest, const MaskType& periodicMask, const std::vector<int> & periodicIndicesMask, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs(); const int numPeriodicDispDofs = dimDomain * numGlobalDofs;
    const int numLocalDofs = _config.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    std::vector<TripletType> tripletListMasked;
    tripletListMasked.reserve( ( aol::Sqr( numLocalDofs * dimDomain ) + 2 * numLocalDofs * dimDomain * numAffineDofs + aol::Sqr( numAffineDofs ) ) * _config.getInitializer().getNumElements() + dimDomain * numPeriodicDispDofs );
    
    //Mixed Part
       std::vector<TripletType> tripletListMixed;
       assembleTripletListMixedPart( tripletListMixed, Factor );
       // Boundary Mask
       for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
            const int colIndex = tripletListMixed[iter].col(); 
            const int colNodeIdx = colIndex % numGlobalDofs;
            const int rowIndex = tripletListMixed[iter].row();
            if( periodicMask[colNodeIdx] ){
            //Periodic Boundary node! 
                const int colZ = colIndex / numGlobalDofs;
                const int colIndexPeriodic = periodicIndicesMask[colNodeIdx] + colZ * numGlobalDofs;
                tripletListMasked.push_back( TripletType( colIndexPeriodic, rowIndex, tripletListMixed[iter].value() ) );
            }else {
                tripletListMasked.push_back( TripletType( colIndex, rowIndex, tripletListMixed[iter].value() ) );
            }
       }
    
    Dest.setFromTriplets( tripletListMasked.begin(), tripletListMasked.end() );
  }
  
  
  
  template <typename BlockMatrixType>
  void assembleNeumann ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs(); const int numNeumannDispDofs = dimDomain * numGlobalDofs;
    const int numLocalDofs = _config.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    
    std::vector<TripletType> tripletList;
    tripletList.reserve( ( aol::Sqr( numLocalDofs * dimDomain ) + 2 * numLocalDofs * dimDomain * numAffineDofs + aol::Sqr( numAffineDofs ) ) * _config.getInitializer().getNumElements() + dimDomain * numNeumannDispDofs );
    
    //second derivative of energy part
    assembleTripletList ( tripletList, Factor );
    
    //constraint : int u_i^Neumann = 0
        typename ConfiguratorType::VectorType constraintVec ( numGlobalDofs ); constraintVec.setZero();
        QuocFEMassIntegrator<ConfiguratorType> ( _config ).assembleAdd( constraintVec );
        //insert into matrix
        for( int nodeIdx=0; nodeIdx < numGlobalDofs; ++nodeIdx ){
            for( int comp=0; comp<dimDomain; ++comp){
                tripletList.push_back( TripletType( nodeIdx + comp * numGlobalDofs,  numNeumannDispDofs + comp,   constraintVec[nodeIdx] ) );
                tripletList.push_back( TripletType( numNeumannDispDofs + comp,  nodeIdx + comp * numGlobalDofs,   constraintVec[nodeIdx] ) );
            }
        }
    Dest.setFromTriplets( tripletList.begin(), tripletList.end() );
  }
  
  
  template <typename BlockMatrixType>
  void assembleMixedNeumann ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    const int numGlobalDofs = _config.getNumGlobalDofs(); const int numNeumannDispDofs = dimDomain * numGlobalDofs;
    const int numLocalDofs = _config.getNumLocalDofs();
    GlobalAffineSymGradBaseFuncSet globAffBfs; const int numAffineDofs = globAffBfs.numBaseFuncs;
    
    //Mixed Part
    std::vector<TripletType> tripletListMixed;
    tripletListMixed.reserve( ( aol::Sqr( numLocalDofs * dimDomain ) + 2 * numLocalDofs * dimDomain * numAffineDofs + aol::Sqr( numAffineDofs ) ) * _config.getInitializer().getNumElements() + dimDomain * numNeumannDispDofs );
    assembleTripletListMixedPart( tripletListMixed, Factor );
    
    // transpose 
    std::vector<TripletType> tripletListTransposed;
    tripletListTransposed.reserve( tripletListMixed.size() );
    for( unsigned int iter=0; iter < tripletListMixed.size(); ++iter ){
        const int colIndex = tripletListMixed[iter].col(); 
        const int colNodeIdx = colIndex % numGlobalDofs;
        const int rowIndex = tripletListMixed[iter].row();
        tripletListTransposed.push_back( TripletType( colIndex, rowIndex, tripletListMixed[iter].value() ) );
    }
    
    Dest.setFromTriplets( tripletListTransposed.begin(), tripletListTransposed.end() );
  }
  
  
  
 template <typename BlockMatrixType>
  void assembleAffine ( BlockMatrixType &Dest, const RealType Factor = 1.0 ) const {
    //Affine Part
    std::vector<TripletType> tripletListAffine;
    assembleTripletListAffinePart( tripletListAffine, Factor );
    Dest.setFromTriplets( tripletListAffine.begin(), tripletListAffine.end() );
  }
  
protected:
  // barton-nackman
  inline Imp& asImp() { return static_cast<Imp&> ( *this ); }
  inline const Imp& asImp() const { return static_cast<const Imp&> ( *this ); }
  
  const ConfiguratorType &_config;
  
};



} //end namespace

#endif
