#ifndef __QUOCMBASEFUNCTIONSETS_H
#define __QUOCBASEFUNCTIONSETS_H

#include <general.h>

namespace quocFE {

//! Inteface
template <typename DataTypeContainer, class QuadRuleType >
class QuocBaseFunctionSetInterface  {
protected:
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
public:
  QuocBaseFunctionSetInterface( ) {}

  int numQuadPoints( ) const { return QuadRuleType::numQuadPoints;}
  inline RealType getWeight ( int QuadPoint ) const { return _quadRule.getWeight ( QuadPoint );}
  inline const DomVecType& getRefCoord ( int QuadPoint ) const { return _quadRule.getRefCoord ( QuadPoint );}

protected:
  QuadRuleType _quadRule;
};



// template <typename RealType, int Dim, class QuadRuleType>
// class QuocBaseFunctionSet {
// };

template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocBaseFunctionSet1D : public QuocBaseFunctionSetInterface< DataTypeContainer, QuadType >  {

  typedef typename DataTypeContainer::RealType RealType;

  static RealType _b1   ( const RealType &c ) { return ( 1. - c ); }
  static RealType _b2   ( const RealType &c ) { return c[0]; }
  
  static RealType _dx_b1   ( const RealType &/*c*/ ) { return - 1.; }
  static RealType _dx_b2   ( const RealType &/*c*/ ) { return 1.; }
  

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealType &RefCoord );
  BASIS_FUNC_TYPE _basis[2];
  BASIS_FUNC_TYPE _deriv_basis[2];
  const RealType _hx;
  
public:
  QuocBaseFunctionSet1D( const RealType hx ) : _hx ( hx ){
    _basis[0] = _b1;
    _basis[1] = _b2;

    _deriv_basis[0] = _dx_b1;
    _deriv_basis[1] = _dx_b2;
  }
  
  enum { numBaseFuncs = 2 };

  void setElement ( const ElementType &/*El*/ ) { }

  RealType evaluate ( int BaseFuncNum, const RealType &RefCoord ) const { return _basis[BaseFuncNum] ( RefCoord );}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  void evaluateGradient ( int BaseFuncNum, const RealType &RefCoord, RealType &Gradient ) const {Gradient = _deriv_basis[BaseFuncNum] ( RefCoord ) / _hx;}
  inline RealType evaluateGradient ( int BaseFuncNum, int QuadPoint ) const { return RealType( _deriv_basis[BaseFuncNum] ( this->_quadRule.getRefCoord( QuadPoint ) ) / _hx );}
};



//! The basefunctionset for bilinear elements in 2d
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocBaseFunctionSet2D : public QuocBaseFunctionSetInterface< DataTypeContainer, QuadType >{
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
  // b1 (0,0)
  // b2 (1,0)
  // b3 (0,1)
  // b4 (1,1)
  
  static RealType _b1   ( const DomVecType &c ) { return ( 1.0 - c[0] ) * ( 1.0 - c[1] ); }
  static RealType _b2   ( const DomVecType &c ) { return c[0] * ( 1.0 - c[1] ); }
  static RealType _b3   ( const DomVecType &c ) { return ( 1.0 - c[0] ) * c[1]; }
  static RealType _b4   ( const DomVecType &c ) { return c[0]*c[1]; }
  
  static RealType _dx_b1   ( const DomVecType &c ) { return c[1] - 1.0; }
  static RealType _dx_b2   ( const DomVecType &c ) { return 1.0 - c[1]; }
  static RealType _dx_b3   ( const DomVecType &c ) { return -c[1]; }
  static RealType _dx_b4   ( const DomVecType &c ) { return c[1]; }

  static RealType _dy_b1   ( const DomVecType &c ) { return c[0] - 1.0; }
  static RealType _dy_b2   ( const DomVecType &c ) { return -c[0]; }
  static RealType _dy_b3   ( const DomVecType &c ) { return 1.0 - c[0]; }
  static RealType _dy_b4   ( const DomVecType &c ) { return c[0]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &c );
  BASIS_FUNC_TYPE _basis[4];
  BASIS_FUNC_TYPE _deriv_basis[2][4];
  const RealType _hx, _hy;

public:
  QuocBaseFunctionSet2D ( const RealType hx, const RealType hy ) : _hx ( hx ), _hy ( hy ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
  }

  enum { numBaseFuncs = 4 };
  
  RealType evaluate ( int BaseFuncNum, const DomVecType &c ) const { return _basis[BaseFuncNum] ( c );}
  inline const RealType evaluate( int BaseFuncNum, int QuadPoint ) const { return evaluate( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );}
  
  void evaluateGradient ( int BaseFuncNum, const DomVecType &c, DomVecType &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
  }

  inline DomVecType evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return DomVecType( _deriv_basis[0][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hx, 
                       _deriv_basis[1][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hy );
  }

};


//! The basefunctionset for trilinear elements in 3d.
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocBaseFunctionSet3D : public QuocBaseFunctionSetInterface< DataTypeContainer, QuadType >{

  typedef typename DataTypeContainer::RealType RealType;  
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
  static inline RealType _b1_1d ( RealType x ) { return 1. - x; }
  static inline RealType _b2_1d ( RealType x ) { return x; }

  static inline RealType _d_b1_1d ( RealType ) { return -1.; }
  static inline RealType _d_b2_1d ( RealType ) { return 1.; }

  // b1 (0,0,0)
  // b2 (1,0,0)
  // b3 (0,1,0)
  // b4 (1,1,0)
  // b5 (0,0,1)
  // b6 (1,0,1)
  // b7 (0,1,1)
  // b8 (1,1,1)
  
  static RealType _b1   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b2   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b3   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b4   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b5   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b6   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b7   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b8   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }

  static RealType _dx_b1 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b2 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b3 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b4 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b5 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b6 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b7 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b8 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dy_b1 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b2 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b3 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b4 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b5 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b6 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b7 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b8 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dz_b1 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b2 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b3 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b4 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b5 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b6 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b7 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b8 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &c );
  BASIS_FUNC_TYPE _deriv_basis[3][8];
  BASIS_FUNC_TYPE _basis[8];

  const RealType _hx, _hy, _hz;

public:
  QuocBaseFunctionSet3D( const RealType hx, const RealType hy, const RealType hz ) : _hx ( hx ), _hy ( hy ), _hz (hz) {
    
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;
    _deriv_basis[2][6] = _dz_b7;
    _deriv_basis[2][7] = _dz_b8;

//     this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 8 };
  
  RealType evaluate ( int BaseFuncNum, const DomVecType &c ) const {
    return _basis[BaseFuncNum] ( c );
  }
  
  inline const RealType evaluate( int BaseFuncNum, int QuadPoint ) const {
    return evaluate( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }

//   inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
//     return BaseFunctionSetInterface<RealType, DomVecType, DomVecType, 8, QuadRuleType, QuocBaseFunctionSet2D<RealType, qc::QC_3D, QuadRuleType> >::evaluate ( BaseFuncNum, QuadPoint );
//   }
//   
  void evaluateGradient ( int BaseFuncNum, const DomVecType &c, DomVecType &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( c ) / _hz;
  }

  inline DomVecType evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return DomVecType( _deriv_basis[0][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hx, 
                       _deriv_basis[1][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hy,
                       _deriv_basis[2][BaseFuncNum]( this->_quadRule.getRefCoord( QuadPoint ) ) / _hz
                     );
  }

};



//! ====================================================================================================================================
//  ====================================================================================================================================
//! ====================================================================================================================================
//                           Cached Intefaces
//  =====================================================================================================================================
//! ====================================================================================================================================
//  =====================================================================================================================================


template <typename DataTypeContainer, class QuadRuleType, int numBaseFuncs, typename Imp >
class QuocCachedBaseFunctionSetInterface  {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
public:
  QuocCachedBaseFunctionSetInterface( ) {}

  int numQuadPoints( ) const { return QuadRuleType::numQuadPoints;}
  inline RealType getWeight ( int QuadPoint ) const {return _quadRule.getWeight ( QuadPoint );}
  inline const DomVecType& getRefCoord ( int QuadPoint ) const { return _quadRule.getRefCoord ( QuadPoint );}

  //! read the cached value of the basis function with number BaseFuncNum at the given quadrature point
  inline RealType evaluate ( int BaseFuncNum, int QuadPoint ) const { return basisQuadValues[BaseFuncNum][QuadPoint];}
  inline const DomVecType& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const { return basisQuadGradients[BaseFuncNum][QuadPoint];}
  inline RealType evaluate ( int BaseFuncNum, const DomVecType &RefCoord ) const { return asImp().evaluate ( BaseFuncNum, RefCoord );}
  inline void evaluateGradient ( int BaseFuncNum, const DomVecType& RefCoord, DomVecType& Gradient ) const {asImp().evaluateGradient ( BaseFuncNum, RefCoord, Gradient );}

protected:
  Imp &asImp() { return static_cast<Imp&> ( *this ); }

  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  void initializeQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ )
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        basisQuadValues[b][i] = evaluate ( b,  _quadRule.getRefCoord ( i ) );
        evaluateGradient ( b, _quadRule.getRefCoord ( i ), basisQuadGradients[b][i] );
      }
  }

  /**** cache the values of the basis functions at the quadrature points ****/
  RealType     basisQuadValues   [numBaseFuncs][QuadRuleType::numQuadPoints];
  DomVecType   basisQuadGradients[numBaseFuncs][QuadRuleType::numQuadPoints];
  QuadRuleType _quadRule;
};



template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocCachedBaseFunctionSet1D : public QuocCachedBaseFunctionSetInterface< DataTypeContainer, QuadType, 2, QuocCachedBaseFunctionSet1D<DataTypeContainer,QuadType,ElementType> >  {

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;

  static RealType _b1   ( const RealType &c ) { return ( 1. - c ); }
  static RealType _b2   ( const RealType &c ) { return c; }
  
  static RealType _dx_b1   ( const RealType &/*c*/ ) { return - 1.; }
  static RealType _dx_b2   ( const RealType &/*c*/ ) { return 1.; }
  

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealType &RefCoord );
  BASIS_FUNC_TYPE _basis[2];
  BASIS_FUNC_TYPE _deriv_basis[2];
  const RealType _hx;
  
public:
  QuocCachedBaseFunctionSet1D( const RealType hx ) : _hx ( hx ){
    _basis[0] = _b1;
    _basis[1] = _b2;

    _deriv_basis[0] = _dx_b1;
    _deriv_basis[1] = _dx_b2;
    
    this->initializeQuadCache( );
  }
  
  enum { numBaseFuncs = 2 };

  void setElement ( const ElementType &/*El*/ ) { }
  
//   RealType evaluate ( int BaseFuncNum, const RealType &c ) const { return _basis[BaseFuncNum] ( c );}
  RealType evaluate ( int BaseFuncNum, const DomVecType &c ) const { return _basis[BaseFuncNum] ( c[0] );}
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return QuocCachedBaseFunctionSetInterface<DataTypeContainer, QuadType, 2, QuocCachedBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> >::evaluate ( BaseFuncNum, QuadPoint );
  }

  void evaluateGradient ( int BaseFuncNum, const DomVecType &RefCoord, DomVecType &Gradient ) const {Gradient[0] = _deriv_basis[BaseFuncNum] ( RefCoord[0] ) / _hx;}

  inline const DomVecType & evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return QuocCachedBaseFunctionSetInterface<DataTypeContainer, QuadType, 2, QuocCachedBaseFunctionSet1D<DataTypeContainer, QuadType, ElementType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }
  
};


//! The basefunctionset for bilinear elements in 2d
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocCachedBaseFunctionSet2D 
  : public QuocCachedBaseFunctionSetInterface< DataTypeContainer, QuadType, 4, QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> >{
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
  static RealType _b1   ( const DomVecType &c ) { return ( 1.0 - c[0] ) * ( 1.0 - c[1] ); }
  static RealType _b2   ( const DomVecType &c ) { return c[0] * ( 1.0 - c[1] ); }
  static RealType _b3   ( const DomVecType &c ) { return ( 1.0 - c[0] ) * c[1]; }
  static RealType _b4   ( const DomVecType &c ) { return c[0]*c[1]; }
  
  static RealType _dx_b1   ( const DomVecType &c ) { return c[1] - 1.0; }
  static RealType _dx_b2   ( const DomVecType &c ) { return 1.0 - c[1]; }
  static RealType _dx_b3   ( const DomVecType &c ) { return -c[1]; }
  static RealType _dx_b4   ( const DomVecType &c ) { return c[1]; }

  static RealType _dy_b1   ( const DomVecType &c ) { return c[0] - 1.0; }
  static RealType _dy_b2   ( const DomVecType &c ) { return -c[0]; }
  static RealType _dy_b3   ( const DomVecType &c ) { return 1.0 - c[0]; }
  static RealType _dy_b4   ( const DomVecType &c ) { return c[0]; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &c );
  BASIS_FUNC_TYPE _basis[4];
  BASIS_FUNC_TYPE _deriv_basis[2][4];
  const RealType  _hx, _hy;

public:
  QuocCachedBaseFunctionSet2D ( const RealType hx, const RealType hy ) : _hx ( hx ), _hy ( hy ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    
    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 4 };
  
  RealType evaluate ( int BaseFuncNum, const DomVecType &c ) const { return _basis[BaseFuncNum] ( c );}
  
  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return QuocCachedBaseFunctionSetInterface<DataTypeContainer, QuadType, 4, QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateGradient ( int BaseFuncNum, const DomVecType &c, DomVecType &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
  }
  
  inline const DomVecType& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return QuocCachedBaseFunctionSetInterface<DataTypeContainer, QuadType, 4, QuocCachedBaseFunctionSet2D<DataTypeContainer, QuadType, ElementType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

};


//! The basefunctionset for trilinear elements in 3d.
template <typename DataTypeContainer, typename QuadType, typename ElementType>
class QuocCachedBaseFunctionSet3D 
  : public QuocCachedBaseFunctionSetInterface< DataTypeContainer, QuadType, 8, QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> >{

  typedef typename DataTypeContainer::RealType RealType;  
  typedef typename DataTypeContainer::DomVecType DomVecType;
    
  static inline RealType _b1_1d ( RealType x ) { return 1. - x; }
  static inline RealType _b2_1d ( RealType x ) { return x; }

  static inline RealType _d_b1_1d ( RealType ) { return -1.; }
  static inline RealType _d_b2_1d ( RealType ) { return 1.; }

  // b1 (0,0,0)
  // b2 (1,0,0)
  // b3 (0,1,0)
  // b4 (1,1,0)
  // b5 (0,0,1)
  // b6 (1,0,1)
  // b7 (0,1,1)
  // b8 (1,1,1)
  
  static RealType _b1   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b2   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b3   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b4   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b1_1d ( c[2] ); }
  static RealType _b5   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b6   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b1_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b7   ( const DomVecType &c ) { return _b1_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }
  static RealType _b8   ( const DomVecType &c ) { return _b2_1d ( c[0] ) * _b2_1d ( c[1] ) * _b2_1d ( c[2] ); }

  static RealType _dx_b1 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b2 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b3 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b4 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dx_b5 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b6 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b7 ( const DomVecType &c ) { return _d_b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dx_b8 ( const DomVecType &c ) { return _d_b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dy_b1 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b2 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b3 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b4 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b1_1d ( c[2] ); }
  static RealType _dy_b5 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b6 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b1_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b7 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }
  static RealType _dy_b8 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_d_b2_1d ( c[1] ) *_b2_1d ( c[2] ); }

  static RealType _dz_b1 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b2 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b3 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b4 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b1_1d ( c[2] ); }
  static RealType _dz_b5 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b6 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b1_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b7 ( const DomVecType &c ) { return _b1_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }
  static RealType _dz_b8 ( const DomVecType &c ) { return _b2_1d ( c[0] ) *_b2_1d ( c[1] ) *_d_b2_1d ( c[2] ); }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const DomVecType &c );
  BASIS_FUNC_TYPE _deriv_basis[3][8];
  BASIS_FUNC_TYPE _basis[8];

  const RealType _hx, _hy, _hz;

public:
  QuocCachedBaseFunctionSet3D( const RealType hx, const RealType hy, const RealType hz ) : _hx ( hx ), _hy ( hy ), _hz (hz) {
    
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;
    _basis[3] = _b4;
    _basis[4] = _b5;
    _basis[5] = _b6;
    _basis[6] = _b7;
    _basis[7] = _b8;

    _deriv_basis[0][0] = _dx_b1;
    _deriv_basis[0][1] = _dx_b2;
    _deriv_basis[0][2] = _dx_b3;
    _deriv_basis[0][3] = _dx_b4;
    _deriv_basis[0][4] = _dx_b5;
    _deriv_basis[0][5] = _dx_b6;
    _deriv_basis[0][6] = _dx_b7;
    _deriv_basis[0][7] = _dx_b8;

    _deriv_basis[1][0] = _dy_b1;
    _deriv_basis[1][1] = _dy_b2;
    _deriv_basis[1][2] = _dy_b3;
    _deriv_basis[1][3] = _dy_b4;
    _deriv_basis[1][4] = _dy_b5;
    _deriv_basis[1][5] = _dy_b6;
    _deriv_basis[1][6] = _dy_b7;
    _deriv_basis[1][7] = _dy_b8;

    _deriv_basis[2][0] = _dz_b1;
    _deriv_basis[2][1] = _dz_b2;
    _deriv_basis[2][2] = _dz_b3;
    _deriv_basis[2][3] = _dz_b4;
    _deriv_basis[2][4] = _dz_b5;
    _deriv_basis[2][5] = _dz_b6;
    _deriv_basis[2][6] = _dz_b7;
    _deriv_basis[2][7] = _dz_b8;

    this->initializeQuadCache( );
  }

  enum { numBaseFuncs = 8 };
  
  RealType evaluate ( int BaseFuncNum, const DomVecType &c ) const {
    return _basis[BaseFuncNum] ( c );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return QuocCachedBaseFunctionSetInterface<DataTypeContainer, QuadType, 8, QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> >::evaluate ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateGradient ( int BaseFuncNum, const DomVecType &c, DomVecType &Gradient ) const {
    Gradient[0] = _deriv_basis[0][BaseFuncNum] ( c ) / _hx;
    Gradient[1] = _deriv_basis[1][BaseFuncNum] ( c ) / _hy;
    Gradient[2] = _deriv_basis[2][BaseFuncNum] ( c ) / _hz;
  }
  
  inline const DomVecType& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return QuocCachedBaseFunctionSetInterface<DataTypeContainer, QuadType, 8, QuocCachedBaseFunctionSet3D<DataTypeContainer, QuadType, ElementType> >::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

};


















///////////////////////////////////////// global affine basefunction set 


//! The basefunctionset for bilinear elements in 2d
template <typename DataTypeContainer>
class QuocGlobalAffineSymGradBaseFunctionSet2D {
    
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::DomVecType DomVecType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
    
  // b1(x) = (x,0)
  // b2(x) = (0,y)
  // b3(x) = (y,x)
  
  static void _b1  ( const DomVecType &c, DomVecType &dest ) { dest[0] = c[0]; dest[1] = 0.0; }
  static void _b2  ( const DomVecType &c, DomVecType &dest ) { dest[0] = 0.0;  dest[1] = c[1]; }
  static void _b3  ( const DomVecType &c, DomVecType &dest ) { dest[0] = c[1]; dest[1] = c[0]; }
  
  static void _symGrad_b1  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,0) = 1.0; }
  static void _symGrad_b2  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(1,1) = 1.0; }
  static void _symGrad_b3  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,1) = 1.0; dest(1,0) = 1.0; }
  
  static RealType _div_b1  ( const DomVecType &c ) { return 1.; }
  static RealType _div_b2  ( const DomVecType &c ) { return 1.; }
  static RealType _div_b3  ( const DomVecType &c ) { return 0.; }

  typedef void ( *BASIS_FUNC_TYPE ) ( const DomVecType &c, DomVecType &dest );
  BASIS_FUNC_TYPE _basis[3];
  typedef void ( *BASIS_FUNC_TYPE_SYMGRAD ) ( const DomVecType &c, DerivativeVectorValuedType &dest );
  BASIS_FUNC_TYPE_SYMGRAD _symGrad_basis[3];
  typedef RealType ( *BASIS_FUNC_TYPE_DIV ) ( const DomVecType &c );
  BASIS_FUNC_TYPE_DIV _div_basis[3];

public:
  QuocGlobalAffineSymGradBaseFunctionSet2D ( ) {
    _basis[0] = _b1; _basis[1] = _b2; _basis[2] = _b3;
    _symGrad_basis[0] = _symGrad_b1; _symGrad_basis[1] = _symGrad_b2; _symGrad_basis[2] = _symGrad_b3;
    _div_basis[0] = _div_b1; _div_basis[1] = _div_b2; _div_basis[2] = _div_b3;
  }

  enum { numBaseFuncs = 3 };
  
  void evaluate ( int BaseFuncNum, const DomVecType &c, DomVecType &dest ) const { _basis[BaseFuncNum] ( c, dest );}
  void evaluateSymGrad ( int BaseFuncNum, const DomVecType &c, DerivativeVectorValuedType &dest ) const { _symGrad_basis[BaseFuncNum] ( c, dest );}
  RealType evaluateDiv ( int BaseFuncNum, const DomVecType &c ) const { return _div_basis[BaseFuncNum] (c);}
  
};


template <typename DataTypeContainer>
class QuocGlobalAffineSymGradBaseFunctionSet3D{

  typedef typename DataTypeContainer::RealType RealType;  
  typedef typename DataTypeContainer::DomVecType DomVecType;
  typedef typename DataTypeContainer::DerivativeVectorValuedType DerivativeVectorValuedType;
  
  static void _b1   ( const DomVecType &c, DomVecType & dest ) { dest[0] = c[0]; dest[1] = 0.0; dest[2] = 0.0; }
  static void _b2   ( const DomVecType &c, DomVecType & dest ) { dest[0] = 0.0; dest[1] = c[1]; dest[2] = 0.0; }
  static void _b3   ( const DomVecType &c, DomVecType & dest ) { dest[0] = 0.0; dest[1] = 0.0; dest[2] = c[2]; }
  static void _b4   ( const DomVecType &c, DomVecType & dest ) { dest[0] = c[1]; dest[1] = c[0]; dest[2] = 0.0; }
  static void _b5   ( const DomVecType &c, DomVecType & dest ) { dest[0] = c[2]; dest[1] = 0.0; dest[2] = c[0]; }
  static void _b6   ( const DomVecType &c, DomVecType & dest ) { dest[0] = 0.0; dest[1] = c[2]; dest[2] = c[1]; }
  
  static void _symGrad_b1  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,0) = 1.0; }
  static void _symGrad_b2  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(1,1) = 1.0; }
  static void _symGrad_b3  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(2,2) = 1.0; }
  static void _symGrad_b4  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,1) = 1.0; dest(1,0) = 1.0; }
  static void _symGrad_b5  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(0,2) = 1.0; dest(2,0) = 1.0; }
  static void _symGrad_b6  ( const DomVecType &c, DerivativeVectorValuedType &dest ) { dest.setZero(); dest(1,2) = 1.0; dest(2,1) = 1.0; }
  
  static RealType _div_b1  ( const DomVecType &c ) { return 1.; }
  static RealType _div_b2  ( const DomVecType &c ) { return 1.; }
  static RealType _div_b3  ( const DomVecType &c ) { return 1.; }
  static RealType _div_b4  ( const DomVecType &c ) { return 0.; }
  static RealType _div_b5  ( const DomVecType &c ) { return 0.; }
  static RealType _div_b6  ( const DomVecType &c ) { return 0.; }

  typedef void ( *BASIS_FUNC_TYPE ) ( const DomVecType &c, DomVecType & dest );
  BASIS_FUNC_TYPE _basis[6];
  typedef void ( *BASIS_FUNC_TYPE_SYMGRAD ) ( const DomVecType &c, DerivativeVectorValuedType &dest );
  BASIS_FUNC_TYPE_SYMGRAD _symGrad_basis[6];
  typedef RealType ( *BASIS_FUNC_TYPE_DIV ) ( const DomVecType &c );
  BASIS_FUNC_TYPE_DIV _div_basis[6];

public:
  QuocGlobalAffineSymGradBaseFunctionSet3D( ) {
    
    _basis[0] = _b1; _basis[1] = _b2; _basis[2] = _b3; _basis[3] = _b4; _basis[4] = _b5; _basis[5] = _b6;
    _symGrad_basis[0] = _symGrad_b1; _symGrad_basis[1] = _symGrad_b2; _symGrad_basis[2] = _symGrad_b3; _symGrad_basis[3] = _symGrad_b4; _symGrad_basis[4] = _symGrad_b5; _symGrad_basis[5] = _symGrad_b6;
    _div_basis[0] = _div_b1; _div_basis[1] = _div_b2; _div_basis[2] = _div_b3; _div_basis[3] = _div_b4; _div_basis[4] = _div_b5; _div_basis[5] = _div_b6;
  }

  enum { numBaseFuncs = 6 };
  
  void evaluate ( int BaseFuncNum, const DomVecType &c, DomVecType &dest ) const {_basis[BaseFuncNum] ( c, dest );}
  void evaluateSymGrad ( int BaseFuncNum, const DomVecType &c, DerivativeVectorValuedType &dest ) const { _symGrad_basis[BaseFuncNum] ( c, dest );}
  RealType evaluateDiv ( int BaseFuncNum, const DomVecType &c ) const { return _div_basis[BaseFuncNum] (c);}


};


}//end namespace quocFE


#endif //__QUOCMESHBASEFUNCTIONSETS_H
