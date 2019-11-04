#ifndef __QUOCQUADRATURE_H
#define __QUOCQUADRATURE_H


namespace quocFE {

template <typename RealType>
class SimpsonQuadrature1D {
public:
  enum { numQuadPoints = 3 };

  SimpsonQuadrature1D(){}
  
  inline const RealType getRefCoord ( int quadpoint ) const  {
    static const RealType c[3] = { 0.0, 0.5, 1.0 };
    return c[quadpoint];
  }

  inline RealType getWeight ( int quadpoint ) const { 
      static const RealType w[3] = { 1./6., 2./3., 1./6. };
      return w[quadpoint];
   }
};

template <typename RealType, typename DomVecType>
class SimpsonQuadrature1DVec {
public:
  enum { numQuadPoints = 3 };

  SimpsonQuadrature1DVec(){
        _points[0][0] = 0.0;
        _points[1][0] = 0.5;
        _points[2][0] = 1.0;
        _weights[0] = 1./6.; _weights[1] = 2./3.; _weights[2] = 1./6.;
  }
  
  inline const DomVecType &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}
   
   protected:
  DomVecType _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
};


template <typename RealType, typename DomVecType>
class SimpsonQuadrature2D {
public:

  enum { numQuadPoints = SimpsonQuadrature1D<RealType>::numQuadPoints * SimpsonQuadrature1D<RealType>::numQuadPoints };

  SimpsonQuadrature2D( ) {

    SimpsonQuadrature1D<RealType> quad1D;

    for ( int k = 0, i = 0; i < SimpsonQuadrature1D<RealType>::numQuadPoints; i++ )
      for ( int j = 0; j < SimpsonQuadrature1D<RealType>::numQuadPoints; j++, k++ ) {
        _points[k][0] = quad1D.getRefCoord( i );
        _points[k][1] = quad1D.getRefCoord ( j );
        _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j );
      }
  }

  inline const DomVecType &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}

protected:
  DomVecType _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
};





template <typename RealType, typename DomVecType>
class SimpsonQuadrature3D {
public:

  enum { numQuadPoints = SimpsonQuadrature1D<RealType>::numQuadPoints * SimpsonQuadrature1D<RealType>::numQuadPoints * SimpsonQuadrature1D<RealType>::numQuadPoints };

  SimpsonQuadrature3D( ) {

    SimpsonQuadrature1D<RealType> quad1D;

    for ( int k = 0, i = 0; i < SimpsonQuadrature1D<RealType>::numQuadPoints; i++ ) {
      for ( int j = 0; j < SimpsonQuadrature1D<RealType>::numQuadPoints; j++ ) {
        for ( int l = 0; l < SimpsonQuadrature1D<RealType>::numQuadPoints; l++, k++ ) {
          _points[k][0] = quad1D.getRefCoord ( i );
          _points[k][1] = quad1D.getRefCoord ( j );
          _points[k][2] = quad1D.getRefCoord ( l );
          _weights[k] = quad1D.getWeight ( i ) * quad1D.getWeight ( j ) * quad1D.getWeight ( l );
        }
      }
    }
  }

  inline const DomVecType &getRefCoord ( int QuadPoint ) const { return _points[ QuadPoint ];}
  inline RealType getWeight ( int QuadPoint ) const { return _weights[ QuadPoint ];}

protected:
  DomVecType _points [ numQuadPoints ];
  RealType   _weights[ numQuadPoints ];
};

}



#endif
