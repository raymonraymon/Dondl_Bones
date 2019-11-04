#ifndef __MULTIPLELOAD__H
#define __MULTIPLELOAD__H


//======================================================================================================================================
//================================= weightFunctionMultipleLoad  ===============================================================================
//======================================================================================================================================
    
// w (sum_i beta_i/( x_i^p) )^(1/q)
template <typename RealType>
class WeightFunctionMultipleLoad {
protected:
    const RealType _w, _p, _q;
    bool _qIsOne;
    const std::vector<RealType> &_weights;
public:
    
    
    WeightFunctionMultipleLoad ( ) : _w( 1. ), _p ( 1. ), _q( 1. ), _weights{} {
        _qIsOne = true;
    }
    
    WeightFunctionMultipleLoad ( const RealType w, const RealType p, const RealType q, const std::vector<RealType> &weights ) : _w( w ), _p( p ), _q( q ), _weights( weights ) {
        _qIsOne = false;
        if( std::abs( _q - 1.0 ) < 2.e-16 ) _qIsOne = true;
    }
    
    RealType evaluate ( const std::vector<RealType> &x ) const {
        RealType aux = 0.0;
        for( int i=0; i<x.size(); ++i ) aux += _weights[i] / std::pow(x[i],_p);
        if( _qIsOne ) return _w * aux;
        return _w * std::pow( aux, 1. / _q );
    }
    RealType evaluateDerivative ( const std::vector<RealType> &x, const int direction ) const {
        RealType outerDeriv = 0.0;
        if( _qIsOne ){
            outerDeriv = 1.0;
        }else{
            RealType tmp = 0.0;
            for( int i=0; i<x.size(); ++i ) tmp += _weights[i] / std::pow(x[i],_p);
            outerDeriv = std::pow( tmp, 1. / _q - 1. ) / _q;
        }
        return -1.0 * _w * outerDeriv * _p * _weights[direction] / std::pow(x[direction],_p+1);   
    }
    
    string description ( ) const {
        std::ostringstream beta;
        beta << _w
             << " \\left( \\sum_i \\frac{\\beta_i}{ x_i^{"
             << _p
             << " }} \\right)^{\\frac{1}{"
             << _q 
             << " }}, \\beta = (";
        for( int i=0; i < _weights.size(); ++i ){
             beta << _weights[i];
             if( i < _weights.size() - 1 ) beta << ", ";
        }
        beta << ")";
        std::string betaString = beta.str();
//         return aol::strprintf("%.2f \\left( \\sum_i \\frac{\\beta_i}{ x_i^{%.2f}} \\right)^{\\frac{1}{%.2f}}, %s", _w, _p, _q, betaString.c_str()  ).c_str();
        return betaString;
    }
};


#endif
