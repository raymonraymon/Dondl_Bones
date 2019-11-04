#ifndef __GENERALLIB_H
#define __GENERALLIB_H

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC system_header
#endif

// C standard libraries
#include <cmath>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <initializer_list>

// STL
#include <complex>
#include <limits>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include <typeinfo>

#include <chrono>
#include <unistd.h> //for sleep function

using namespace std;

#include <stdexcept>
#include <stdint.h>

#include <functional>

//Eigen
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <Eigen/UmfPackSupport>
#include <Eigen/CholmodSupport>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>




namespace aol {

//! table of ansi color codes
namespace color {
const string reset       = "\033[0;0m";
const string invert      = "\033[0;7m";
const string black       = "\033[0;30m";
const string red         = "\033[0;31m";
const string green       = "\033[0;32m";
const string brown       = "\033[0;33m";
const string blue        = "\033[0;34m";
const string purple      = "\033[0;35m";
const string cyan        = "\033[0;36m";
const string light_grey  = "\033[0;37m";
const string dark_grey   = "\033[1;30m";
const string light_red   = "\033[1;31m";
const string light_green = "\033[1;32m";
const string yellow      = "\033[1;33m";
const string light_blue  = "\033[1;34m";
const string pink        = "\033[1;35m";
const string light_cyan  = "\033[1;36m";
const string white       = "\033[1;37m";
const string beep        = "\007";
const string error       = beep + red;
const string ok          = green;
}





//! Give back formatted string, analogously to sprintf, but save the long way 'round with char arrays.
string strprintf(const char * format, ...) {
  // declare variable argument list
  va_list az;
  // copy from my input into this list(second argument is nothing really used, but the last named argument of myself)
  va_start(az, format);
  // give this argument list variable to vscprintf instead of my own arguments (that is in what functions like vsprintf differ from functions like sprintf)
  const int sizeNeeded = vsnprintf ( NULL, 0, format, az ) + 1;
  // restore stack into clean state:
  va_end(az);

  char *buffer = new char[sizeNeeded];

  va_start(az, format);
  vsprintf (buffer, format, az);
  va_end(az);

  // automatic return type conversion into string:
  string ret = buffer;
  delete[] buffer;
  return ret;
}

// Returns minimum/maximium
template<class T> inline T Min ( const T a, const T b ) { return ( ( a < b ) ? a : b ); }
template<class T> inline T Max ( const T a, const T b ) { return ( ( a < b ) ? b : a ); }
template<class T> inline T Clamp ( const T Value, const T Min, const T Max ) { return ( aol::Max ( aol::Min ( Value, Max ), Min ) ); } // Returns Value clamped into [Min,Max].
template<class T> inline T Sqr (const T a) { return a * a; }
template<class T> inline T Cub (const T a) { return a * a * a; }

template<typename RealType, typename VectorType>
void thresholdVector ( const VectorType & vec, VectorType & thresholdVec, RealType lowerBound, RealType upperBound, RealType threshold ){
    for( int i=0; i<vec.size(); ++i ){
      if( vec[i] < threshold ) thresholdVec[i] = lowerBound;
      else thresholdVec[i] = upperBound;
    }
}

template<typename RealType, typename Matrix>
RealType ddProd( const Matrix &MatA, const Matrix &MatB ){
    RealType out = 0.0;
    for( int i=0; i<MatA.rows(); ++i ) 
        for( int j=0; j<MatB.cols(); ++j )
            out += MatA(i,j) * MatB(i,j);
        
    return out;
}

} // namespace aol

#endif

