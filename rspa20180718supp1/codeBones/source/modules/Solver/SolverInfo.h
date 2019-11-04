#ifndef __SOLVERINFO_H
#define __SOLVERINFO_H

#include <general.h>


template<typename DataTypeContainer>
class SolverInfo{
    
protected:
    typedef typename DataTypeContainer::RealType  RealType;
    
    RealType _error;
    int _numIterations;
    string _solverStatus;
    RealType _runningTime;

public:
  SolverInfo( ) {}
     
  void setError( const RealType error ) { _error = error; }
  RealType getError ( ) const { return _error; }
  void setNumIterations( const int num ) { _numIterations = num; }
  int getNumIterations ( ) const { return _numIterations; }
  void setSolverStatus( const string status ) { _solverStatus = status; }
  string getSolverStatus ( ) const { return _solverStatus; }
  void setRunningTime ( const RealType time ) { _runningTime = time; }
  RealType getRunningTime ( ) const { return _runningTime; }
};  



#endif
