#ifndef __LINEARSYSTEMSOLVER_H
#define __LINEARSYSTEMSOLVER_H

#include <general.h>
#include <loadAndSave.h>
#include <unsupported/Eigen/IterativeSolvers>
#include <EigenBiCGSTAB_New.h>
#include <EigenGMRES_New.h>



template<typename DataTypeContainer, typename IterativeLinearSolverType = typename DataTypeContainer::IterativeLinearSolverType>
class IterativeLinearSystemSolver{
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  MatrixType;
    
    IterativeLinearSolverType _linearSolver;
    
    const string _errorMessage;
    const string _saveDirectory;
    const bool _debug;
    
public:
     
  IterativeLinearSystemSolver( const string &errorMessage, const string &saveDirectory = "", const int outputLevel = 0, const bool debug = true  ) : 
  _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ), _debug ( debug ) { }   
  
  //
  bool prepareSolver( const MatrixType& systemMatrix ) {
      _linearSolver.compute( systemMatrix );
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solveForPreparedMatrix( VectorType& solution, const VectorType& rhs  ) {
    solution = _linearSolver.solve( rhs );
    if(_linearSolver.info()!=Eigen::Success ) {
        cout << "error comes from solve in " << _errorMessage << endl;
        cout << "#iterations:     " << _linearSolver.iterations() << endl;
        cout << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error()      << endl;
        return false;  
    }
            
    return true;
  }
  
  
  bool solveForPreparedMatrixWithGuess( VectorType& solution, const VectorType& rhs, const VectorType& guess ) {
      
    solution = _linearSolver.solveWithGuess( rhs, guess );
    if(_linearSolver.info()!=Eigen::Success ) {
        cout << "error comes from solve in " << _errorMessage << endl;
        cout << "#iterations:     " << _linearSolver.iterations() << endl;
        cout << "estimated error: " << std::setprecision( 32 ) << _linearSolver.error()      << endl;
        return false;  
    }
    return true;
  }
  
  
  // solve Ax = b for x
  void solve( const MatrixType& systemMatrix, VectorType& solution, const VectorType& rhs, const RealType tolerance = Eigen::NumTraits<RealType>::epsilon(), const RealType maxItersFac = 1.   ) {
      bool succesPrepareSolver = prepareSolver( systemMatrix );
      const RealType newTolerance = sqrt( tolerance * tolerance * systemMatrix.cols() );
      _linearSolver.setTolerance( newTolerance );
      const int maxIters = static_cast<int> ( maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
      _linearSolver.setMaxIterations(maxIters);
      bool success = solveForPreparedMatrix( solution, rhs );
  }
  
  void solveWithGuess( const MatrixType& systemMatrix, VectorType& solution, const VectorType& rhs, const VectorType &guess, const RealType tolerance = Eigen::NumTraits<RealType>::epsilon(), const RealType maxItersFac = 1.   ) {
      bool succesPrepareSolver = prepareSolver( systemMatrix );
      const RealType newTolerance = sqrt( tolerance * tolerance * systemMatrix.cols() );
      _linearSolver.setTolerance( newTolerance );
      const int maxIters = static_cast<int> ( maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
      _linearSolver.setMaxIterations(maxIters);
      bool success = solveForPreparedMatrixWithGuess( solution, rhs, guess );
  }
   
};








template<typename DataTypeContainer, typename DirectLinearSolverType = typename DataTypeContainer::DirectLinearSolverType>
class DirectLinearSystemSolver{
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  MatrixType;
    
    DirectLinearSolverType _linearSolver;
    
    const string _errorMessage;
    const string _saveDirectory;
    
public:
     
  DirectLinearSystemSolver( const string &errorMessage, const string &saveDirectory = "" ) : 
  _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) { }   
  
  
  void analyzePattern( const MatrixType& systemMatrix ) {
      _linearSolver.analyzePattern(systemMatrix);
  }
  
  //
  bool prepareSolver( const MatrixType& systemMatrix ) {
      _linearSolver.compute(systemMatrix);
      if(_linearSolver.info()!=Eigen::Success) {
          cout << "prepare linear solver failed" << endl;
          return false;
      }
      return true;
  }
  
  //
  bool solve( VectorType& solution, const VectorType& rhs ) {
    bool succes = true;
    solution = _linearSolver.solve( rhs );
    return succes;
  }
  
  // solve Ax = b for x
  void solve( const MatrixType& systemMatrix, VectorType& solution, const VectorType& rhs ) {
      prepareSolver( systemMatrix );
      bool success = solve( solution, rhs );
  }
   
};  









#endif //__LINEARSYSTEMSOLVER_H
