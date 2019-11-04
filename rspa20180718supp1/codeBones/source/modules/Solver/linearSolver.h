#ifndef __LINEARSOLVER_H
#define __LINEARSOLVER_H

#include <general.h>
#include <loadAndSave.h>
#include <unsupported/Eigen/IterativeSolvers>
#include <EigenBiCGSTAB_New.h>
#include <EigenGMRES_New.h>


enum LinearSolverMethod {
    //direct solver
    UmfPackLU = 1,
    CholmodSupernodalLLT = 2,
    //iterative solver
    ConjugateGradient = 11,
    BiCGSTAB = 21,
    BiCGSTAB_WithGuess = 22,
    BiCGSTAB_ILU = 23,
    BiCGSTAB_NEW = 24,
    EigenGMRES = 31,
    EigenGMRES_WithGuess = 32,
    EigenGMRES_ILU = 33,
    EigenGMRES_NEW = 34
};





template<typename DataTypeContainer>
class LinearSystemSolver{
    
protected:
    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType  VectorType;
    typedef typename DataTypeContainer::SparseMatrixType  MatrixType;

    int _solverType;

    Eigen::UmfPackLU<MatrixType> _LUSolver;
    Eigen::CholmodDecomposition<MatrixType> _CholmodSolver;
    Eigen::ConjugateGradient<MatrixType, Eigen::Lower|Eigen::Upper> _CGSolver;
    Eigen::BiCGSTAB<MatrixType> _BiCGSolver;
    Eigen::BiCGSTAB_NEW<MatrixType> _BiCGSolverNew;
    Eigen::GMRES<MatrixType> _GMRESSolver;
    Eigen::GMRES_NEW<MatrixType> _GMRESSolverNew;
    
//     const RealType _tolerance;
//     const Eigen::Index _maxIterations;
    
    const string _errorMessage;
    const string _saveDirectory;

public:
  LinearSystemSolver( LinearSolverMethod linearSolver = UmfPackLU ) : _solverType(linearSolver) {}
     
  LinearSystemSolver( LinearSolverMethod linearSolver, const string &errorMessage, const string &saveDirectory  ) : 
  _solverType(linearSolver), _errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) {}   
  
//   LinearSystemSolver( LinearSolverMethod linearSolver, const RealType tolerance, const Eigen::Index maxIterations, const Eigen::Index maxIterations, const string &errorMessage, const string &saveDirectory  ) : 
//   _solverType(linearSolver), _tolerance(tolerance), _maxIterations(maxIterations),_errorMessage( errorMessage ), _saveDirectory ( saveDirectory ) {}   
  
  //
  void prepareSolver( const MatrixType& systemMatrix ) {
    switch( _solverType ){
       case UmfPackLU:{
           _LUSolver.compute(systemMatrix);
           if( _LUSolver.info()!=Eigen::Success ){
//             _LUSolver.umfpackReportInfo();
               throw std::invalid_argument ( aol::strprintf ( "prepare UmfPackLU solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
               return;
           }
           break;
       }
      
      case CholmodSupernodalLLT:{
          _CholmodSolver.compute(systemMatrix);
          if(_CholmodSolver.info()!=Eigen::Success){
              throw std::invalid_argument ( aol::strprintf ( "prepare Cholmod solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
          }
          break;   
      }
      
      case ConjugateGradient: {
            _CGSolver.compute(systemMatrix);
            if(_CGSolver.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "prepare Cg solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
            }
            break;
      }
      
      case BiCGSTAB : {
            _BiCGSolver.compute(systemMatrix);
            if(_BiCGSolver.info()!=Eigen::Success) {
              cout << "error comes from prepare matrix in " << _errorMessage << endl;
              throw std::invalid_argument ( aol::strprintf ( "prepare BICGSTAB solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
            }
            break;
      }
      
     case BiCGSTAB_NEW : {
            _BiCGSolverNew.compute(systemMatrix);
            if(_BiCGSolverNew.info()!=Eigen::Success) {
              cout << "error comes from prepare matrix in " << _errorMessage << endl;
              throw std::invalid_argument ( aol::strprintf ( "prepare BICGSTAB solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
            }
            break;
      }
      
      case EigenGMRES : {
          _GMRESSolver.compute( systemMatrix );
          if(_GMRESSolver.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "prepare GMRES solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
            }
          break;
      }
      
      case EigenGMRES_NEW : {
          _GMRESSolverNew.compute( systemMatrix );
          if(_GMRESSolverNew.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "prepare GMRES solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
            }
          break;
      }
      
      default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel %s. In File %s at line %d.", _solverType,  __FILE__, __LINE__ ).c_str() );
          break;
    }
  }
  
  //
  bool solve( VectorType& solution, const VectorType& rhs, const RealType tolerance = Eigen::NumTraits<RealType>::epsilon()  ) {
    bool succes = true;
    switch( _solverType ){

       case UmfPackLU:{
           solution = _LUSolver.solve( rhs );
           break;
       }
      
      case CholmodSupernodalLLT:{
           solution = _CholmodSolver.solve( rhs );         
           break;   
      }
      
      case ConjugateGradient: {
            solution = _CGSolver.solve( rhs );
            break;
      }
      
      case BiCGSTAB : {
            solution = _BiCGSolver.solve( rhs );
            if(_BiCGSolver.info()!=Eigen::Success) {
                cout << "error comes from solve in " << _errorMessage << endl;
                std::cout << "#iterations:     " << _BiCGSolver.iterations() << std::endl;
                std::cout << "estimated error: " << _BiCGSolver.error()      << std::endl;
                return false;  
            }
            break;
      }
      
     case BiCGSTAB_NEW : {
         _BiCGSolverNew.setTolerance( tolerance );
         solution = _BiCGSolverNew.solve( rhs );
            if(_BiCGSolverNew.info()!=Eigen::Success ) {
                cout << "error comes from solve in " << _errorMessage << endl;
                std::cout << "#iterations:     " << _BiCGSolverNew.iterations() << std::endl;
                std::cout << "estimated error: " << _BiCGSolverNew.error()      << std::endl;
                return false;  
            }
            break;
      }
      
      case EigenGMRES : {
            solution = _GMRESSolver.solve( rhs );
            if(_GMRESSolver.info()!=Eigen::Success) {
                 cout << "error comes from " << _errorMessage << endl;
                 std::cout << "#iterations:     " << _GMRESSolver.iterations() << std::endl;
                 std::cout << "estimated error: " << _GMRESSolver.error()      << std::endl;
                return false;
            } 
            break;
      }
      
      case EigenGMRES_NEW : {
            _GMRESSolverNew.setTolerance( tolerance );
            solution = _GMRESSolverNew.solve( rhs );
            if( (_GMRESSolverNew.info()!=Eigen::Success) && (_GMRESSolverNew.error() > 1.e+3 * _GMRESSolverNew.tolerance() ) ) {
                 cout << "error comes from " << _errorMessage << endl;
                 std::cout << "#iterations:     " << _GMRESSolverNew.iterations() << std::endl;
                 std::cout << "estimated error: "  << std::setprecision(25) << _GMRESSolverNew.error()      << std::endl;
                return false;
            } 
            break;
      }
      
      default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel %s. In File %s at line %d.", _solverType,  __FILE__, __LINE__ ).c_str() );
          break;
    }
    return succes;
  }
  
  
  
  // solve Ax = b for x
  void solve( const MatrixType& systemMatrix, VectorType& solution, const VectorType& rhs, const RealType tolerance = Eigen::NumTraits<RealType>::epsilon()  ) {
      prepareSolver( systemMatrix );
      bool success = solve( solution, rhs, tolerance );
  }
  
  void solveWithGuess( const MatrixType& systemMatrix, VectorType& solution, const VectorType& rhs, const VectorType &guess  ) {
      prepareSolver( systemMatrix );
      solveWithGuess( solution, rhs, guess );
  }
   
};  





















template<typename DataTypeContainer>
void solveLinearSystem ( const typename DataTypeContainer::SparseMatrixType & systemMatrix , 
                         typename DataTypeContainer::VectorType &solution, 
                         const typename DataTypeContainer::VectorType &rhs,
                         const LinearSolverMethod linearSolver = UmfPackLU ) {

    typedef typename DataTypeContainer::RealType RealType;
    typedef typename DataTypeContainer::VectorType VectorType;
    typedef typename DataTypeContainer::SparseMatrixType SparseMatrixType;
                             
    switch( linearSolver ){

      case UmfPackLU:{
          Eigen::UmfPackLU<SparseMatrixType> solver;
          solver.compute(systemMatrix);
          if(solver.info()!=Eigen::Success) {
              if(solver.info() == Eigen::NumericalIssue ) cout << "The provided data did not satisfy the prerequisites." << endl;
              if(solver.info() == Eigen::NoConvergence ) cout << "Iterative procedure did not converge." << endl;
              if(solver.info() == Eigen::InvalidInput ) cout << "The inputs are invalid, or the algorithm has been improperly called. When assertions are enabled, such errors trigger an assert." << endl;
              throw std::invalid_argument ( aol::strprintf ( "UmfPackLU solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
          }
          solution = solver.solve( rhs );
          break;
      }
      
      
      case CholmodSupernodalLLT:{
          Eigen::CholmodSupernodalLLT<SparseMatrixType> solver;
          solver.compute(systemMatrix);
          if(solver.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "CholmodSupernodalLLT solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
          }
          solution = solver.solve( rhs );         
          break;   
      }
      
      case ConjugateGradient: {
            Eigen::ConjugateGradient<SparseMatrixType, Eigen::Lower|Eigen::Upper> solver;
            solver.compute(systemMatrix);
            if(solver.info()!=Eigen::Success) {
              throw std::invalid_argument ( aol::strprintf ( "Cg solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
              return;
            }
            solution = solver.solve( rhs );
            //std::cout << "#iterations:     " << cg.iterations() << std::endl;
            //std::cout << "estimated error: " << cg.error()      << std::endl;
            break;
      }
      
      case BiCGSTAB : {
            Eigen::BiCGSTAB<SparseMatrixType> solver;
            solver.compute(systemMatrix);
            if(solver.info()!=Eigen::Success) {
                 throw std::invalid_argument ( aol::strprintf ( "BiCGStab solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                 return;
            }
            solution = solver.solve( rhs );
            break;
      }
      
    case BiCGSTAB_ILU : {
            Eigen::BiCGSTAB<SparseMatrixType,Eigen::IncompleteLUT<RealType,int>> solver;
            solver.compute(systemMatrix);
            if(solver.info()!=Eigen::Success) {
                 throw std::invalid_argument ( aol::strprintf ( "BiCGStab solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
                 return;
            }
            solution = solver.solve( rhs );
            break;
      }
      
      
    case BiCGSTAB_NEW : {
        Eigen::BiCGSTAB_NEW<SparseMatrixType> solver;
        const RealType tolerance = 3.e-16;
        const RealType newTolerance = sqrt( tolerance * tolerance * systemMatrix.cols() );
        solver.setTolerance( newTolerance );
        const RealType maxItersFac = 0.25;
        const int maxIters = static_cast<int> ( maxItersFac * static_cast<RealType> ( systemMatrix.cols() ) );
        solver.setMaxIterations(maxIters);
        
        solver.compute(systemMatrix);
        if(solver.info()!=Eigen::Success) {
            throw std::invalid_argument ( aol::strprintf ( "BiCGStab solver failed in File %s at line %d.", __FILE__, __LINE__ ).c_str() );
            return;
        }
        solution = solver.solve( rhs );
        VectorType residual = systemMatrix * solution;
        residual -= rhs;
        cout << "error = " << residual.norm() << endl;
        break;
    }
      
      case EigenGMRES : {
          Eigen::GMRES<SparseMatrixType > solver(systemMatrix);
          solution = solver.solve( rhs );
          std::cout << "#iterations:     " << solver.iterations() << std::endl;
          std::cout << "estimated error: " << solver.error()      << std::endl;
          break;
      }
      
     case EigenGMRES_ILU : {
          Eigen::GMRES<SparseMatrixType, Eigen::IncompleteLUT<RealType,int>  > solver(systemMatrix);
          solution = solver.solve( rhs );
          std::cout << "#iterations:     " << solver.iterations() << std::endl;
          std::cout << "estimated error: " << solver.error()      << std::endl;
          break;
      }
      
      default:
          throw std::invalid_argument( aol::strprintf ( "Wrong channel %s. In File %s at line %d.", linearSolver,  __FILE__, __LINE__ ).c_str() );
          break;
    }
}






#endif
