#ifndef __DERIVATIVETESTER_H
#define __DERIVATIVETESTER_H

#include <general.h>
#include <loadAndSave.h>
#include <energyDefines.h>

#include <tikzPlotter.h>

// Compare DE(testPoint)(testDirection) with diffQuotient 1/h ( E(testPoint + h testDirection ) - E(testPoint) )
template<typename DataTypeContainer >
class DerivativeTester {

protected:

  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::SparseMatrixType MatrixType;
    
  const aol::NonlinearEnergyOp<DataTypeContainer> &_energyOp;
  const RealType _stepSize;
  const VectorType _testPoint;

  const int _outputLevel;
  const RealType _tolerance;
  
public:
  DerivativeTester ( const aol::NonlinearEnergyOp<DataTypeContainer> &E, const RealType stepSize, const VectorType & testPoint, const int outputLevel = 1, const RealType tolerance = 1.e-5  ) 
  : _energyOp ( E ), _stepSize ( stepSize ), _testPoint ( testPoint ), _outputLevel ( outputLevel ), _tolerance ( tolerance ) {}        
   
  RealType testFirstDerivative_SingleDirection ( const VectorType & testDirection ) const {
    
    // evalualte J(m)
    RealType energy;
    _energyOp.evaluateEnergy( _testPoint, energy );
    
    //evaluate J(m + h testDirection )
    RealType energyShifted;
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection;    
    _energyOp.evaluateEnergy( shiftedPoint, energyShifted );
    RealType diffQuotient = ( energyShifted - energy ) / (_stepSize ) ;
    if(_outputLevel > 1 ) {
     cout << "E(p+ht) = " << energyShifted << endl;
     cout << "E(p)    = " << energy << endl;
     cout << "approxDerivative = " << diffQuotient << endl;
    }
    return diffQuotient;
  }
  
  void testFirstDerivative_AllDirections( ) const {
      
      unsigned numDirections = _testPoint.size();
      
      //evaluate DJ(m)
      VectorType derivative ( numDirections );
      _energyOp.evaluateJacobian( _testPoint, derivative );
      
      // diffQuotient
      VectorType approxDerivative ( numDirections );
      
      VectorType testDirection( numDirections );
      cout << "num directions to test = " << numDirections << endl;
      aol::printVector<VectorType> ( _testPoint, 10, "test point" );
      for( unsigned i=0; i<numDirections; ++i ){
        if( _outputLevel > 2 )  cout << "direction = " << i << endl;
        testDirection.setZero(); testDirection[i] = 1.;
        approxDerivative[i] = testFirstDerivative_SingleDirection ( testDirection );
      }
      
      VectorType diff = derivative - approxDerivative;
      if( _outputLevel > 0 ){
          aol::printVector<VectorType> ( derivative, 10, "gradient" );
          aol::printVector<VectorType> ( approxDerivative, 10, "diffQuotient" );
          aol::printVector<VectorType> ( diff, 10, "grad - diffQuotient" );
      }
      for( int i=0; i<numDirections; ++i ){
       if( std::abs( diff[i] ) > _tolerance ){
           cout << std::fixed << std::setprecision( 12 ) << "error in direction " << i << ": " << std::abs( diff[i] ) << endl 
                                                         << " deriv  = " << derivative[i] << endl 
                                                         << " approx = " << approxDerivative[i] << endl;
           if( std::abs( derivative[i] ) > 1.e-15 ){
              cout << std::fixed << std::setprecision( 12 ) << " relativ error in direction " << i << ": " << std::abs( diff[i] / derivative[i] ) << endl;
           }
       }
      }
  }
  
    
  RealType testSecondDerivative_SingleDirection ( const VectorType & testDirection1, const VectorType & testDirection2 ) const {
      
    // evalualte DJ(m)
    VectorType derivative ( testDirection1.size() );
    _energyOp.evaluateJacobian( _testPoint, derivative );
    
    //evaluate DJ(m + h testDirection )
    VectorType derivativeShifted ( testDirection1.size() );
    VectorType shiftedPoint = _testPoint + _stepSize * testDirection2;    
    _energyOp.evaluateJacobian( shiftedPoint, derivativeShifted );
    RealType diffQuotient = ( derivativeShifted.dot(testDirection1) - derivative.dot(testDirection1) ) / (_stepSize ) ;
    
    return diffQuotient;
  }
  
  void testSecondDerivative_AllDirections( ) const {
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      MatrixType hessian ( numDirections, numDirections ); 
      Eigen::MatrixXd approxHessian (numDirections, numDirections );
      _energyOp.evaluateHessian( _testPoint, hessian );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSecondDerivative_SingleDirection ( testDirection1, testDirection2 );
          }
      }
      
      Eigen::MatrixXd hessianDense =  Eigen::MatrixXd( hessian );
      Eigen::MatrixXd diff = hessianDense - approxHessian;
      if( _outputLevel > 0 ){
          cout << "hessian = " << endl << Eigen::MatrixXd( hessian ) << endl << endl;
          cout << "diffQuotient = " << endl << approxHessian << endl;
          cout << endl << endl << endl << "hessian - diffQuotient = " << endl << diff << endl;
      }
      for( unsigned i=0; i<numDirections; ++i ){
          for( unsigned j=0; j<numDirections; ++j ){
            if( std::abs( diff(i,j) ) > 1.e-5 ) 
                cout << aol::color::red << "error in direction (" << i << " , " << j << " ) :\t" << "hessian = " << hessianDense(i,j) << "\t approx = " << approxHessian(i,j) << aol::color::reset << endl; 
            else{ 
              if( std::abs( hessianDense(i,j) ) > 1.e-16  ) 
                cout << aol::color::green << "correct in direction (" << i << " , " << j << " ) :\t" << "hessian = " << hessianDense(i,j) << "\t approx = " << approxHessian(i,j) << aol::color::reset << endl;
            }
          }
      }
   }
   
   
   void testSecondDerivative_AllDirections( const string &saveDirectory, const string &fileName ) const {
      
      unsigned numDirections = _testPoint.size();
      //evaluate D2J(m)
      MatrixType hessian ( numDirections, numDirections ); 
      Eigen::MatrixXd approxHessian (numDirections, numDirections );
      _energyOp.evaluateHessian( _testPoint, hessian );
      VectorType testDirection1( numDirections ), testDirection2( numDirections );
      for( unsigned i=0; i<numDirections; ++i ){
          testDirection1.setZero();
          testDirection1[i] = 1.;
          for( unsigned j=0; j<numDirections; ++j ){
              testDirection2.setZero();
              testDirection2[j] = 1.;
              approxHessian(i,j) = testSecondDerivative_SingleDirection ( testDirection1, testDirection2 );
          }
      }
      
      Eigen::MatrixXd hessianDense =  Eigen::MatrixXd( hessian );
      Eigen::MatrixXd diff = hessianDense - approxHessian;
      
      //plot with latex
      std::ofstream out ( aol::strprintf( "%s/%s", saveDirectory.c_str(), fileName.c_str()  ) );
      
      TikzPlotterHelperClass<RealType> tikzHelper;
        
      tikzHelper.generateIncludes( out, 400, 100, 0.5, 0.5, 0.5, 0.5 );
      tikzHelper.generateBeginDocument( out );

      out << "\\setcounter{MaxMatrixCols}{" << numDirections << "}" << endl;
      out << "\\begin{align}" << endl
          << "\\begin{pmatrix}" << endl;
          
     for( unsigned i=0; i<numDirections; ++i ){
          for( unsigned j=0; j<numDirections; ++j ){
            if( std::abs( diff(i,j) ) > 1.e-5 ) 
                out << "{\\color{red} " << hessianDense(i,j) << " } ";
            else{
              if( std::abs( hessianDense(i,j) ) > 1.e-16  ) 
                out << "{\\color{green} " << hessianDense(i,j) << " } ";
              else
                out << "{\\color{black} " << hessianDense(i,j) << " } ";
            }
            if (j < numDirections - 1) out << " & ";
          }
          out << " \\\\" << endl;
      }
          
      out << "\\end{pmatrix}" << endl;
      out << "\\end{align}" << endl;
        
      tikzHelper.generateEndDocument( out );
      
      std::ofstream BashFilePDF ( aol::strprintf ( "%s/%s.sh", saveDirectory.c_str (), "generatePDF"  ) );
      BashFilePDF << "cd " << saveDirectory.c_str() << endl;
      BashFilePDF << "lualatex " << fileName.c_str() << endl;
      string systemCommand = "bash " + saveDirectory + "/generatePDF.sh";
      cout << "systemCommand = " << systemCommand << endl;
      bool failed;
      failed= ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
      if ( failed ) cerr << "programm returned an error." << endl;
      
   }
                      
};

#endif //__DERIVATIVETESTER_H
