#ifndef __QUOCDEFINES_H
#define __QUOCDEFINES_H

#include <general.h>
#include <loadAndSave.h>

#include <EigenBiCGSTAB_New.h>
#include <EigenGMRES_New.h>
#include <EigenUmfPackSupport_New.h>
#include <EigenProjectedPCG_New.h>


namespace quocFE {

struct QuocDataTypeContainerBase {
public :
  typedef double          RealType;
  
  typedef Eigen::Matrix<RealType, 2, 2> Matrix22;
  typedef Eigen::Matrix<RealType, 3, 2> Matrix32;
  typedef Eigen::Matrix<RealType, 3, 3> Matrix33;
 
  typedef Eigen::VectorXd  VectorType;
  typedef std::vector<bool>  MaskType;
  typedef Eigen::MatrixXd FullMatrixType;
  

  typedef Eigen::SparseMatrix<RealType, Eigen::RowMajor,int> SparseMatrixType;
  typedef Eigen::Triplet<RealType> TripletType;
  typedef aol::BoostParser ParameterParserType;
  
  typedef Eigen::ConjugateGradient_NEW<SparseMatrixType, Eigen::Lower|Eigen::Upper> IterativeLinearSolverType;
  typedef Eigen::UmfPackLU_NEW<SparseMatrixType> DirectLinearSolverType; //for long int
};
  



struct Quoc1DDataTypeContainer : public QuocDataTypeContainerBase {
public :
  typedef Eigen::Matrix<int, 1, 1> IntVecType;
  typedef Eigen::Matrix<RealType, 1, 1 > DomVecType;
  typedef RealType Point1DType;
  typedef RealType PointType;
  typedef Eigen::Vector2i Indices1DType;

};

struct Quoc2DDataTypeContainer : public QuocDataTypeContainerBase {
public :
  typedef Eigen::Vector2i IntVecType;
  typedef Eigen::Vector2d DomVecType;
  typedef RealType DomVecTypeBoundary;
  typedef Eigen::Vector2d Point2DType;
  typedef Eigen::Vector2d PointType;
  typedef Eigen::Vector4i Indices2DType;
  
  typedef Matrix22 DerivativeVectorValuedType;
  

};


struct Quoc3DDataTypeContainer : public QuocDataTypeContainerBase {
public :
  typedef Eigen::Vector3i IntVecType;
  typedef Eigen::Vector3d DomVecType;
  typedef Eigen::Vector2d DomVecTypeBoundary;
  typedef Eigen::Vector3d Point3DType;
  typedef Eigen::Vector3d PointType;
  typedef Eigen::Matrix< int, 8, 1 > Indices3DType;
  
  typedef Matrix33 DerivativeVectorValuedType;
  
  typedef Eigen::Matrix<RealType, 6, 6 > VoigtMatrixType;
  typedef Eigen::Matrix<RealType,6,1> VoigtVecType; 
  
};



enum QuocBoundaryType {
      NOBOUNDARY = 0,
      ALL = 100,
      //2D and 3D
      LEFT = 1,
      RIGHT = 2,
      TOP = 3,
      BOTTOM = 4,
      //only 3D
      FRONT = 5,
      BACK = 6
};



} //end namespace


#endif //__QUOCMESHDEFINES_H
