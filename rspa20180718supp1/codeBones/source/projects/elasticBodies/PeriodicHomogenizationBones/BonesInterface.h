#ifndef __PERIODICHOMOGENIZATIONBONESINTERFACE_H
#define __PERIODICHOMOGENIZATIONBONESINTERFACE_H

#include <general.h>
#include <VTKPlotter.h>
#include <VTKHandler.h>
#include "BonesMatOptDefines.h"
#include "BonesMatOptEnergies.h"
#include "BonesOptDeformSolver.h"



namespace shapeOptBonePolymerPeriodicHomogenization{


    
template< typename MatOptConfigurator, MaterialTypeBonePolymer MaterialType>
class BonesInterface{
    
protected:
  typedef typename MatOptConfigurator::ConfiguratorType ConfiguratorType;
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::DomVecType DomVecType;
  typedef typename ConfiguratorType::PointType PointType;
  typedef typename ConfiguratorType::VectorType VectorType;
  typedef typename ConfiguratorType::QuadRuleType QuadRuleType;
  typedef typename ConfiguratorType::DerivativeVectorValuedType DerivativeVectorValuedType;
  typedef typename ConfiguratorType::DTContainer DataTypeContainer;
  
  const MatOptConfigurator &_matOptConf;
  const Material<RealType> & _MaterialBone,  &_MaterialPolymer;
  RealType _mu, _lambda;
  const RealType _factorVoidMaterial;
  
public :
    
    BonesInterface ( const MatOptConfigurator &matOpConf ) :
      _matOptConf ( matOpConf ),
      _MaterialBone ( matOpConf._materialInfo._MaterialBone ), 
      _MaterialPolymer( matOpConf._materialInfo._MaterialPolymer ),
      _factorVoidMaterial( matOpConf._factorVoidMaterial ) {
          switch( MaterialType ){
              case BONE :
                _mu = _MaterialBone.getMu(); _lambda = _MaterialBone.getLambda();
              break;
            case POLYMER :
                _mu = _MaterialPolymer.getMu(); _lambda = _MaterialPolymer.getLambda();
                break;
            default:
            throw std::invalid_argument( aol::strprintf ( "Wrong channel. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );
            break;
         }
      }
    
    
    void evaluateStress( const string fileName, const VectorType &material, const VectorType &dispPeriodic, const VectorType &dispAffine, std::vector<RealType> &stressVec, const bool isBlock ) const{
        
       PointType lengthSingleCell; for( int i=0; i< lengthSingleCell.size(); ++i ) lengthSingleCell[i] = _matOptConf._conf.getInitializer().getWidth(i);
        
       vtkHandler vtkHan;
       std::vector<PointType> pointVec;
       vtkHan.template getPoints<PointType>( fileName, pointVec );
       stressVec.resize ( pointVec.size() );
       
       //affine part
       const QuocDiscreteFunctionDefaultAffine<DataTypeContainer, ConfiguratorType::dimDomain> dispAffineDFD ( dispAffine );
       const RealType div_Affine = dispAffineDFD.getDiv();
       const DerivativeVectorValuedType GradSym_Affine = dispAffineDFD.getSymGrad();
       
       QuocDiscreteFunctionDefault<ConfiguratorType> materialDFD ( _matOptConf._conf, material );
       QuocDiscreteVectorFunctionDefault<ConfiguratorType> dispDFD ( _matOptConf._conf, dispPeriodic );
       for( int p=0; p<pointVec.size(); ++p ){
           PointType pointVecSingleCell;
           if( isBlock ){
            for( int i=0; i< pointVecSingleCell.size(); ++i ){
                const RealType relativeCoord = pointVec[p][i] / lengthSingleCell[i];
                const int fac = static_cast<int> ( relativeCoord );
                pointVecSingleCell[i] = relativeCoord - fac;
            }
           }else{
             pointVecSingleCell = pointVec[p];   
           }
           int elementNumber; PointType LocalCoord;
           _matOptConf._conf.getLocalCoords ( pointVecSingleCell, elementNumber, LocalCoord );
           const RealType chi = _matOptConf.template approxCharFct_material<MaterialType> ( materialDFD.evaluate( _matOptConf._conf.getInitializer().getElement(elementNumber), LocalCoord )  );
           const RealType materialfactor = chi + _factorVoidMaterial * (1. - chi);
           RealType div_Periodic = dispDFD.evaluateDivergence( _matOptConf._conf.getInitializer().getElement(elementNumber), LocalCoord );
           DerivativeVectorValuedType GradSym_Periodic;
           dispDFD.evaluateSymmetrizedGradient( _matOptConf._conf.getInitializer().getElement(elementNumber), LocalCoord, GradSym_Periodic );
           DerivativeVectorValuedType stress;
           stress.setZero();
           stress += 2. * _mu * ( GradSym_Periodic + GradSym_Affine );
           for( int i=0; i< _matOptConf._conf.dimDomain; ++i ) stress(i,i) += _lambda * ( div_Periodic + div_Affine );
           
           RealType vonMisesStress = 0.;
           for( int i=0; i< _matOptConf._conf.dimDomain; ++i ){
               vonMisesStress += aol::Sqr( stress(i,i) );
               for( int j=i+1; j< _matOptConf._conf.dimDomain; ++j ){
                     vonMisesStress -= stress(i,i) * stress(j,j);
                     vonMisesStress += 3 * aol::Sqr( stress(i,j) );
               }
           }
           stressVec[p] = std::sqrt( vonMisesStress );
       }
    }


    void plotStress( const string inputFileNameVTK, const string outputFileNameVTK, const string name, const std::vector<RealType> &stressVec ) const {
        string file_extension_input  = boost::filesystem::extension(inputFileNameVTK);
        vtkSmartPointer<vtkDataSet> data;
        if( file_extension_input == ".vtp" ){
            auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(inputFileNameVTK.c_str());
            reader->Update();
            data->DeepCopy( reader->GetOutput() );
        }
        if( file_extension_input == ".vtk" ){
            
            auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
            reader->SetFileName(inputFileNameVTK.c_str());
            reader->Update();

            if( reader->IsFileUnstructuredGrid() ){
              data = vtkSmartPointer<vtkUnstructuredGrid>::New();
            }
            if( reader->IsFilePolyData() ){
              data = vtkSmartPointer<vtkPolyData>::New();
            }
            if( reader->IsFileStructuredPoints() ){
              data = vtkSmartPointer<vtkStructuredPoints>::New();   
            }
            
            data->DeepCopy(reader->GetOutput());
        }
        
        
        vtkSmartPointer<vtkDoubleArray> stressVecVTK = vtkSmartPointer<vtkDoubleArray>::New();
        stressVecVTK->SetName( name.c_str() );
        stressVecVTK->SetNumberOfTuples(data->GetNumberOfPoints());
        for( int i=0; i<stressVec.size(); ++i) stressVecVTK->SetValue( i, stressVec[i] );
        data->GetPointData()->AddArray(stressVecVTK);
        
        string file_extension_output  = boost::filesystem::extension(outputFileNameVTK);
        if( file_extension_output == ".vtp" ){
            auto vtkwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            vtkwriter->SetFileName( outputFileNameVTK.c_str() );
            vtkwriter->SetInputData(data);
            vtkwriter->SetDataModeToAscii();
            vtkwriter->Update();
            vtkwriter->Write();
        }
        if( file_extension_output == ".vtk" ){
            auto vtkwriter = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
            vtkwriter->SetFileName( outputFileNameVTK.c_str() );
            vtkwriter->SetInputData(data);
            vtkwriter->SetFileTypeToASCII();
            vtkwriter->Update();
            vtkwriter->Write();
        }
        

    }
    
};
    

}//end namespace

#endif //__PERIODICHOMOGENIZATIONBONESINTERFACE_H
