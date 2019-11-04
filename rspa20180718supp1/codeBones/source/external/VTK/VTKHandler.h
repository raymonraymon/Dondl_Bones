#ifndef __VTKHANDLER_H
#define __VTKHANDLER_H

#include <general.h>
#include <loadAndSave.h>


#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkImageData.h>
#include <vtkLookupTable.h>
#include <vtkPLYReader.h>
#include <vtkPNGWriter.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkRectilinearGrid.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridVolumeMapper.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXMLCompositeDataReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

 
class vtkHandler {
    
public :
    vtkHandler( ) {}
    
    template<class TReader> 
    void ReadAnXMLFile(const string &inputFileName, vtkSmartPointer<vtkDataSet> &dataSet ) const{
        auto reader = vtkSmartPointer<TReader>::New();
        reader->SetFileName(inputFileName.c_str());
        reader->Update();
        reader->GetOutput()->Register(reader);
        dataSet = vtkDataSet::SafeDownCast(reader->GetOutput());
    }
    
    void getVTKDataSet( const string &inputFileName, vtkSmartPointer<vtkDataSet> &dataSet ) const{
        string extension  = boost::filesystem::extension(inputFileName);
        if      (extension == ".vtu") ReadAnXMLFile<vtkXMLUnstructuredGridReader> (inputFileName, dataSet );
        else if (extension == ".vtp") ReadAnXMLFile<vtkXMLPolyDataReader> (inputFileName, dataSet );
        else if (extension == ".vts") ReadAnXMLFile<vtkXMLStructuredGridReader> (inputFileName, dataSet );
        else if (extension == ".vtr") ReadAnXMLFile<vtkXMLRectilinearGridReader> (inputFileName, dataSet );
        else if (extension == ".vti") ReadAnXMLFile<vtkXMLImageDataReader> (inputFileName, dataSet );
        else if (extension == ".vtk") ReadAnXMLFile<vtkDataSetReader> (inputFileName, dataSet );
        else throw invalid_argument( aol::strprintf ( "Unknown VTK-FileExtension. In File %s at line %d.", __FILE__, __LINE__ ).c_str() );    
     }
    
    
    template<typename PointType>
    void getPoints( const string &inputFileName, std::vector<PointType> &pointVec ) const {
        
        string file_extension  = boost::filesystem::extension(inputFileName);
        
        // Get all data from the file
        vtkSmartPointer<vtkDataSet> data;
        this->getVTKDataSet( inputFileName, data );
        
        vtkIdType idNumPointsInFile = data->GetNumberOfPoints();
  
        PointType dummyPoint; const int dim = dummyPoint.size(); 
        for(int i = 0; i < idNumPointsInFile; i++){
            double point[dim];
            data->GetPoint(i, point);
            PointType p; for( int j=0;j<dim; ++j ) p[j] = point[j];
            pointVec.push_back( p );
        }
        
    }
    
    
    template<typename Indices3DType>
    void getCells( const string &inputFileName, std::vector<Indices3DType> &triangleVec ) const {
        
        string file_extension  = boost::filesystem::extension(inputFileName);
        
        // Get all data from the file
        vtkSmartPointer<vtkDataSet> data;
        this->getVTKDataSet( inputFileName, data );
        
        vtkIdType idNumCellsInFile = data->GetNumberOfCells();
  
        for(vtkIdType i = 0; i < idNumCellsInFile; i++){
            vtkCell* cell = data->GetCell(i);
            
            vtkIdList *listOfPointIdx = cell->GetPointIds();
            Indices3DType trianglePointIdx;
            for( int i=0; i < listOfPointIdx->GetNumberOfIds(); ++i ){
                trianglePointIdx[i] = listOfPointIdx->GetId(i);
            }
            triangleVec.push_back( trianglePointIdx );
        }
        
    }
    
  
    //Attention: does not work for legacy file format, only detects first array
    template<typename VecType>
    void getPointDataVec ( const string &inputFileName, const string &arrayName, std::vector<VecType> &dataVec ) const {
        
        string file_extension  = boost::filesystem::extension(inputFileName);
        
        // Get all data from the file
        vtkSmartPointer<vtkDataSet> dataSet;
        this->getVTKDataSet( inputFileName, dataSet );
        
        vtkIdType idNumPointsInFile = dataSet->GetNumberOfPoints(); 
        dataVec.resize( idNumPointsInFile );

        vtkSmartPointer<vtkPointData> pd = dataSet->GetPointData();
        vtkSmartPointer<vtkDataArray> dataArray = dataSet->GetPointData()->GetArray(arrayName.c_str());
            if(dataArray){
                    vtkIdType numVectors = dataArray->GetNumberOfTuples();
                    int numComponents = dataArray->GetNumberOfComponents();
                    if( numVectors > 1 ){
                        for (vtkIdType tupleIdx = 0; tupleIdx < numVectors; ++tupleIdx)
                            for( int compIdx = 0; compIdx < numComponents; ++compIdx ){
                                dataVec[tupleIdx][compIdx] =  dataArray->GetComponent(tupleIdx, compIdx);         
                            }
                    }
            }
        
  }

};


#endif
