#ifndef __VTKINTERFACETOSURFACE_H
#define __VTKINTERFACETOSURFACE_H

#include <iostream>
#include <string>


#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkLookupTable.h>
#include <vtkMarchingCubes.h>
#include <vtkMath.h>
#include <vtkPLYReader.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkPostScriptWriter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredPointsReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridVolumeMapper.h>
#include <vtkVersion.h>
#include <vtkVoxelContoursToSurfaceFilter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkLoopSubdivisionFilter.h>
#include <vtkTriangleFilter.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>



 
class vtkInterfaceToSurface {
    
public :
    vtkInterfaceToSurface( ) {}
    
    //! \note  inputFile has to be vtkImageData (ie. either structuredpoint or uniformgrid)
    //         outputFile is polydata (ie. either .vtk with polydata or .vtp)
    void getSurface( const string &inputFileName, const string &outputFileName, const float threshold = 0.5, const bool extractLargest = false ){
            
        auto reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
        reader->SetFileName(inputFileName.c_str());
        reader->Update();
         
        // using marching cubes algorithm
        auto mc = vtkSmartPointer<vtkMarchingCubes>::New();
        mc->SetInputConnection(reader->GetOutputPort());
        mc->ComputeNormalsOn();
//         mc->ComputeGradientsOn();
        mc->ComputeScalarsOff();
        mc->SetValue(0, threshold);
            
        // To remain largest region
        auto confilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        confilter->SetInputConnection(mc->GetOutputPort());
        confilter->SetExtractionModeToLargestRegion();

        //save 
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New(); //.vtk
        writer->SetFileName(outputFileName.c_str());
        writer->SetFileTypeToASCII();
        if (extractLargest) writer->SetInputConnection(confilter->GetOutputPort());
        else                writer->SetInputConnection(mc->GetOutputPort());
        writer->Write();

    }
    
    
    
    //! \note inputfile has to be of type polydata 
    //! \todo check if input file is polydata
    void loopSubdivision( const string inputFileNameVTK, const string outputFileNameVTK, const int numberOfSubdivisions ) const{
        
        string file_extension_input  = boost::filesystem::extension(inputFileNameVTK);

        vtkSmartPointer<vtkPolyData> originalMesh;
        if( file_extension_input == ".vtp" ){
            auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
            reader->SetFileName(inputFileNameVTK.c_str());
            reader->Update();
            // Subdivision filters only work on triangles
            auto triangles =  vtkSmartPointer<vtkTriangleFilter>::New();
            triangles->SetInputConnection(reader->GetOutputPort());
            triangles->Update();
            originalMesh = triangles->GetOutput();
        }
        if( file_extension_input == ".vtk" ){
               auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
               reader->SetFileName(inputFileNameVTK.c_str());
               reader->Update();
               
               // Subdivision filters only work on triangles
               auto triangles =  vtkSmartPointer<vtkTriangleFilter>::New();
               triangles->SetInputConnection(reader->GetOutputPort());
               triangles->Update();
               originalMesh = triangles->GetOutput();
        }
    
        auto subdivisionFilter = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
        subdivisionFilter->SetNumberOfSubdivisions(numberOfSubdivisions);
        subdivisionFilter->SetInputData(originalMesh);
        subdivisionFilter->Update();   
        
        string file_extension_output  = boost::filesystem::extension(outputFileNameVTK);
        cout << "file_extension_output = " << file_extension_output.c_str() << endl;
        if( file_extension_output == ".vtp" ){
            auto vtkwriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
            vtkwriter->SetFileName( outputFileNameVTK.c_str() );
            vtkwriter->SetInputData(subdivisionFilter->GetOutput());
            vtkwriter->SetDataModeToAscii(); //vtkwriter->SetDataModeToBinary();
            vtkwriter->Update();
            vtkwriter->Write();
        }
        if( file_extension_output == ".vtk" ){
            auto vtkwriter = vtkSmartPointer<vtkPolyDataWriter>::New();
            vtkwriter->SetFileName( outputFileNameVTK.c_str() );
            vtkwriter->SetInputData(subdivisionFilter->GetOutput());
            vtkwriter->SetFileTypeToASCII();
            vtkwriter->Update();
            vtkwriter->Write();
        }
        
    }
    
};


#endif
