#ifndef __VTKPLOTTER_H
#define __VTKPLOTTER_H

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataReader.h>
#include <vtkDataSetMapper.h>
#include <vtkFloatArray.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkLight.h>
#include <vtkLightKit.h>
#include <vtkLookupTable.h>
#include <vtkPiecewiseFunction.h>
#include <vtkPLYReader.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPostScriptWriter.h>
#include <vtkProperty.h>
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


#include <iostream>
#include <string>
 
class vtkPlotter {
    
public :
    
    vtkPlotter( const string saveDirectory ) {}
    
    vtkPlotter( ) {}
    
    template<typename ParameterParserType>
    void plotToPngWithParserInfo( const string &inputFileName, const string &outputFileName, const ParameterParserType &parser ){
             
      // Get all data from the file
        auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(inputFileName.c_str());
        reader->Update();

      // Create scalarbar and a lookup table to share between the mapper and the scalarbar
        vtkSmartPointer<vtkDataSet> grid;
        if( reader->IsFileUnstructuredGrid() ) grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        if( reader->IsFilePolyData() ) grid = vtkSmartPointer<vtkPolyData>::New();
        if( reader->IsFileStructuredPoints() )  grid = vtkSmartPointer<vtkStructuredPoints>::New();   
        //if( IsFileStructuredGrid() )  grid = vtkSmartPointer<vtkStructuredGrid>::New();
        //if( IsFileRectilinearGrid() )   grid = vtkSmartPointer<vtkRectilinearGrid>::New();

        grid->DeepCopy(reader->GetOutput());
        
        auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
        mapper->SetInputData(grid);
        
      //actor
        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetSpecular( parser.template get<double> ( "Actor.Specular" ) ); 
        actor->GetProperty()->SetSpecularPower( parser.template get<double> ( "Actor.SpecularPower" ) ); 
        double DiffuseColor[3]; parser.template getFixSizeVector<double,double[3]> ( "Actor.DiffuseColor", DiffuseColor );
        actor->GetProperty()->SetDiffuseColor( DiffuseColor ); 
        
      //light
         auto light = vtkSmartPointer<vtkLightKit>::New();
        //light->SetPosition(lightPosition);         
        light->SetKeyLightWarmth(parser.template get<double>("Light.KeyWarmth")); light->SetKeyLightIntensity(parser.template get<double>("Light.KeyIntensity"));
        light->SetKeyLightAngle(parser.template get<double>("Light.KeyElevation"),parser.template get<double>("Light.KeyAzimuth"));
        light->SetFillLightWarmth(parser.template get<double>("Light.FillWarmth")); light->SetKeyToFillRatio(parser.template get<double>("Light.FillKeyRatio"));
        light->SetFillLightAngle(parser.template get<double>("Light.FillElevation"),parser.template get<double>("Light.FillAzimuth"));
        light->SetBackLightWarmth(parser.template get<double>("Light.BackWarmth")); light->SetKeyToBackRatio(parser.template get<double>("Light.BackKeyRatio"));
        light->SetBackLightAngle(parser.template get<double>("Light.BackElevation"),parser.template get<double>("Light.BackAzimuth"));
        light->SetHeadLightWarmth(parser.template get<double>("Light.HeadWarmth")); light->SetKeyToHeadRatio(parser.template get<double>("Light.HeadKeyRatio"));
        
      //scalar bar
        auto scalarBar = vtkSmartPointer<vtkScalarBarActor>::New(); 
        if( parser.template get<bool>("ScalarBar.useScalarData") ){
            
             // color transfer
            auto ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
            
             const int numRGBPoints = parser.template get<int> ( "ScalarBar.numRGBPoints" );
             for( int i=1; i<=numRGBPoints; ++i ){
                 double RGBValue = parser.template get<double> ( aol::strprintf("ScalarBar.RGBValue%d", i ).c_str() );
                 double RGBPoint[3]; parser.template getFixSizeVector<double,double[3]> ( aol::strprintf("ScalarBar.RGBPoint%d", i ).c_str(), RGBPoint );
                 ctf->AddRGBPoint(RGBValue, RGBPoint[0], RGBPoint[1], RGBPoint[2]);
             }
             
             switch( parser.template get<int> ( "ScalarBar.ColorSpaceType" ) ){
                 case 1 : ctf->SetColorSpaceToDiverging(); break;
                 case 2 : ctf->SetColorSpaceToRGB(); break;
                 case 3 : { 
                     ctf->SetColorSpaceToHSV(); 
                     ctf->HSVWrapOff(); //ctf->HSVWrapOn();
                 }break;
                 case 4 : ctf->SetColorSpaceToLab(); break;
                 default : throw invalid_argument( aol::strprintf ( "Unknown ColorSpaceType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
             }
            
             
             switch( parser.template get<int> ( "ScalarBar.ScaleType" ) ){
                 case 1 : ctf->SetScaleToLinear(); break;
                 case 2 : ctf->SetScaleToLog10(); break;
                 default : throw invalid_argument( aol::strprintf ( "Unknown ScaleType. In File %s at line %d.", __FILE__, __LINE__ ).c_str() ); 
             }
             
            mapper->SetLookupTable( ctf );
            
            scalarBar->SetLookupTable( vtkLookupTable::SafeDownCast(mapper->GetLookupTable()) );
            scalarBar->SetTitle("");
            scalarBar->SetNumberOfLabels(4); 
            scalarBar->SetLookupTable( ctf );
            
        }
        
        
      //box at boundary
        auto readerBox = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        vtkSmartPointer<vtkDataSet> gridBox;
        auto mapperBox = vtkSmartPointer<vtkDataSetMapper>::New();
        auto actorBox = vtkSmartPointer<vtkActor>::New();
        if( parser.template get<int> ("BoundaryBox.useBoundaryBox") ){
            readerBox->SetFileName( parser.template get<string>("BoundaryBox.fileName").c_str() );
            readerBox->Update();
            if( reader->IsFileUnstructuredGrid() ) gridBox = vtkSmartPointer<vtkUnstructuredGrid>::New();
            if( reader->IsFilePolyData() ) gridBox = vtkSmartPointer<vtkPolyData>::New();
            if( reader->IsFileStructuredPoints() )  gridBox = vtkSmartPointer<vtkStructuredPoints>::New();   
            gridBox->DeepCopy(readerBox->GetOutput()); 
            mapperBox->SetInputData(gridBox);
            actorBox->SetMapper(mapperBox);
        }

      //Renderer
        auto renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor);
        if( parser.template get<int> ("BoundaryBox.useBoundaryBox") ) renderer->AddActor(actorBox);
        light->AddLightsToRenderer(renderer);//for vtklight: renderer->AddLight(light);
        if( parser.template get<bool>("ScalarBar.useScalarBar")  ) renderer->AddActor2D(scalarBar);
        
        double Position[3]; parser.template getFixSizeVector<double,double[3]> ( "Camera.Position", Position );
        renderer->GetActiveCamera()->SetPosition(Position);
        double FocalPoint[3]; parser.template getFixSizeVector<double,double[3]> ( "Camera.FocalPoint", FocalPoint );
        renderer->GetActiveCamera()->SetFocalPoint(FocalPoint);
        double ViewUp[3]; parser.template getFixSizeVector<double,double[3]> ( "Camera.ViewUp", ViewUp );
        renderer->GetActiveCamera()->SetViewUp(ViewUp);
        renderer->GetActiveCamera()->SetParallelScale( parser.template get<double> ( "Camera.ParallelScale" ) );
        renderer->GetActiveCamera()->Roll( parser.template get<double> ( "Camera.Roll" ) );
        renderer->GetActiveCamera()->Azimuth( parser.template get<double> ( "Camera.Azimuth" ) );
        renderer->GetActiveCamera()->Elevation( parser.template get<double> ( "Camera.Elevation" ) );

        renderer->ResetCamera();
        renderer->GetActiveCamera()->Dolly( parser.template get<double> ( "Camera.Dolly" ) );
        renderer->ResetCameraClippingRange();
        double Background[3]; parser.template getFixSizeVector<double,double[3]> ( "Renderer.Background", Background );
        renderer->SetBackground(Background);
        renderer->UseFXAAOn(); //anti aliasing, renderWindow->SetAAFrames(10) is deprecated
        auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        
        auto renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);
        if( parser.template get<int>( "Window.FullScreenRendering" ) == 1 ) renderWindow->FullScreenOn();
        else renderWindow->SetSize( parser.template get<int>("Window.imageWidth"), parser.template get<int>("Window.imageHeight") );
        
        
        renderWindow->Render();
        // renderWindowInteractor->Start();
        renderWindow->OffScreenRenderingOn();
        renderWindow->Render();
        
      // Screenshot  
        auto windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);         
        
        auto writer = vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName( outputFileName.c_str() );
        writer->SetCompressionLevel( parser.template get<int>("Writer.compressionLevel") ); //from 0-9: 0 no compression
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
    }
    
};

#endif
