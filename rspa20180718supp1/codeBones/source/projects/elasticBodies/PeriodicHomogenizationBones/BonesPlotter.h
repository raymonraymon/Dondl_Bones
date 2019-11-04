#ifndef __PERIODICHOMOGENIZATIONBONESMULTIPLELOADPLOTTER_H
#define __PERIODICHOMOGENIZATIONBONESMULTIPLELOADPLOTTER_H

#include <general.h>
#include <VTKPlotter.h>
#include <tikzPlotter.h>
#include <SolverInfo.h>
#include <imageHandler.h>
#include "BonesMatOptDefines.h"
#include "BonesMatOptEnergies.h"
#include "BonesOptDeformSolver.h"



namespace shapeOptBonePolymerPeriodicHomogenization{


    
    
template< typename DataTypeContainer >
class VTKPlotterMaterialOptimizationBonesMultipleLoad{
private:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::PointType PointType;
  typedef typename DataTypeContainer::ParameterParserType   ParameterParserType;
  
  const ParameterParserType &_parser;
  const string _saveDirectory;
  const int _numLoads;

  
public :
    
    VTKPlotterMaterialOptimizationBonesMultipleLoad ( const ParameterParserType &parser, const int numLoads ) :
    _parser ( parser ), _numLoads( numLoads ), _saveDirectory( parser.template get<string>( "saving.saveDirectory" ).c_str() ) { }
    

   void plotAll( const string name ) const{
       
       vtkPlotter plotter( _saveDirectory );
       ParameterParserType parserVTK( "../../../../ParameterParser/shapeDesignBones/vtkPlotInterfaceCell.ini" );
       ImageHandler imHandler;
       
       for( int loadIdx=0; loadIdx< _numLoads; ++loadIdx ){
            plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf ( "%s/Deformation/%s_DeformationBone_Dir%d_Total.vtk", _saveDirectory.c_str(), name.c_str(), loadIdx ), 
                               aol::strprintf ( "%s/Deformation/%s_DeformationBone_Dir%d_Total.png", _saveDirectory.c_str(), name.c_str(), loadIdx ).c_str(),
                               parserVTK );
            plotter.template plotToPngWithParserInfo<ParameterParserType>(
                               aol::strprintf ( "%s/Deformation/%s_DeformationPolymer_Dir%d_Total.vtk", _saveDirectory.c_str(), name.c_str(), loadIdx ), 
                               aol::strprintf ( "%s/Deformation/%s_DeformationPolymer_Dir%d_Total.png", _saveDirectory.c_str(), name.c_str(), loadIdx ).c_str(),
                               parserVTK );
       }
       plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                           aol::strprintf ( "%s/%s_Undeformed.vtk", _saveDirectory.c_str(), name.c_str() ),  
                           aol::strprintf ( "%s/%s_Undeformed_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                           parserVTK );
       
       imHandler.trimAndResizeImage( aol::strprintf ( "%s/%s_Undeformed_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     aol::strprintf ( "%s/%s_Undeformed.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       
    }
    
    void plotInterface( const string name ) const{
       
       vtkPlotter plotter( _saveDirectory );
       ParameterParserType parserVTK( "../../../../ParameterParser/shapeDesignBones/vtkPlotInterfaceCell.ini" );
       ImageHandler imHandler;
       
       plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf ( "%s/Interface/%s_InterfaceBone_Subdiv2.vtk", _saveDirectory.c_str(), name.c_str() ), 
                               aol::strprintf ( "%s/Interface/%s_InterfaceBone_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                               parserVTK );
       imHandler.trimAndResizeImage( aol::strprintf ( "%s/Interface/%s_InterfaceBone_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     aol::strprintf ( "%s/Interface/%s_InterfaceBone_Subdiv2.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf ( "%s/Interface/%s_InterfacePolymer_Subdiv2.vtk", _saveDirectory.c_str(), name.c_str() ), 
                               aol::strprintf ( "%s/Interface/%s_InterfacePolymer_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                               parserVTK );
       imHandler.trimAndResizeImage( aol::strprintf ( "%s/Interface/%s_InterfacePolymer_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     aol::strprintf ( "%s/Interface/%s_InterfacePolymer_Subdiv2.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf ( "%s/Interface/%s_BlockInterfaceBone_Subdiv2.vtk", _saveDirectory.c_str(), name.c_str() ), 
                               aol::strprintf ( "%s/Interface/%s_BlockInterfaceBone_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                               parserVTK );
       imHandler.trimAndResizeImage( aol::strprintf ( "%s/Interface/%s_BlockInterfaceBone_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     aol::strprintf ( "%s/Interface/%s_BlockInterfaceBone_Subdiv2.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf ( "%s/Interface/%s_BlockInterfacePolymer_Subdiv2.vtk", _saveDirectory.c_str(), name.c_str() ), 
                               aol::strprintf ( "%s/Interface/%s_BlockInterfacePolymer_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                               parserVTK );
       imHandler.trimAndResizeImage( aol::strprintf ( "%s/Interface/%s_BlockInterfacePolymer_Subdiv2_old.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     aol::strprintf ( "%s/Interface/%s_BlockInterfacePolymer_Subdiv2.png", _saveDirectory.c_str(), name.c_str() ).c_str(),
                                     parserVTK.template get<int> ("Window.resizeWidth"), parserVTK.template get<int> ("Window.resizeHeight")  );
       
    }
    
    
    void plotVonMisesStressesOnInterface( const string name, const int numLoads ) const{
       vtkPlotter plotter( _saveDirectory );
       ParameterParserType parserVTK( "../../../../ParameterParser/shapeDesignBones/vtkPlotInterfaceCell.ini" );
       ImageHandler imHandler;
       
       ParameterParserType parserVTKStress ( parserVTK );
       parserVTKStress.set( "ScalarBar.useScalarData", parserVTK.template get<int> ( "ScalarBarStress.useScalarData" )  );
       parserVTKStress.set( "ScalarBar.useScalarBar", parserVTK.template get<int> ( "ScalarBarStress.useScalarBar" ) );
       parserVTKStress.set( "ScalarBar.ColorSpaceType", parserVTK.template get<int> ( "ScalarBarStress.ColorSpaceType" ) );
       parserVTKStress.set( "ScalarBar.ScaleType", parserVTK.template get<int> ( "ScalarBarStress.ScaleType" ) );
//        parserVTKStress.set( "Light.useLight", 0 ); //TODO e.g. for cylinder
       const int scaleRGBPointsToRange = parserVTK.template get<int> ( "ScalarBarStress.scaleRGBPointsToRange" );
       parserVTKStress.set( "ScalarBar.scaleRGBPointsToRange", scaleRGBPointsToRange );
       const int numRGBPoints = parserVTK.template get<int> ( "ScalarBarStress.numRGBPoints" );
       parserVTKStress.set( "ScalarBar.numRGBPoints", numRGBPoints );
       for( int i=1; i<=numRGBPoints; ++i ){
            parserVTKStress.set( aol::strprintf("ScalarBar.RGBValue%d", i ).c_str(), parserVTK.template get<double>( aol::strprintf("ScalarBarStress.RGBValue%d", i ).c_str() ) );
            PointType RGBPoint; parserVTK.template getFixSizeVector<double,PointType> ( aol::strprintf("ScalarBarStress.RGBPoint%d", i ).c_str(), RGBPoint );
            parserVTKStress.template setFixSizeVector<PointType>( aol::strprintf("ScalarBar.RGBPoint%d", i ).c_str(), RGBPoint, false );
       }
       
     
       for( int loadIdx=0; loadIdx<numLoads; ++loadIdx ){
           
            plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfaceBone_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfaceBone_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                               parserVTKStress, true, "VonMisesStresses", 1 );
            imHandler.trimAndResizeImage( aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfaceBone_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfaceBone_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
            plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfacePolymer_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfacePolymer_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                               parserVTKStress, true, "VonMisesStresses", 1  );
            imHandler.trimAndResizeImage( aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfacePolymer_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_InterfacePolymer_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
            plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfaceBone_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfaceBone_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                               parserVTKStress, true, "VonMisesStresses", 1  );
            imHandler.trimAndResizeImage( aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfaceBone_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfaceBone_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
            plotter.template plotToPngWithParserInfo<ParameterParserType>( 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfacePolymer_Direction%d.vtp", _saveDirectory.c_str(), loadIdx ).c_str(), 
                               aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfacePolymer_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                               parserVTKStress, true, "VonMisesStresses", 1  );
            imHandler.trimAndResizeImage( aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfacePolymer_Direction%d_old.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          aol::strprintf( "%s/StressOnInterface/SolMaterial_VonMises_BlockInterfacePolymer_Direction%d.png", _saveDirectory.c_str(), loadIdx ).c_str(),
                                          parserVTKStress.template get<int> ("Window.resizeWidth"), parserVTKStress.template get<int> ("Window.resizeHeight")  );
        }
    }
    
//    void plotAll( ) const{
//      this->plotAll( "InitMaterial" );
//      this->plotAll( "SolMaterial" );
//      this->plotInterface( "SolMaterial" );
//    }
};
    
    
    
    
    
    
    
    
    
    
template< typename DataTypeContainer >
class TikzPlotterMaterialOptimizationBonesMultipleLoad{
private:
  
  typedef typename DataTypeContainer::RealType RealType;
  typedef typename DataTypeContainer::VectorType VectorType;
  typedef typename DataTypeContainer::ParameterParserType  ParameterParserType;
  
  TikzPlotterHelperClass<RealType> _tikzHelper;
  const ParameterParserType &_parser;
  const int _precision;
  
public :
    
    TikzPlotterMaterialOptimizationBonesMultipleLoad ( const ParameterParserType &parser, const int precision = 5 ) : _parser ( parser ), _precision ( precision ) {}

    template<typename MatOptConf>
    void plotParameters( std::ofstream &out, const MatOptConf &matOpConf ) const{
       typedef typename MatOptConf::ConfiguratorType ConfiguratorType;
       
       out << endl << endl
           << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large Parameters},colbacktitle=yellow,coltitle=black]" << endl
           << "\\begin{minipage}{0.9\\textwidth}" << endl
           << "  \\begin{minipage}{0.35\\textwidth}" << endl
           << "    \\begin{tabular}{ | c | c | }" << endl
           << "       \\hline" << endl
           << "       \\multicolumn{2}{|c|}{Mesh} \\\\ \\hline" << endl
           << "       num vertices & $ " << matOpConf._conf.getInitializer().getNumVertices() << " $\\\\" << endl
           << "                    & $ ( = ";
           if( ConfiguratorType::dimDomain == 2 ) 
               out <<  matOpConf._conf.getInitializer().getNumDofs( 0 ) << " \\times " << matOpConf._conf.getInitializer().getNumDofs( 1 ) << " ) $ \\\\ \\hline" << endl;      
           if( ConfiguratorType::dimDomain == 3 ) 
               out <<  matOpConf._conf.getInitializer().getNumDofs( 0 ) << " \\times " << matOpConf._conf.getInitializer().getNumDofs( 1 ) << " \\times " << matOpConf._conf.getInitializer().getNumDofs( 2 ) << " ) $ \\\\ \\hline" << endl;   
       out << "       num elements & $ " << matOpConf._conf.getInitializer().getNumElements() << " $\\\\ \\hline" << endl
           << "    \\end{tabular}" << endl
           << "  \\end{minipage}" << endl
           << "  %" << endl
           << "  \\begin{minipage}{0.35\\textwidth}" << endl
           << "    \\begin{tabular}{ | c | c | }" << endl
           << "      \\hline" << endl
           << "      \\multicolumn{2}{|c|}{Parameters} \\\\ \\hline" << endl
           << "      $\\varepsilon$ (interface) & $ " << matOpConf._epsInterfaceLength << "$ \\\\" << endl
           << "                                 & ( $= " << matOpConf._epsFactor << "\\cdot h$ ) \\\\ \\hline" << endl
           << "      reg parameter  $\\max_\\alpha$ & " << _parser.template get<RealType> ("MaterialOptimization.regParameterMaxFunction") << " \\\\ \\hline" << endl
           << "      $c_{\\textrm{compl}}$ & $ " << matOpConf._factorComplianceCost  << " $ \\\\ \\hline" << endl
           << "      $c_{\\textrm{interface}}$ & $ " << matOpConf._factorInterfaceCost  << " $ \\\\ \\hline" << endl
           << "      $c_{\\textrm{doubleWell}}$ & $ " << matOpConf._factorDoubleWell  << " $ \\\\ \\hline" << endl
           << "    \\end{tabular}" << endl
           << "  \\end{minipage}" << endl
           << "  %" << endl
           << "  \\begin{minipage}{0.325\\textwidth}" << endl
           << "    \\begin{tabular}{ | c | c | c | }" << endl
           << "      \\hline" << endl
           << "      \\multicolumn{3}{|c|}{Material} \\\\ \\hline" << endl
           << "                  & Bone & Polymer \\\\ \\hline" << endl
           << "       $E$        & " << matOpConf._materialInfo._MaterialBone.getElastModulus() << " & "
                                     << matOpConf._materialInfo._MaterialPolymer.getElastModulus() << " \\\\ \\hline" << endl
           << "       $\\nu$     & " << matOpConf._materialInfo._MaterialBone.getPoissonRatio() << " & "
                                     << matOpConf._materialInfo._MaterialPolymer.getPoissonRatio() << " \\\\ \\hline" << endl
           << "       $\\mu$     & " << matOpConf._materialInfo._MaterialBone.getMu() << " & "
                                     << matOpConf._materialInfo._MaterialPolymer.getMu() << " \\\\ \\hline" << endl
           << "       $\\lambda$ & " << matOpConf._materialInfo._MaterialBone.getLambda() << " & "
                                     << matOpConf._materialInfo._MaterialPolymer.getLambda() << " \\\\ \\hline" << endl
           << "     \\end{tabular}" << endl
           << "   \\end{minipage}" << endl
           << "\\end{minipage}" << endl;
           
           
      //AFFINE DISPLACEMENTS

        //read loads
           const int numAffineSymGradDofs = matOpConf._conf.numAffineSymGradDofs;
           const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
           std::vector<VectorType> affineDispBone, affineDispPolymer;
           for( int i=1; i <= numLoads; ++i ){
               VectorType affineDisp ( numAffineSymGradDofs );
               _parser.template getFixSizeVector<RealType,VectorType> ( aol::strprintf("AffineDisp.Load%d", i ).c_str(), affineDisp );
               affineDispBone.push_back ( affineDisp ); affineDispPolymer.push_back ( affineDisp );
           }
           
           
        //! \todo here only for 3D case
        out << endl << "\\bigskip" << endl << endl
            << "\\begin{minipage}{\\textwidth}" << endl            
            << " \\begin{tabular}{ | c | c | c | }" << endl
            << "  \\hline" << endl
            << "  \\multicolumn{3}{|c|}{Makroskopic affine displacements $u_\\text{aff}(x) = Ax$," << endl
            << "                        $\\qquad$" << endl 
            << "                        $\\xi = \\begin{pmatrix} \\xi_1, & \\ldots, & \\xi_6 \\\\ \\end{pmatrix}" << endl
            << "                        \\Rightarrow A = \\begin{pmatrix} \\xi_1 & \\xi_4 & \\xi_5 \\\\ \\xi_4 & \\xi_2 & \\xi_6 \\\\ \\xi_5 & \\xi_6 & \\xi_3 \\\\ \\end{pmatrix} $ } \\\\ \\hline" << endl;
        for( int loadIdx = 0; loadIdx < numLoads; ++loadIdx ){ 
          out << "    Load" << loadIdx << endl
              << "    &" << endl
              << "    $\\xi^{B} = \\begin{pmatrix} ";
          for( int i=0; i<numAffineSymGradDofs; ++i ) out << (affineDispBone[loadIdx])[i] << " & ";
          out << "\\\\ \\end{pmatrix} $" << endl
              << "    &" << endl
              << "    $\\xi^{P} = \\begin{pmatrix} ";
           for( int i=0; i<numAffineSymGradDofs; ++i ) out << (affineDispPolymer[loadIdx])[i] << " & ";
           out << "\\\\ \\end{pmatrix} $" << endl
               << "    \\\\ \\hline" << endl;
        }
        out << " \\end{tabular}" << endl
            << "\\end{minipage}" << endl;
        
        
        out << endl << "\\bigskip" << endl << endl
            << "\\begin{minipage}{\\textwidth}" << endl            
            << " \\begin{tabular}{ | c | }" << endl
            << "  \\hline" << endl
            << "$ g^{B}(x) = " << matOpConf.template getWeightFunctionLoad<BONE>().description().c_str() << " $ \\\\ \\hline" << endl
            << "$ g^{P}(x) = " << matOpConf.template getWeightFunctionLoad<POLYMER>().description().c_str() << " $ \\\\ \\hline" << endl
            << " \\end{tabular}" << endl
            << "\\end{minipage}" << endl;
            
            
            
            
        out << "\\end{tcolorbox}" << endl;

    }
    
    
    void plotConvergence( std::ofstream &out, const SolverInfo<DataTypeContainer> &solverInfo ) const{
            
        out << "   \\begin{minipage}{0.35\\textwidth}" << endl
            << "     \\begin{tabular}{ | c | c | }" << endl
            << "       \\hline" << endl   
            << "       \\multicolumn{2}{|c|}{Convergence} \\\\ \\hline" << endl
            << "       solver status & " << solverInfo.getSolverStatus().c_str() << " \\\\ \\hline" << endl
            << "       error & $ " << solverInfo.getError() << " $ \\\\ \\hline" << endl
            << "       num iterations & $ " << solverInfo.getNumIterations() << " $ \\\\ \\hline" << endl
            << "     \\end{tabular}" << endl
            << "   \\end{minipage}" << endl;
    }
    
    
    template<typename MatOptConf>
    void plotResult( std::ofstream &out, const MatOptConf &matOpConf, const MaterialOptimizationMultipleLoadEnergyInfo<MatOptConf> &energyInfo,
                    const string fileName, const string infostring ) const{

        out << setprecision(_precision);
                        
        out << "\\bigskip" << endl << endl 
            << "\\begin{tcolorbox}[width=\\textwidth,colback={white},title={\\centering \\large " << infostring.c_str() << "},colbacktitle=yellow,coltitle=black]" << endl
            << "  \\begin{minipage}{\\textwidth}" << endl;
                out << "     % image " << endl
                    << "     \\begin{minipage}{0.33\\textwidth}" << endl
                    << "        \\includegraphics[width=\\textwidth]{" << fileName.c_str() << "}" << endl
                    << "     \\end{minipage}" << endl;
                out << "     % total energy " << endl
                    << "     \\begin{minipage}{0.33\\textwidth}" << endl
                    << "        \\begin{tabular}{ | c | c | }" << endl 
                    << "          \\hline" << endl
                    << "          \\multicolumn{2}{|c|}{Cost Functional} \\\\ \\hline" << endl
                    << "          Compliance & $ " << energyInfo._complianceEnergy << " $ \\\\ \\hline" << endl
                    << "          Interface & $ " << energyInfo._interfaceEnergy << " $ \\\\ \\hline" << endl
                    << "          Total & $ " << matOpConf._factorComplianceCost * energyInfo._complianceEnergy + matOpConf._factorInterfaceCost * energyInfo._interfaceEnergy << " $ \\\\ \\hline" << endl
                    << "          Volume Bone & $ " << energyInfo._volumeBone << " $ \\\\ \\hline" << endl
                    << "          Volume Polymer & $ " << energyInfo._volumePolymer << " $ \\\\ \\hline" << endl
                    << "          Barycenter & $ ( ";
//                     out << energyInfo._barycenterBone.transpose();
                    for( int i = 0; i < energyInfo._barycenterBone.size(); ++i ){
                        out << std::setprecision(2) << energyInfo._barycenterBone[i];
                        if( i < energyInfo._barycenterBone.size() - 1 ) out << ", ";
                    }
                    out << setprecision(_precision);
                    out << " ) $ \\\\ \\hline" << endl
                    << "        \\end{tabular}" << endl             
                    << "      \\end{minipage}" << endl; 
                out << "     % compliance in detail" << endl
                    << "     \\begin{minipage}{0.3\\textwidth}" << endl
                    << "       \\begin{tabular}{ | c | c | c | }" << endl
                    << "         \\hline" << endl
                    << "         & $C_{*}^B \\xi : \\xi$ & $C_{*}^P \\xi : \\xi$  \\\\ \\hline" << endl;
                    for( int i=0; i<energyInfo._numLoads; ++i ){
                        out << "         Load $" << i 
                        << " $ & $ " << energyInfo._complianceBone[i]  << " $ "
                        << "   & $ " << energyInfo._compliancePolymer[i]  << " $ "
                        << " \\\\ \\hline" << endl;
                    }
                    out << "        $ g^{\\rho}(x) = $"
                    << " &  $ " << energyInfo._complianceWeightBone << "$" 
                    << "           &  $" << energyInfo._complianceWeightPolymer << "$ \\\\ \\hline" << endl;
                    out << "       \\end{tabular}" << endl
                        << "     \\end{minipage}" << endl;
        out << "  \\end{minipage}" << endl;
        out << "\\end{tcolorbox}" << endl;

   }
    
    

   template<typename MatOptConf>
   void plotAll( const MatOptConf &matOpConf,
                 const MaterialOptimizationMultipleLoadEnergyInfo<MatOptConf> & energyInfoInit, const MaterialOptimizationMultipleLoadEnergyInfo<MatOptConf> & energyInfoSolution,
                 const SolverInfo<DataTypeContainer> &solverInfo ) const{
        
        std::ofstream out ( aol::strprintf ( "%s/Results.tex", _parser.template get<string>("saving.saveDirectory").c_str()  ).c_str ()  );
        
        _tikzHelper.generateIncludes( out, 21, 29.7, 0.5, 0.5, 0.5, 0.5 );
        _tikzHelper.generateBeginDocument( out );

        _tikzHelper.generateComment( out, "Parameters" );
        this->plotParameters<MatOptConf>( out, matOpConf );
        
        _tikzHelper.generateComment( out, "Initial Material" );
        string fileNameInit = aol::strprintf( "InitMaterial_Undeformed.png" ).c_str();
        this->plotResult<MatOptConf>( out, matOpConf, energyInfoInit, fileNameInit, "initial material" );
        
        _tikzHelper.generateComment( out, "Solution Material" );
        const string fileNameSolution = aol::strprintf( "SolMaterial_Undeformed.png" ).c_str();
        this->plotResult<MatOptConf>( out, matOpConf, energyInfoSolution, fileNameSolution, "solution material" );
        
        this->plotConvergence( out, solverInfo );
        
         _tikzHelper.generateEndDocument( out );
    }
    
    
    
    
    
    template<typename MatOptConf>
    void plotResultSingle( std::ofstream &out, const MatOptConf &matOpConf, const MaterialOptimizationMultipleLoadEnergyInfo<MatOptConf> &energyInfo,
                    const string fileName, const string infostring, const bool plotWeightFunction = true ) const{
        
        out << " \\begin{minipage}{0.75\\textwidth}" << endl
            << "        \\includegraphics[width=\\textwidth]{" << fileName.c_str() << "}" << endl
            << " \\end{minipage}" << endl;
          
        out << endl << "  \\bigskip" << endl << endl;
            
        out << " \\begin{minipage}{0.75\\textwidth}" << endl
        << "        \\begin{tabular}{ | c | c | }" << endl 
        << "          \\hline" << endl
        << "          \\multicolumn{2}{|c|}{Cost Functional} \\\\ \\hline" << endl
        << "          Compliance & $ " << energyInfo._complianceEnergy << " $ \\\\ \\hline" << endl
        << "          Interface & $ " << energyInfo._interfaceEnergy << " $ \\\\ \\hline" << endl
        << "          Total & $ " << matOpConf._factorComplianceCost * energyInfo._complianceEnergy + matOpConf._factorInterfaceCost * energyInfo._interfaceEnergy << " $ \\\\ \\hline" << endl
        << "        \\end{tabular}" << endl
        << "      \\end{minipage}" << endl; 
                     
        out << endl << "  \\bigskip" << endl << endl;
        
        //Compliance
        out << "     \\begin{minipage}{0.9\\textwidth}" << endl
            << "       \\begin{tabular}{ | c | c | c | }" << endl
            << "         \\hline" << endl
            << "         & $C_{*}^B \\xi : \\xi$ & $C_{*}^P \\xi : \\xi$  \\\\ \\hline" << endl;
        for( int i=0; i<energyInfo._numLoads; ++i ){
        out << "         Load $" << i << " $ & $ " << energyInfo._complianceBone[i]  << " $ "
                                      << " & $ " << energyInfo._compliancePolymer[i]  << " $ "
                                      << " \\\\ \\hline" << endl;
        }
        if( plotWeightFunction ){
            out << "        $ g^{\\rho}(x) = "
                << " $ & $ " << energyInfo._complianceWeightBone << "$" 
                << "           &  $" << energyInfo._complianceWeightPolymer << "$ \\\\ \\hline" << endl;
        }
        out << "       \\end{tabular}" << endl
            << "     \\end{minipage}" << endl
            << "     %" << endl;
   }
    
    
   template<typename MatOptConf>
   void plotInit( const MatOptConf &matOpConf,
                  const MaterialOptimizationMultipleLoadEnergyInfo<MatOptConf> & energyInfoInit ) const{
        
        std::ofstream out ( aol::strprintf ( "%s/ResultsInit.tex", _parser.template get<string>("saving.saveDirectory").c_str()  ).c_str ()  );
        
        const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
        
        _tikzHelper.generateIncludes( out, 6, 8 + static_cast<RealType> ( numLoads ) / 2., 0.0, 0.0, 0.0, 0.0 );
        _tikzHelper.generateBeginDocument( out );
        
        _tikzHelper.generateComment( out, "Initial Material" );
        string fileNameInit = aol::strprintf( "InitMaterial_Undeformed.png" ).c_str();
        this->plotResultSingle<MatOptConf>( out, matOpConf, energyInfoInit, fileNameInit, "initial material", false );
        
         _tikzHelper.generateEndDocument( out );
    }
    
    
    
    void plotResults_CompareDesigns ( const string fileName, const int numDesigns ) const {
        
        std::ofstream out ( aol::strprintf ( "%s/%s.tex", _parser.template get<string>("saving.saveDirectory").c_str (), fileName.c_str ()  ) );
        
        const int numLoads = _parser.template get<int> ( "AffineDisp.numLoads" );
        
        _tikzHelper.generateIncludes( out, 10 + 5 * numDesigns, 16 + numLoads, 0.5, 0.5, 0.5, 0.5 );
        out << "\\usepackage{graphicx}" << endl
            << "\\usepackage{array}" << endl
            << "\\usepackage{colortbl}" << endl;
            
        _tikzHelper.generateBeginDocument( out );

        
        //plot Parameters (use from Design0)
        _tikzHelper.generateComment( out, "Parameters" );
        out << "\\begin{minipage}{0.9\\textwidth}" << endl
            << " \\include{Design-0/Parameter_0}" << endl
            << " \\end{minipage}" << endl << endl << endl;
        
        _tikzHelper.generateComment( out, "Initial Material" );
        out << "\\begin{minipage}{\\textwidth}" << endl
            << "\\begin{tabular}{ "; for(int i=0; i<numDesigns; ++i ) out << " | c ";
        out << "| }" << endl
            << "  \\hline" << endl;
        for( int i=0; i<numDesigns; ++i ){
             if( i > 0 ) out << "  &" << endl;
             string fileNamePDF = "Design-" + aol::strprintf( "%d", i ) + "/ResultsInit_0.pdf"; 
             out << "  \\begin{minipage}{" << 1. / static_cast<RealType> ( numDesigns + 2 ) << "\\textwidth}" << endl
                 << "      \\includegraphics[width=\\linewidth, height=10cm]{" << fileNamePDF.c_str() << "}" << endl
                 << "  \\end{minipage}" << endl;
        }
        out << "  \\\\ \\hline" << endl;
        
        out << "\\end{tabular}" << endl 
            << "\\end{minipage}" << endl;
            
     _tikzHelper.generateEndDocument( out );

   }
   
   template<typename MatOptConf>
   void generateTexCodeForPaper ( const string fileName, const MaterialOptimizationMultipleLoadEnergyInfo<MatOptConf> &energyInfo ) const {
       
       std::ofstream out ( aol::strprintf ( "%s/%s.tex", _parser.template get<string>("saving.saveDirectory").c_str (), fileName.c_str ()  ) );
       out << setprecision(_precision);
       
       out << "\\begin{figure}[!htbp]" << endl
           << "  \\begin{tabular}{ c  c  c }" << endl
           << "     \\begin{minipage}{0.25\\textwidth}" << endl
           << "        \\begin{tabular}{ | c | c | c | }" << endl
           << "           \\hline" << endl
           << "              & B & P \\\\ \\hline" << endl;
           
       for( int i=0; i<energyInfo._numLoads; ++i ){
           out << "         Load $" << i + 1 
           << " $ & $ " << energyInfo._complianceBone[i]  << " $ "
           << "   & $ " << energyInfo._compliancePolymer[i]  << " $ "
           << " \\\\ \\hline" << endl;
       }
       out << "         vol & " << energyInfo._volumeBone << " & " << energyInfo._volumePolymer << " \\\\ \\hline" << endl;
       
       out << "        \\end{tabular} " << endl
           << "      \\end{minipage}" << endl;
           
       out << "  & " << endl
           << "  \\begin{minipage}{0.25\\textwidth} " << endl
           << "    {\\includegraphics[width=0.9\\textwidth]{images/3DCompression_Cell_Bone.png}}" << endl
           << "  \\end{minipage}" << endl;
           
       out << "  & " << endl
           << "  \\begin{minipage}{0.25\\textwidth} " << endl
           << "    {\\includegraphics[width=0.9\\textwidth]{images/3DCompression_Block_Bone.png}}" << endl
           << "  \\end{minipage}" << endl;

       out << " \\\\ " << endl
           << "\\end{tabular}" << endl;
       
   }

};

}//end namespace

#endif //__PERIODICHOMOGENIZATIONBONESMULTIPLELOADPLOTTER_H
