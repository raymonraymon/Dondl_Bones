#include <general.h>
#include <loadAndSave.h>

#include <quocDefines.h>
#include <quocConfigurators.h>
#include <quocHandler.h>

#include "BonesMatOptSolver.h"
#include "BonesPlotter.h"


using namespace quocFE;
using namespace shapeOptBonePolymerPeriodicHomogenization;


typedef double RealType;
typedef aol::BoostParser ParameterParserType;
typedef Quoc3DDataTypeContainer  DataTypeContainer;


//--------------------------------------------------------
int main(int argc, char ** argv) {
    
    aol::BoostParser parser;
    parser = aol::BoostParser ( "../../../../ParameterParser/shapeDesignBones/BonesAffinePeriodic.ini", "../../../../ParameterParser/counter.txt", "/CompareInitiatlizationMaterialOptBones" );
    const string saveDirectory = parser.template get<string> ("saving.saveDirectory" );
    cout << "saveDirectory is " << saveDirectory.c_str() << endl;    

    //option: optimize or only init
    bool optimize = false;
    if( parser.template get<int>("CompareDesigns.optimize") == 1 ) optimize = true;
    
    const int numDesignsToCompare = parser.template get<int>( "CompareDesigns.numDesignsToCompare" );
    
    const string numThreads = aol::strprintf( "%d", parser.template get<int>( "BASH.numThreads" ) );
    cout << "using " << numThreads.c_str() << " threads for each design" << endl;
    
    for( int designType = 0; designType < numDesignsToCompare; designType++){ 
    
       cout << endl << endl << "===================================================" << endl << "designType " << designType << endl;
        
       const string saveDirectoryDesign = aol::strprintf( "%s/Design-%d", saveDirectory.c_str(), designType );
       boost::filesystem::create_directory ( saveDirectoryDesign ); 
       
       string systemCommand = "OMP_NUM_THREADS=" + numThreads + " nohup nice xvfb-run -a --server-args \"-screen 0 1920x1080x24\" ";
       systemCommand += "./BonesMatOpt " + aol::strprintf("%s/ParameterParser.ini", saveDirectory.c_str() ) + " " + saveDirectoryDesign; 
       if( !optimize ){
           systemCommand += " 1";
       }else{
           systemCommand += " 0";
       }
       
       if( parser.template get<int> ("CompareDesigns.readDesignsFromFiles") == 1 ){
           string fileName = parser.template get<string> ( aol::strprintf("CompareDesigns.materialFile%d", designType + 1 ) );
           systemCommand += aol::strprintf( " 0 %s", fileName.c_str() );
       }else{
           systemCommand += aol::strprintf( " %d", designType );
       }
       
       
       systemCommand += " >" + saveDirectoryDesign + "/ConsoleOutput.out 2>&1 ";
       cout << "systemCommand = " << systemCommand << endl;
       bool failedMatOptDesign;
       failedMatOptDesign = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
       if ( failedMatOptDesign ) cerr << "programm returned an error." << endl;
       
       int InitOrSolution = 1; //init 1, sol 2, both 3
       if( optimize ) InitOrSolution = 3;
       const int numBlocksToPlot = parser.template get<int> ("saving.numBlocks");
       const int plotInterface = parser.template get<int> ("saving.plotInterface");
       string systemCommandPlot = "nohup nice xvfb-run -a --server-args \"-screen 0 1920x1080x24\" ";
       systemCommandPlot += "./BonesPlotResults " + saveDirectory;
       systemCommandPlot += " " + std::to_string(InitOrSolution);
       systemCommandPlot += " " + std::to_string(numBlocksToPlot);
       systemCommandPlot += " " + std::to_string(plotInterface);
       systemCommandPlot += " >" + saveDirectoryDesign + "/ConsoleOutputPrint.out 2>&1 ";
       cout << "systemCommandPlot = " << systemCommandPlot << endl;
       bool failedPlot = ( system ( systemCommandPlot.c_str() ) != EXIT_SUCCESS );
       if ( failedPlot ) cerr << "plot returned an error." << endl;
       
       std::ofstream BashFilePDFDesign ( aol::strprintf ( "%s/%s.sh", saveDirectoryDesign.c_str (), "generatePDF"  ) );
       BashFilePDFDesign << "cd " << saveDirectoryDesign.c_str() << endl;
       BashFilePDFDesign << "pdflatex ResultsInit.tex" << endl;
       string systemCommandPDF = "bash " + saveDirectoryDesign + "/generatePDF.sh";
       cout << "systemCommandPDF = " << systemCommandPDF << endl;
       bool failedPDF;
       failedPDF = ( system ( systemCommandPDF.c_str() ) != EXIT_SUCCESS );
       if ( failedPDF ) cerr << "pdflatex returned an error." << endl;

    }
    
    TikzPlotterMaterialOptimizationBonesMultipleLoad<DataTypeContainer> tikzPlotter ( parser );
    tikzPlotter.plotResults_CompareDesigns ( "ResultsAll", numDesignsToCompare );
    std::ofstream BashFilePDF ( aol::strprintf ( "%s/%s.sh", saveDirectory.c_str (), "generatePDF"  ) );
    BashFilePDF << "cd " << saveDirectory.c_str() << endl;
    BashFilePDF << "pdflatex ResultsAll.tex" << endl;
    string systemCommandPDF = "bash " + saveDirectory + "/generatePDF.sh";
    cout << "systemCommandPDFAll = " << systemCommandPDF << endl;
    bool failedPDF;
    failedPDF = ( system ( systemCommandPDF.c_str() ) != EXIT_SUCCESS );
    if ( failedPDF ) cerr << "pdflatex returned an error." << endl;

  return 0;
}
