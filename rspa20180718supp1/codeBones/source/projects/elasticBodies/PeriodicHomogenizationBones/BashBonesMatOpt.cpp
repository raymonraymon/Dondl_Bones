#include <general.h>
#include <loadAndSave.h>
 
typedef double RealType;
typedef aol::BoostParser ParameterParserType;

//--------------------------------------------------------
int main(int argc, char ** argv) {
    
  cout << "Usage of programm :" << endl
         << "  1 - renewKerberosTicket (default false)" << endl
         << "  2 - useNice (default false)" << endl 
         << "  3 - dimension ( 2 or 3(default) )"
         << endl;
    
    //!
    bool renewKerberosTicket = false;
    if( argc > 1 ){
        if( stoi(argv[1]) != 0 ) renewKerberosTicket = true;
    };
    bool useNice = false;
    if( argc > 2 ){
        if( stoi(argv[2]) != 0 ) useNice = true;
    }
    
    int dimension = 3;
    if( argc > 3 ){
        dimension = stoi(argv[3]);
    }
    
    aol::BoostParser parser;
    if( dimension == 2 )
        parser = aol::BoostParser ( "../../../../ParameterParser/shapeDesignBones/BonesAffinePeriodic2D.ini", "../../../../ParameterParser/counter.txt", "/MaterialOptBones2D" );
    if( dimension == 3 )
        parser = aol::BoostParser ( "../../../../ParameterParser/shapeDesignBones/BonesAffinePeriodic3D.ini", "../../../../ParameterParser/counter.txt", "/MaterialOptBones3D" );
    const string saveDirectory = parser.template get<string> ("saving.saveDirectory" );
    cout << "saveDirectory is " << saveDirectory.c_str() << endl;    
    
    aol::consoleOutput( "start programm BonesMatOpt" );
    const string numThreads = aol::strprintf( "%d", parser.template get<int>( "BASH.numThreads" ) );
    cout << "using " << numThreads.c_str() << " threads for each design" << endl;
    string systemCommand = "OMP_NUM_THREADS=" + numThreads + " nohup ";
    if( renewKerberosTicket ) systemCommand += "krenew -- ";
    if( useNice ) systemCommand += "nice ";
    systemCommand += "xvfb-run -a --server-args \"-screen 0 1920x1080x24\" ";
    if( dimension == 2 )
        systemCommand += "./BonesMatOpt2D ";
    if( dimension == 3 )
        systemCommand += "./BonesMatOpt3D ";
    systemCommand += aol::strprintf("%s/ParameterParser.ini", saveDirectory.c_str() ) + " " + saveDirectory;
    systemCommand += " >" + saveDirectory + "/ConsoleOutput.out 2>&1 "; //! \note "2>&1" is to avoid warning nohup: redirecting stderr to stdout
    cout << "systemCommand = " << systemCommand << endl;
    bool failedMatOpt = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
    if ( failedMatOpt ) cerr << "programm returned an error." << endl;
    cout << "finished computation" << endl;
 
  return 0;
}
