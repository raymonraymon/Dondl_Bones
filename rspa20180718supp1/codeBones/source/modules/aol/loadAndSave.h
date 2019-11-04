#ifndef __GENERALLOADANDSAVE_H
#define __GENERALLOADANDSAVE_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <general.h>

#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>


namespace aol {
    
 //Output
template<typename VectorType>
void printVector( const VectorType &vec, const int breakNumber = 5, const std::string name = "", const unsigned precision = 5 ){

  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << "-----------" << name.c_str() << "-----------------------" << endl;

  int iter = 0;
  cout.precision(precision);
  for ( int i = 0; i < vec.size(); ++i ){
      if( std::abs(vec[i]) < 1.e-16 ) cout << "0" << " ; \t "; 
      else cout << vec[i] << " ; \t ";
      if( iter == breakNumber - 1 ){
        cout << endl;
        iter = 0;
      }else iter++;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
  cout.precision(5);
}

template<typename VectorType>
void printVectorToFile( const VectorType &vec, const string fileName, const unsigned precision = 5 ){
    
  std::ofstream out ( fileName.c_str() );

  out.precision(precision);
  for ( int i = 0; i < vec.size(); ++i ){
      if( std::abs(vec[i]) < 1.e-16 ) out << "0" << endl; 
      else out << vec[i] << endl;
  }
  out.close();
}


template<typename VectorType>
void loadVectorFromFile(  VectorType &vec, const string fileName ){
  std::ifstream File; 
  File.open ( fileName.c_str() );
  if ( File.fail() ) throw std::invalid_argument( aol::strprintf ( "Cannot open file with name %s. In File %s at line %d.", fileName.c_str(), __FILE__, __LINE__ ).c_str() );
  int iter = 0;
  while ( !File.eof() && iter < vec.size() ) {
      char line[256];
      File.getline ( line, 256 );  
      char * substring = strtok ( line, " " );
      vec[iter] = atof( substring );
      iter++;
  }
  File.close();
}


template<typename SparseMatrixType>
void printSparseMatrixToFile( const SparseMatrixType &mat, const string fileName, const unsigned precision = 5 ){

  std::ofstream out ( fileName.c_str() );
  out.precision(precision);
  
  for (int k=0; k<mat.outerSize(); ++k)
     for (typename SparseMatrixType::InnerIterator it(mat,k); it; ++it){
         out << it.row() << " " << it.col() << " " << it.value() << endl;
  }
  
  out.close();
}


template<typename TripletType, typename SparseMatrixType>
void loadSparseMatrixFromFile(  SparseMatrixType &mat, const string fileName ){
  std::ifstream File; 
  File.open ( fileName.c_str() );
  if ( File.fail() ) throw std::invalid_argument( aol::strprintf ( "Cannot open file with name %s. In File %s at line %d.", fileName.c_str(), __FILE__, __LINE__ ).c_str() );
  int iter = 0;
  std::vector<TripletType> tripletList;
  while ( !File.eof() ) {
      int row; File >> row;
      int col; File >> col;
      double val; File >> val;
      tripletList.push_back( TripletType(row, col, val ) );
  }
  File.close();
  mat.setFromTriplets( tripletList.begin(), tripletList.end() );
}


void consoleOutputStartProgramm( const std::string message) {
    cout << endl << endl << endl  
    << "================================================================================" << endl 
    << " Start Programm: " << message << endl 
    << "================================================================================" << endl << endl << endl;
}

void consoleOutput( const std::string message) {
    cout << endl << endl << endl  
    << "================================================================================" << endl 
    << message << endl 
    << "================================================================================" << endl << endl << endl;
}   
    
    
void consoleOutputFinishProgramm( const std::string message) {
    cout << endl << endl << endl  
    << "================================================================================" << endl 
    << " Finished Programm " << message << endl 
    << "================================================================================" << endl << endl << endl;
}
    
    
    
    

    
    
class BoostParser {

    boost::property_tree::ptree _pt;
    
public:
    
    BoostParser ( ) {}
    
    BoostParser ( const std::string & ParFilename ){ boost::property_tree::ini_parser::read_ini ( ParFilename, _pt );}
    
    BoostParser ( const std::string & ParFilename, const int Counter, const string additionalName = "", bool dump = true ) : BoostParser( ParFilename ) {
        this->addCounterToSaveDirectory( Counter, additionalName );
        if( dump ) this->dumpToFile( "ParameterParser.ini" );
    }
    
    BoostParser ( const std::string & ParFilename, const string CounterFileName, const string additionalName = "", bool dump = true ) : BoostParser( ParFilename ) {
        this->addCounterToSaveDirectory( CounterFileName, additionalName );
        if( dump ) this->dumpToFile( "ParameterParser.ini" );
    }
    
    BoostParser& operator = ( const BoostParser& parser ) {
        // Beware of self-assignment
        if ( this == &parser ) return *this;
        _pt = parser._pt;
        return *this;
    }
    
    void merge( const std::string & ParFilename ){
        boost::property_tree::ptree ptUpdates;
        boost::property_tree::ini_parser::read_ini ( ParFilename, ptUpdates );
        BOOST_FOREACH( auto& update, ptUpdates ){
            _pt.put_child( update.first, update.second );
        }
    }    
    
    template<typename VariableType>
    VariableType get ( const std::string variableName ) const { 
        VariableType val;
        try {
            val = _pt.get<VariableType> ( variableName );
        } catch (...){
            throw std::invalid_argument ( aol::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
        }
        return val;
    }
  
  template<typename VariableType, typename VecType>
  void getFreeSizeVector(const std::string& variableName, VecType &vec) const{
      try {
        std::string s = _pt.template get<std::string>(variableName);   
        std::stringstream ss(s);
        std::string item;
        while(std::getline(ss, item, ',')) {
            vec.push_back(boost::lexical_cast<VariableType>(item));
        }
      } catch (...){
            throw std::invalid_argument ( aol::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
      }
  }
  
  template<typename VariableType, typename VecType>
  void getFixSizeVector(const std::string& variableName, VecType &vec) const{
      try {
        std::string s = _pt.template get<std::string>(variableName);   
        std::stringstream ss(s);
        std::string item;
        int counter = 0;
        while(std::getline(ss, item, ',')) {
            vec[counter] = boost::lexical_cast<VariableType>(item);
            counter++;
        }
      } catch (...){
            throw std::invalid_argument ( aol::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
      }
  }
  
  template<typename VariableType, typename MatType>
  void getFixSizeMatrix(const std::string& variableName, MatType &mat) const{
      try {
        const int numRows = mat.rows();
        const int numCols = mat.cols();
        std::string s = _pt.template get<std::string>(variableName);   
        std::stringstream ss(s);
        std::string item;
        int counter = 0;
        int row=-1, col;
        while(std::getline(ss, item, ',')) {
            if( (counter % numCols) == 0 ){
              row++;
              col = 0;   
            }else{
              col++;   
            }
            mat(row,col) = boost::lexical_cast<VariableType>(item);
            counter++;
        }
      } catch (...){
            throw std::invalid_argument ( aol::strprintf ( "error: cannot find %s in File %s at line %d.", variableName.c_str(),  __FILE__, __LINE__ ).c_str() );
      }
  }
    
    
    
  template<typename VariableType>
  void set ( const std::string variableName, const VariableType & variable ) { _pt.put( variableName, variable );}
    
  template<typename VecType>
  void setFixSizeVector(const std::string& variableName, const VecType &vec, const bool useInteger ) {
      std::string s = "";
      if( useInteger ){
        for( int i=0; i<vec.size(); ++i ){
            if( i < vec.size() - 1 ) s += aol::strprintf ( "%d,", vec[i] ).c_str();
            else s += aol::strprintf ( "%d", vec[i] ).c_str();
        }
      }else{
        for( int i=0; i<vec.size(); ++i ){
            if( i < vec.size() - 1 ) s += aol::strprintf ( "%f,", vec[i] ).c_str();
            else s += aol::strprintf ( "%f", vec[i] ).c_str();
        } 
      }
      _pt.put( variableName, s );
      cout << "setFixSizeVector: " << s.c_str() << endl;
  }
    
  void setDirectoryName( const string additionalName ){
        std::string newFileName = _pt.get < std::string > ( "saving.saveDirectory" );
        newFileName += additionalName;
        boost::filesystem::create_directory ( newFileName );
        _pt.put ( "saving.saveDirectory", newFileName );
  }
  
  string createSubDirectory ( const string additionalName ) {
        std::string newFileName = _pt.get < std::string > ( "saving.saveDirectory" );
        newFileName += "/" + additionalName;
        boost::filesystem::create_directory ( newFileName );
        return newFileName;
  }
    
  void addCounterToSaveDirectory ( const int Counter, const string additionalName = "" ) {
        std::string newFileName = _pt.get < std::string > ( "saving.saveDirectory" );
        newFileName += additionalName;
        newFileName += "-";
        newFileName += std::to_string(Counter);
        boost::filesystem::create_directory ( newFileName );
        _pt.put ( "saving.saveDirectory", newFileName );
  }
    
  int addCounterToSaveDirectory ( const std::string CounterFileName, const string additionalName = "" ) {
        if ( !boost::filesystem::exists( CounterFileName ) ) {
            std::ofstream out ( CounterFileName.c_str() );
            out << 0 << std::endl;
            out.close ( );
        }
        std::fstream counterFile;
        counterFile.open ( CounterFileName.c_str () );
        int counter = 0;
        if ( counterFile.is_open () ) {
            std::string temp;
            std::getline ( counterFile, temp );
            counter = atoi ( temp.c_str () );
            counterFile.seekg ( std::ios::beg );
            counterFile << ++counter;
        }
        else throw std::invalid_argument ( "cannot open counter file for writing in File" +  std::string(__FILE__) + "at line " + std::to_string( __LINE__ ) );
        counterFile.close ();
        
        this->addCounterToSaveDirectory( counter, additionalName );

        return counter;
  }
  
  void dumpToFile ( const string& FileName ) const {
      std::string saveName = _pt.get < std::string > ( "saving.saveDirectory" );
      saveName += "/" + FileName;
      boost::property_tree::write_ini ( saveName, _pt );
  }
    
};
    
    


//! \brief Output stream buffer that sends its content to multiple other stream buffers.
class MultiStreambuf : public streambuf {
public:
  MultiStreambuf() { setp(0, 0); }

  size_t addStreambuf(streambuf * buffer) {
    _streambufPtrList.push_back(buffer);
    return _streambufPtrList.size();
  }
  size_t removeStreambuf(streambuf * buffer) {
    _streambufPtrList.remove(buffer);
    return _streambufPtrList.size();
  }

protected:
  typedef std::char_traits<char> traits_type;
  typedef traits_type::int_type  int_type;

  typedef list<streambuf*> StreambufList;
  StreambufList _streambufPtrList;

  int_type overflow( int_type c) {
    if (!traits_type::eq_int_type(c, traits_type::eof())) {
        StreambufList::iterator iter = _streambufPtrList.begin();
        for (; iter != _streambufPtrList.end(); ++iter)
        (*iter)->sputc(c);
    }
    return traits_type::not_eof(c);
  }
  //! synchronizes all buffers. Returns for failure if at least one participation buffer has failed and returned -1. 
  int_type sync() {
    int ret = 0;
    StreambufList::iterator iter = _streambufPtrList.begin();
    for (; iter != _streambufPtrList.end(); ++iter)
        // if ret has already set to value "-1" (failed), do not overwrite this, but keep.
        if ( ret == -1 )
            (*iter)->pubsync();
        // otherwise, give *iter a chance to indicate failure.
        else
            ret = (*iter)->pubsync();
        return ret;
  }
};

//! \brief Print cout and cerr to a file.
class AdditionalOutputToFile {
public:
  AdditionalOutputToFile ( string Filename ) : 
  _filename ( Filename ), _filestream ( Filename.c_str() ), _previousCoutStreambuf ( cout.rdbuf() ), _previousCerrStreambuf ( cerr.rdbuf() ), _previousClogStreambuf ( clog.rdbuf() ) {
    _multiCoutStreambuf.addStreambuf ( _previousCoutStreambuf ); _multiCerrStreambuf.addStreambuf ( _previousCerrStreambuf ); _multiClogStreambuf.addStreambuf ( _previousClogStreambuf );
    _multiCoutStreambuf.addStreambuf ( _filestream.rdbuf() ); _multiCerrStreambuf.addStreambuf ( _filestream.rdbuf() ); _multiClogStreambuf.addStreambuf ( _filestream.rdbuf() );
    cout.rdbuf ( &_multiCoutStreambuf ); cerr.rdbuf ( &_multiCerrStreambuf ); clog.rdbuf ( &_multiClogStreambuf );
  }
  
  ~AdditionalOutputToFile () {
    cout.flush(); cerr.flush(); clog.flush();
    cout.rdbuf ( _previousCoutStreambuf ); cerr.rdbuf ( _previousCerrStreambuf ); clog.rdbuf ( _previousCerrStreambuf );
    clog << "All output has been written to file " << _filename << endl;
  }

protected:
  string      _filename;
  ofstream    _filestream;
  streambuf * _previousCoutStreambuf; streambuf * _previousCerrStreambuf; streambuf * _previousClogStreambuf;
  MultiStreambuf _multiCoutStreambuf; MultiStreambuf _multiCerrStreambuf; MultiStreambuf _multiClogStreambuf;
};



} // namespace aol


    
#endif // __GENERALLOADANDSAVE_H
