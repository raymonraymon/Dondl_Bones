#ifndef __IPOPTINCLUDES_H
#define __IPOPTINCLUDES_H

#ifdef USE_IPOPT

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <cstddef>
#include <general.h>
#include <energyDefines.h>


#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#undef HAVE_CSTDDEF

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpOrigIpoptNLP.hpp"



string getIpoptStatus ( Ipopt::ApplicationReturnStatus ipoptStatus  ) {
    
    string status = "";
    
    switch( ipoptStatus ){
        case 0 :
            status = "Solve Succeeded";
            break;
            
        case 1 :
            status = "Solved To Acceptable Level";
            break;
            
        case 2 :
            status = "Infeasible Problem Detected";
            break;
            
        case 3 :
            status = "Search Direction Becomes Too Small";
            break;
            
        case 4 :
            status = "Diverging Iterates";
            break;
        
        case 5 :
            status = "User Requested Stop";
            break;
            
        case 6 :
            status = "Feasible Point Found";
            break;
            
        case -1 :
            status = "Maximum Iterations Exceeded";
            break;
            
        case -2 :
            status = "Restoration Failed";
            break;
            
        case -3 :
            status = "Error In Step Computation";
            break;
            
        case -10 :
            status = "Not Enough Degrees Of Freedom";
            break;
            
        case -11 :
            status = "Invalid Problem Definition";
            break;
            
        case -12 :
            status = "Invalid Option";
            break;
            
        case -13 :
            status = "Invalid Number Detected";
            break;
            
        case -100 :
            status = "Unrecoverable Exception";
            break;
            
        case -101 :
            status = "NonIpopt Exception Thrown";
            break;
            
        case -102 :
            status = "Insufficient Memory";
            break;
            
        case -199 :
            status = "Internal Error";
            break;
    }
    
    return status;

}


void outputIpoptStatus ( Ipopt::ApplicationReturnStatus ipoptStatus, const bool outputOnlyIfFailed = false ) {
    string status = getIpoptStatus( ipoptStatus );
    // status: 0 - Solve_Succeeded, 1 - Solved_To_Acceptable_Level
    if( !outputOnlyIfFailed && ( ipoptStatus == 0 || ipoptStatus == 1 ) ) cout << "Ipopt finished with status: " << status.c_str() << endl;
}



#endif //USE_IPOPT

#endif //__IPOPTINCLUDES_H
