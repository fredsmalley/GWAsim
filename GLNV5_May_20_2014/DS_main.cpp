#include <algorithm>
#include <ctime>
#include <fstream>
using std::ifstream;
using std::ofstream;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <string>
using std::string;
using std::to_string;

#include "DS_LocalConstraint.h"
#include "GLNGlobals.h"
#include "kStop.h"
#include "kStop_TestCase.h"
#include "Topology.h"

#if defined(__APPLE__) || defined(__MACH__)  
#include <sys/time.h>
#include <unistd.h>
#elif defined __linux__
#include <sys/time.h>
#include <unistd.h>
#elif defined _WIN32 || defined _WIN64
#else
#error Unknown Operating System Platform
#endif


void GLNExit(int status, const char * message) 
{
    if( message )
        cerr << message << endl;

    exit(status);
}



extern char * optarg;
int getopt (int argc, char *const argv[], const char *opts);



int main(int argc, char * argv[])
{
    srand( static_cast<unsigned>( time(NULL) ) );

	string localConstraint_inputfile;
	string globalConstraint_inputfile;
	
	string outputfile("");
    string mode;
    string SearchDirection("random"); //Default search order to traverse one node's adjacent list in IDDFS
    int depth = 5; 					  //Specifies the maximum depth IDDFS can go

    /***************************************/
    /*Following parameters are used in test mode*/
    double percentage = 1.0;
    bool incompleteDynamicPattern = false;
    string inputGLNFileName;
    int numOfGlobalConstraints = 5;
    Topology mustIncludeParents;
    /***************************************/
    
    int option;	// value returned from getopt()
	while( (option=getopt(argc, argv, "c:d:D:l:g:IM:n:O:T:t:?")) != -1 ) {

		switch(option) {
            
            case 'c': //Will be used in test mode. This percentage indicates that how much the orginal GLN will be kept. If a node has 7 parents and percentage is 0.65, then it will left with 4 parents in the resulting GLN (7*0.65=4.55, numbers after decimal point are truncated).
                percentage = atof(optarg);
                if( percentage<=0 || percentage > 1 )
                {
                    GLNExit(EXIT_FAILURE, "ERROR: '-c' parameter (percentage) must be greater than 0 and less or equal to 1");
                }
                break;

            case 'd':
                depth = atoi(optarg);
                if( depth < 0 ) 
                {
                    GLNExit(EXIT_FAILURE, "ERROR: '-d' parameter (depth) must be non-zero");
                }
                break;
            
            case 'D':
                SearchDirection = optarg;
                transform(SearchDirection.begin(), SearchDirection.end(), SearchDirection.begin(), ::tolower);
                
                if(SearchDirection != "forward" && SearchDirection != "backward" && SearchDirection != "random") 
                {
                    GLNExit(EXIT_FAILURE, "ERROR: '-D' parameter (Direction) must be forward, backward or random\n");
                }
                break;

			case 'l':   // input file for local constraint
				localConstraint_inputfile = optarg;
				break;

			case 'g':	// input file for global constraint
				globalConstraint_inputfile = optarg;
				break;

            case 'I':
                incompleteDynamicPattern = true;
                break;

			case 'M':	// Mode: exhaustive, effcgln, kstop
				mode = optarg;			
                transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

				if(
	  	   mode != "kstop" &&
                   mode != "plot" &&
                   mode != "test" &&
                   mode != "edge" ) 
                {
                        GLNExit(EXIT_FAILURE, "ERROR: '-M' parameter (mode) must be 'exhaustive', 'effcgln', 'kstop' or 'test'.\n");
				}
				break;

            case 'n':
                numOfGlobalConstraints = atoi( optarg );
                if( numOfGlobalConstraints < 2 )
                {
                    GLNExit( EXIT_FAILURE, "number of global constriants must be at least 2.");
                }
                break;

			case 'O':	// Output file
				outputfile = optarg;			
				break;

            case 'T':
                inputGLNFileName = optarg; //Input GLN file name for 'test' mode
                break;

            case 't':
                mustIncludeParents.scan( optarg );
                break;

            default:                
				cerr << "ERROR: invalid option -" << (char) option;
				if(optarg) {
					cerr << ' ' << optarg;
				}
				cerr << "!" << endl;
				GLNExit(EXIT_FAILURE);

		} //end switch
	} //end while


	//*****************************************************************************
    if( mode == "edge")
    {
        GeneralizedLogicalNetwork gln;
        gln.scan( inputGLNFileName + ".gln" );

        gln.saveEdgeList( inputGLNFileName + "-EdgeList.txt");
    }else if(mode == "kstop"){
        DS_LocalConstraint lc(localConstraint_inputfile);
        EffCGLN_GlobalConstraint gc(globalConstraint_inputfile);

        
        string output_file_name(outputfile);
        
        lc.saveAsTopologyFile( output_file_name + "LocalConstraintTopology.txt");
        mustIncludeParents.scan( output_file_name + "LocalConstraintTopology.txt");
        

        vector<GeneralizedLogicalNetwork> GLNList;
        vector<string> GLNOutputtedNames;
        GeneralizedLogicalNetwork tmp;

        int t = 0;
        //If there are too many satisfying candidates, only use the first 100.
        while(t<10)
        {
            //kStop::kStop(local_constraint, global_constraint, traverse_order, deepth_limit_for_IDDFS)
            kStop instance( lc, gc, SearchDirection, depth);

            if ( !instance.FindPath() )
            {
                cerr<<"Using search detph "<<depth<<", no solution could be found with the input local and global constraints. Considering increase the search detph if viable.\n"<<endl;
                return 0;
            }
        
            
            if( !instance.getTrajectoryCollection().CheckTraDeterministic() )
            {
                GLNExit( EXIT_FAILURE, "kStop() produced an impure trajectory collection (It shouldn't).\n Program will be terminated.\n");
            }
            

            instance.saveTra(output_file_name + "-" + to_string(t) + ".tra");

            tmp = ReconstructPureGLN( instance.getGLNConfig( output_file_name, to_string(t) ) );

            if( find( GLNList.begin(), GLNList.end(), tmp ) != GLNList.end() )
                break;
            
            GLNList.push_back( tmp  );
            
            tmp.save( output_file_name + to_string(t) + ".gln");
            GLNOutputtedNames.push_back( output_file_name + to_string(t) + ".gln");
            ++t;
        }//end of while()
        assert( GLNList.size() != 0 );


        //<parent,child>:occurrence
        vector< pair< pair<int,int>, unsigned> > edgeRank( kStop_TestCase::RankEdges( GLNList ) );
        
        ofstream o("Edge_List_Result.txt");

        o<<"<Parent, Child>\tCount"<<endl;
        for( const auto & e: edgeRank)
        {
            o<<GLNList.front().getNodeName(e.first.first-1)<<"\t"<<GLNList.front().getNodeName(e.first.second-1)<<"\t"<<(double)e.second/GLNList.size()*100<<"%"<<endl;
        }
        o.close();


        kStop_TestCase::outputGLNsRankedByEdgeMedian( edgeRank, GLNList, GLNOutputtedNames);
    }else if( mode == "plot" )
    {
        if( inputGLNFileName.empty() ){
            GLNExit(EXIT_FAILURE, "Please specify a gln file with option '-T'\n");
        }

        if( outputfile.empty() ){
            GLNExit(EXIT_FAILURE, "Please specify an output file name for the dot file.\n");
        }


        GeneralizedLogicalNetwork gln;
        gln.scan( inputGLNFileName );
        gln.saveEdgesDotFormat( "TransitionDiagram" + outputfile );
        gln.saveDot( "Topology" + outputfile );

    }else if(mode == "test") {
        if( inputGLNFileName.empty() ){
            GLNExit( EXIT_FAILURE, "No input gln file for 'test' mode.\n");
        }


        kStop_TestCase::test(inputGLNFileName, incompleteDynamicPattern, percentage, numOfGlobalConstraints, SearchDirection, depth);
    }else {
		cerr << "ERROR: Invalid mode \"" << mode << "\"!" << endl;
		GLNExit( EXIT_FAILURE );
	}


	return 0;
}

