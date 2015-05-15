// main.cpp -- The main() function for Generalized Logical Network Modeling
//
// Eric Lance, Chris Lewis, Curtis Luce, Joe Song, Yang Zhang
//
// Created: 2006
// Modified: September 2008. Joe Song. Added minimum number of parents
// Modified: September 7, 2008. Joe Song.
// Modified: September 27, 2008. Joe Song. Made major changes by introducing a class to
//   save all the parameters needed for estimation
// add L for number of permutations under comparison mode
// Modified: July 15, 2011. Joe Song. Revised the help screen.
// Last modified: October 1, 2011. Joe Song. Cleaned up unused variables and extracted some functions
//   out of the main

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include <ctime>
#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm> //Added by Haizhou Wang, Oct 7, 2012. In order to use transform()
using std::transform;
using std::cout;
using std::cerr;

#if defined(__APPLE__) || defined(__MACH__)
#include <sys/time.h>
#include <getopt.h>//Hua modified Nov 10 2014
#include <unistd.h>
#elif defined __linux__
#include <sys/time.h>
#include <unistd.h>
#elif defined _WIN32 || defined _WIN64
#else
#error Unknown Operating System Platform
#endif

#include "GLN_copyright.h"
#include "GLNConfig.h"
#include "GLN.h"
#include "GLNGlobals.h"
#include "EnvGLNRec.h"
#include "ComparativeChisq.h"
#include "InteractionEvaluationController.h"

/*extern "C" {
 char * optarg;
 int getopt (int argc, char *const argv[], const char *opts);
 }*/

extern char * optarg;
int getopt (int argc, char *const argv[], const char *opts);

// -------------------------------------------------------------
// Below is Legacy code used in past publications
vector<int> ScanExternalNodes(const string & externalNodesFile);
vector<string> ScanListFile(const char * listfile);
// Above is legacy code.
// -------------------------------------------------------------

void GLNExit(int status, const char * message)
{
   if( message ) {
      cerr << message << endl;
   }
   
#ifdef ENABLE_MPI
   MPI_Finalize();
#endif
   
   exit(status);
}

int main(int argc, char * argv[])
{
#ifdef ENABLE_MPI
   MPI_Init(&argc, &argv);
#endif
   
   GLNConfig gc;
   
   string mode; // simulation, estimation, generation, and more
   
   // user provided seed for srand
   unsigned int seed;
   
#if defined(__APPLE__) || defined(__MACH__)
   timeval tim;
   gettimeofday(&tim, NULL);
   seed = (unsigned int)tim.tv_sec*1000000+tim.tv_usec;
   // srand(tim.tv_sec*1000000+tim.tv_usec);
#elif defined __linux__
   timeval tim;
   gettimeofday(&tim, NULL);
   seed = tim.tv_sec*1000000+tim.tv_usec;
   // srand(tim.tv_sec*1000000+tim.tv_usec);
#elif defined _WIN32 || defined _WIN64
   seed = (unsigned)time( NULL );
   // srand( (unsigned)time( NULL ) );
#else
#error Unknown Operating System Platform
#endif
   InteractionEvaluationController& evalController = InteractionEvaluationController::getController();
   
   
   // -------------------------------------------------------------
   // Parameters used in more than one modes
   string inputfile;               // definition dependent on the mode
   
   string inputfile1;              // definition dependent on the mode
   string inputfile2;              // definition dependent on the mode
   
   string outputfile;              // definition dependent on the mode
   
   // -------------------------------------------------------------
   // Parameters exclusively for "estimation" mode:
   string filePermutationTable;
   
   // -------------------------------------------------------------
   // Parameters exclusively for "estimation" or "comparison" mode:
   vector<string> input_files;		// multiple trajectory files
   string fileCandidate; //added by Yang Zhang 9/12/2009
   
   string methodCmpParentSets;
   bool readTraj = false;			// whether to read in file names of trajectories
   
   //Added by Haizhou Wang, Sep 5, 2012
   /*
    *Pure reconstruction is turned off by default
    *To turn on pure reconstruction, user needs to specify "-F" option in the command line.
    */
   bool fastRecon = false;
   //End of Haizhou's Code
   
   string localConstraintFileName;
   
   
   // -------------------------------------------------------------
   // Parameters exclusively for "generation" mode:
   int networkSize=0, maxBase=0, minBase=0;	// Related to the generateRandomNetwork function
   
   // -------------------------------------------------------------
   // Parameters exclusively for "simulation" mode:
   int simul_number=1;                 // Number of simulations you want to run
   double levelNoise = 0;
   int init_traj = -1;                 // Number of initial trajectories to produce
   string init_trajectory;             // Name of initial trajectory file
   int numboftrj = 1;                  // Number of trajectories to simulate 08/11/16 by yangzhang
   string houseNoiseModel; //added by Yang Zhang 7/22/2012 add command line parameter C for new house noise model
   
   // -------------------------------------------------------------
   // Parameters exclusively for "evaluation" mode:
   string trajectoryfileforcomparison;
   int output_type=1;	//Standard output or file output? 1 is for stndrd, 2 is for file
   string differencefile;
   string summaryfile;
   
   // -------------------------------------------------------------
   // Parameters for comparative pathway analysis:
   string filePathway; //added by yangzhang 1/26/2009 to only enumerate parents within a given pathway
   
   int permuteStrategy = 0;    // default 0: permute On Time(on row across condition),
   // 1: permute gene-wide(on column across condition),
   // 2: is permute on condition
   
   
   // -------------L-E-G-A-C-Y-------------------------------------
   // Below is Legacy parameters used in past publications
   string externalNodeFile;		// Joe Song 08/04/06
   vector<int> externalNodeIndices(0); 	// Joe Song 08/04/06
   
   string fileExtNodeNames;  		// Joe Song 08/05/06  file of external node names
   string fileIntNodeNames;  		// Joe Song 08/05/06  file of internal node names
   
   vector<string> extNodeNames(0);
   vector<string> intNodeNames(0);
   // Above is legacy code.
   // -------------L-E-G-A-C-Y-------------------------------------
   
   
   /******************************************************
    ****** Takes in arguments via a command line **********
    *******************************************************/
   
   int option;	// value returned from getopt()
   
   while( (option=getopt(argc, argv, "1:2:A:a:B:b:C:c:D:d:E:e:Ff:G:g:H:hI:i:J:j:K:k:L:l:M:m:N:n:O:o:P:p:Q:R:r:S:s:T:t:U:v:Vw:X:x:Y:Zz:?")) != -1 ) {
      
      switch(option) {
            
            // ----------------------------------------------------------------------
            // The following is obsolete command arguments used for publication of
            //   Song et al., (2009) EURASIP J Bioinformatics and Systems Biology
         case 'e':
            externalNodeFile = optarg;  // Joe Song 08/04/06
            externalNodeIndices = ScanExternalNodes(externalNodeFile);
            break;								// Joe Song 08/04/06
            
         case 'E':
            extNodeNames = ScanListFile(optarg);
            break;
            
         case 'I':
            intNodeNames = ScanListFile(optarg);
            break;
            
            // The above options are legacy options not used any longer
            // ----------------------------------------------------------------------
            
            
            
            
            // ----------------------------------------------------------------------
         case '1':   // First input file
            inputfile1 = optarg;
            break;
         case '2':	// Second input file
            inputfile2 = optarg;
            break;
            
         case 'a':
            gc.m_adjustedpValFile = optarg;
            break;
            
         case 'A':	// Alpha for the pvalue
            gc.m_alpha = atof(optarg);
            break;
            
         case 'B':	// Maximum base of a node
            maxBase = atoi(optarg);
            break;
            
         case 'b':	// Minimum base of a node
            minBase = atoi(optarg);	//base must be greater than or equal to 2
            if(minBase < 2){
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-b' parameter (min base) must be > 1!"
                       );
            }//end if
            break;
            
         case 'C':
            houseNoiseModel = optarg;
            break;
            
         case 'c':
            trajectoryfileforcomparison = trajectoryfileforcomparison + optarg;
            break;
            
         case 'D':	// Export network to Dot format
            // networkDotFormat = true;
            gc.m_fileDOT = optarg;
            break;
            
            //Added by Haizhou Wang, Jan 23, 2013
         case 'd':  //File name for local cosntraints, to be used in GLN reconstruction mode
            localConstraintFileName = optarg;
            break;
            
            //Added by Haizhou Wang, Sep 5, 2012
         case 'F': //Stands for F(ast) pure reconstruction
            fastRecon = true;
            break;
            //End of Haizhou's Code
            
         case 'f':  // noise level
            levelNoise = atof(optarg);
            break;
            
         case 'g':	//Min number of parents.  Joe Song.  September 5, 2008
            
            gc.m_min_parents_allowed = atoi(optarg);
            //changed by Yang Zhang 4.18.2012, min # of parents can be 0
            if(gc.m_min_parents_allowed < 0){
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-g' parameter (min # parents) must be >= 0!"
                       );
            }//end if
            
            break;
            
         case 'G':
            gc.m_fileIntxCmp = optarg;
            break;
            
         case 'H': //added by Yang Zhang 9/12/2009
            fileCandidate = fileCandidate + optarg;
            break;
            
         case 'i':
            // file name of the initial trajectories, mandatory
            // when the offset is "greater than" -1, i.e., -2, -3, ...
            init_trajectory = optarg;
            break;
            
         case 'h':
         case '?':
            manpage(argv[0]);
            GLNExit(EXIT_SUCCESS);
            
         case 'j':
            gc.m_fileTopology = optarg;
            break;
            
         case 'J':
            gc.m_min_Markov_order = atoi(optarg);
            
            if(gc.m_min_Markov_order > 0){
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-J' parameter (min offset) must be <=0!"
                       );
            }
            break;
            
         case 'k': // added by Yang Zhang 2/19/2010. 1: same parents only (the old one); 2: allow different parents
            switch (atoi(optarg)) {
               case 1:
                  gc.m_acceptableTopologies = SAME;
                  break;
               case 2:
                  gc.m_acceptableTopologies = DIFFERENT;
                  break;
                  
               default:
                  break;
            }
            break;
            
         case 'K'://holds the parents time offset, if not specified it will be zero
            gc.m_max_Markov_order = atoi(optarg);
            
            if(gc.m_max_Markov_order > 0) {
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-K' parameter (max offset) must be <= 0!"
                       );
            }
            break;
            
         case 'L':
            gc.m_numPermutations = atoi(optarg);
            if(gc.m_numPermutations < 0){
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-L' parameter (number of permutations) must be >= 0!"
                       );
            }
            break;
            
         case 'l': // permutation strategy, added by Yang Zhang 9.4.2010
            permuteStrategy = atoi(optarg);
            if(permuteStrategy < 0 || permuteStrategy > 2){
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-l' parameter (Permutation strategy) must be between 0 and 2!"
                       );
            }
            break;
            
         case 'M':	// Modeling mode
            mode = optarg;
            
            //To deal with lowercase or uppercase mode name (Haizhou Wang, Oct 7, 2012)
            transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
            
            if(mode != "generate"
               && mode != "simulation"
               && mode != "evaluate"
               && mode != "estimation"
               && mode != "comparison"
               && mode != "randomstates"
               && mode != "forPoster"
               && mode != "alcohol"
               && mode != "test"
               && mode != "compareglns"
               ) {
               cerr << "ERROR: '-M' parameter (mode) must be " << endl
               << "  generate, simulation, estimation, randomstates, test,"
               << "  comparison, CompareGLNs, or compareGLNs!"
               << endl;
               GLNExit(EXIT_FAILURE);
            }
            
            break;
            
         case 'm':
            differencefile = optarg;
            break;
            
         case 'n':
            summaryfile = optarg;
            break;
            
         case 'N':
            gc.m_startAtNode = atoi(optarg);
            if(gc.m_startAtNode < 1) {
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-N' parameter (start node) must be 1 or bigger!"
                       );
            }
            break;
            
         case 'O':	// Output file
            outputfile = optarg;
            break;
            
         case 'o':	//Output type
            output_type = atoi(optarg);
            if(output_type != 1 && output_type != 2) {
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-o' parameter (output type) must be 1 or 2!"
                       );
            }
            break;
            
         case 'p':	// Max number of parents
            
            gc.m_max_parents_allowed = atoi(optarg);
            
            if(gc.m_max_parents_allowed <= 0){
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-p' parameter (max # parents) must be > 0!"
                       );
            }//end if
            
            else if (gc.m_max_parents_allowed > 10){
               cout << "Warning: Finding 10 parents may take a long time!" << endl;
               
            }//end if
            break;
            
         case 'P':	//Pvalue mode
            
            gc.m_pValueMode= atoi(optarg);
            
            if( gc.m_pValueMode < 1 || (gc.m_pValueMode > 10 && gc.m_pValueMode < 101 && gc.m_pValueMode!=15) || gc.m_pValueMode > 107) {
               // Tyler: Changed to include 6 (G-test) 9/8/2011
               // Tyler: Changed again to include 101-107(conditional modes) 8/1/2012
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-P' parameter (p-value mode) not supported!"
                       );
            }
            evalController.setPvalMode(gc.m_pValueMode);
            break;
            
         case 'Q':
            gc.m_dataType = optarg;
            break;
            
         case 'r':
            seed = atoi(optarg);
            break;
            
         case 'R':
            inputfile = inputfile + optarg;
            break;
            
         case 'S':	//Simulation amount must be greater than 0
            simul_number = atoi(optarg);
            
            if(simul_number < 0){
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-S' parameter (# simulation steps) must be >= 0!"
                       );
            }//end if
            break;
            
         case 's':	//Size of the network
            networkSize = atoi(optarg);
            
            if(networkSize <= 1){
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-s' parameter (network size) must be > 1!"
                       );
            }//end if
            break;
            
         case 'T':	//read in a file that contains a list of trajectory files
            readTraj = true;
            inputfile = inputfile + optarg;
            break;
            
         case 't'://case for the number of initial trajectories the user wants
            //if not provided no initial trajectories will be provided
            init_traj = atoi(optarg);
            if(init_traj < 1){
               
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-t' parameter (number of initial trajectories) must be 1 or more!");
            }
            break;
            
         case 'U':
            filePermutationTable = optarg;
            break;
            
         case 'v':
            gc.m_adjustedDOTFile = optarg;
            break;
            
         case 'V':
            
            cout << "__________________________________________________________________" << endl;
            cout << GLNVersion << endl;
            cout << Copyright << endl;
            cout << Affiliation << endl;
            
            GLNExit(EXIT_SUCCESS);
            
            break;
            
         case 'w':
            gc.setMarginalMethod(optarg);
            break;
            
         case 'X': // Specify a directory (instead of a file) to store the results.
            gc.m_recordResultFile = optarg;
            break;
            
         case 'x': // number of trajectories for simulation
            numboftrj = atoi(optarg);
            if(numboftrj < 0){
               GLNExit(EXIT_FAILURE,
                       "ERROR: '-x' parameter (# trajectories) must be >= 0!"
                       );
            }//end if
            
            break;
            
         case 'Y':
            methodCmpParentSets = optarg;
            break;
            
         case 'Z':
            gc.m_allowSelfCycle = true;
            break;
            
         case 'z':
            filePathway = filePathway + optarg;
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
   
   evalController.setGLNConfig(gc);
   
   srand(seed);
   
   //*****************************************************************************
   
   if (mode == "estimation" || mode == "comparison") {
      
      if( readTraj ) {
         
         // Read the file which contains the name of other trajectory files
         input_files = ScanListFile(inputfile.c_str());
         
      } else {
         
         if(mode == "estimation") {
            
            input_files.push_back(inputfile);
            
         } else if(mode == "comparison") {
            
            input_files.push_back(inputfile1);
            input_files.push_back(inputfile2);
            
         }
      }
      
      // Input trajectory files
      
      gc.m_listTrajFiles = input_files;
      gc.m_fileGLN = outputfile;
      
      gc.m_permStatTable.scan(filePermutationTable.c_str());
      
      //------------------------------------------------------------
      //added by yangzhang 10.27.2009 to get genenames
      TrajectoryCollection trajCol;
      trajCol.scan(input_files[0]);
      
      //------------------------------------------------------------
      if(! filePathway.empty() ) { // Pathway analysis initialization
         if(filePathway.find(",")==string::npos)
         {
            gc.m_pathways.resize(1);
            gc.m_pathways[0].setNodeNames(trajCol.getIntNodeNames());
            gc.m_pathways[0].setNumConditions((int)input_files.size());
            gc.m_pathways[0].scan(filePathway);
            gc.m_pathways[0].setFile(filePathway);
         }else //individual pathways for each condition
         {
            while ( filePathway.length() > 0)
            {
               int m = filePathway.find(',');
               string tempFile;
               if(m>0)
               {
                  tempFile = filePathway.substr(0,m);
                  filePathway = filePathway.substr(m+1);
               }else
               {
                  cerr << "Incorrect parameter format, must end with comma!" << endl;
                  exit(EXIT_FAILURE);
               }
               Topology topo;
               topo.setNodeNames(trajCol.getIntNodeNames());
               topo.setNumConditions((int)input_files.size());
               topo.scan(tempFile);
               topo.setFile(tempFile);
               gc.m_pathways.push_back(topo);
               
               
            }
            
         }
      }
      
      //------------------------------------------------------------
      if(! fileCandidate.empty() ) { // Parent candidate initialization
         //added by YangZhang 9/12/2009
         if(fileCandidate.find(",")==string::npos)
         {
            gc.m_candidateTopologys.resize(1);
            gc.m_candidateTopologys[0].setNodeNames( trajCol.getIntNodeNames() );
            gc.m_candidateTopologys[0].setNumConditions( (int)input_files.size() );
            gc.m_candidateTopologys[0].scan(fileCandidate.c_str());
            gc.m_candidateTopologys[0].setFile(fileCandidate);
         }else //individual candidate topology for each condition
         {
            while ( fileCandidate.length() > 0)
            {
               int m = fileCandidate.find(',');
               string tempFile;
               if(m>0)
               {
                  tempFile = fileCandidate.substr(0,m);
                  fileCandidate = fileCandidate.substr(m+1);
               }else
               {
                  cerr << "Incorrect parameter format, must end with comma!" << endl;
                  exit(EXIT_FAILURE);
               }
               Topology topo;
               topo.setNodeNames(trajCol.getIntNodeNames());
               topo.setNumConditions((int)input_files.size());
               topo.scan(tempFile);
               topo.setFile(tempFile);
               gc.m_candidateTopologys.push_back(topo);
               
               
            }
            
         }
      }
      
      //Added by Haizhou Wang, Jan 24, 2012
      if( !localConstraintFileName.empty() )
      {
         gc.m_mustIncludeParents.setNodeNames( trajCol.getIntNodeNames() );
         gc.m_mustIncludeParents.setNumConditions( (int)input_files.size() );
         
         gc.m_mustIncludeParents.scan( localConstraintFileName );
         gc.m_mustIncludeParents.setFile(localConstraintFileName);
      }
      
      //added by YangZhang 9/15/2010
      gc.m_BP = trajCol.getBeParent();
      
      gc.m_HP = trajCol.getHaveParent();
      
      gc.setMethodCmpParSets(methodCmpParentSets, mode);
      
      if ( mode == "estimation" ) {
         //Added by Haizhou Wang, July 24, 2012
         EnvGLNRec envman;
         envman.initialize(gc);
         
         //Added by Haizhou Wang, June 18, 2012
         GeneralizedLogicalNetwork result;
         
         //Modified by Haizhou Wang, Sep 5, 2012
         /*
          *By default, ReconstructGLN() will be called on the input trajectory (collection) object
          *If '-F' is specified, program will try to invoke ReconstructPureGLN() first if the input is pure.
          *Otherwise, it will fall back to call ReconstructGLN().
          */
         if( fastRecon ) // '-F' option is specified by user
         {
            if( envman.getTraCol().CheckTraDeterministic() )
               result = ReconstructPureGLN(gc);
            else
            {
               cerr<<"Input trajectory (collection) file is not pure. Can not do pure reconstruction. Going to call ReconstructGLN()\n"<<endl;
               result = ReconstructGLN(gc);
            }
         }
         else
         {
            result = ReconstructGLN(gc);
         }
         result.save(outputfile);
         //End of Haizhou's Code
         
      } else if (mode == "comparison") {
         
         if(gc.m_numPermutations == 0) {
            
            if (gc.m_methodCmpxParCmp == BY_EACH_COND
                && gc.m_acceptableTopologies == DIFFERENT) { // 2) {
               
               BuildTopologyThenCompareGLNs(gc);
               
            } else {
               CompareGLNs(gc);
            }
            
         } else if(gc.m_numPermutations > 0){
            
            //ComparePathways(gc, permuteStrategy, input_files);
            CompareGLNsResample(gc);
            
         }
      }
   }
   
   else if(mode == "test") {
      
      if (! inputfile1.empty() && ! inputfile2.empty()) {
         input_files.push_back(inputfile1);
         input_files.push_back(inputfile2);
         applyHeteroChisqTest(input_files, gc.m_pValueMode, gc.m_methodMarginal);
      }else if(! inputfile1.empty() && inputfile2.empty()){
         applyFunctionalChisqTest(inputfile1, gc.m_pValueMode);//Hua Apr 14 2014, apply test to one table
      }else {
         GLNTest();
      }
      
   } else if(mode == "generate") {
      
      cout << "Generalized logical network random generation ... " << endl;
      
      GeneralizedLogicalNetwork gln;
      
      gln.generateRandomNetwork(networkSize, minBase, maxBase,
                                gc.m_min_parents_allowed,
                                gc.m_max_parents_allowed,
                                gc.m_min_Markov_order,
                                gc.m_max_Markov_order,
                                init_trajectory, init_traj);
      
      gln.eliminateFictitiousParents();
      
      gln.save(outputfile);
      
   } else if (mode == "simulation") {
      
      cout << "Generalized logical network simulation ... " << endl;
      
      GeneralizedLogicalNetwork gln;
      gln.scan(inputfile);
      
      cout << "The Markovian order is:  " << gln.getMaxMarkovOrder() << endl;
      
      gln.simulate(init_trajectory, outputfile, numboftrj, simul_number, levelNoise, houseNoiseModel);
      
   } else if(mode == "evaluate") {
      
      cout << "Generalized logical network evaluation ... " << endl;
      
      GeneralizedLogicalNetwork gln;
      gln.scan(inputfile);
      cout << "The Markovian order is:  " << gln.getMaxMarkovOrder() << endl;
      gln.evaluate (output_type, trajectoryfileforcomparison, outputfile, 
                    gln.getMaxMarkovOrder(), differencefile, summaryfile);
      
   } else if(mode == "randomstates") {
      
      cout << "Generalized logical network random generation with random initial states ..." << endl;
      
      GeneralizedLogicalNetwork gln;
      
      gln.initialStates(inputfile, init_trajectory, init_traj, 
                        gc.m_min_Markov_order, gc.m_max_Markov_order);
      
   } else if(mode == "compareglns") {
      
      GeneralizedLogicalNetwork gln1, gln2;
      gln1.scan(inputfile1);
      gln2.scan(inputfile2);
      
      gln1.compareAllTruthTables(gln2, gc.m_fileIntxCmp);
      
   } else {
      
      cerr << "ERROR: Invalid mode \"" << mode << "\"!" << endl;
      GLNExit( EXIT_FAILURE );        
   }
   
   cout << "done" << endl;
   
#ifdef ENABLE_MPI
   
   MPI_Finalize();
   
#endif
   
   return 0;
   
}

//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>

// _CrtDumpMemoryLeaks();

