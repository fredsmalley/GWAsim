// GLN_Test.cpp -- This program will call estimate and also simulate the data
// 
// Curtis Luce, Joe Song
// Created: 2007
// Last modified: September 27, 2008. Joe Song. Name Changed from test.cpp to GLN_Test.cpp
//Oct.29.2009, Yang Zhang, add code to test both reconstruction and comparative GLN
//Jan.20. 2013 Yang Zhang Add test case for fisher exact tests 

#include <ctime>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>

using std::stringstream;
using std::cout;
using std::endl;
using std::cerr;
using std::fstream;
using std::ios;
using std::ofstream;

#include <algorithm>
#include <string>
using std::string;

#include <vector>
using std::vector;

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "GLN.h"
#include "EnvGLNRec.h"
#include "GLNConfig.h"
#include "TrajectoryCollection.h"
#include "Topology.h"

#include "GLNGlobals.h"

void GLNTest()
{
    clock_t begin_clock = clock();

#ifdef ENABLE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
#endif
        
        void test_chisq_upper_bound();
        test_chisq_upper_bound();
        
        void test_ChisqPvalue();
        test_ChisqPvalue();
        
        void test_Multinomial();
        test_Multinomial();
        
        // cout << "BEGIN test_Hetero() ... " << endl;
        void test_Hetero();
        test_Hetero();
        // cout << "END test_Hetero()." << endl << endl;
        
        // cout << "BEGIN test_TwoSampleChisqTest() ..." << endl;
        void test_TwoSampleChisqTest();
        test_TwoSampleChisqTest();
        // cout << "END test_TwoSampleChisqTest()." << endl << endl;
        
        // cout << "BEGIN test_KSampleChisqTest() ... " << endl;
        void test_KSampleChisqTest();
        test_KSampleChisqTest();
        // cout << "END test_KSampleChisqTest()." << endl << endl;
        
        void test_CompareGLNsResample();
        test_CompareGLNsResample();        
		
#ifdef ENABLE_MPI
    }
#endif
    
    void test_gln(void);
    test_gln();
    
    void test_comparative_analysis();
    test_comparative_analysis();
    
#ifdef ENABLE_MPI
    if(rank == 0) {
#endif
        
        void test_GLN_exhaust_simulation();
        test_GLN_exhaust_simulation();
        
        void test_gln_fictitious() ;
        test_gln_fictitious();
        
        //void test_fisher_exact();
        //test_fisher_exact();
        
        void test_EMT_fisher();
        test_EMT_fisher();
        
        //void test_fisher_comparative();
        //test_fisher_comparative();

        // void test_exact_hetero();
        // test_exact_hetero();
        
        void test_exact_hetero_3();
        test_exact_hetero_3();

        void test_exact_comparative(short version);
        //test_exact_comparative(1);
        //test_exact_comparative(2);
        test_exact_comparative(3);
        
        clock_t end = clock();
        double elapsed_secs = ((double) end - begin_clock) / CLOCKS_PER_SEC;
        
        unsigned minutes = floor(elapsed_secs / 60.0);
        unsigned seconds = elapsed_secs - minutes * 60;
        cout << "Test completed in " << minutes << " minutes "
        << seconds << " seconds." << endl;
        
#ifdef ENABLE_MPI
    }
#endif
    
}

void test_GLN_exhaust_simulation()
{
    cout << ">>>> Testing GLN simulation ..." << endl;
    
    GeneralizedLogicalNetwork gln;
    
    gln.generateRandomNetwork(3, 2, 3, 1, 3, -1, -1, "", -1);

    gln.setNodeType(0, 'e');
    
    string fileTCF = "test_sim_exhaust.txt";
    gln.simulate("", fileTCF, 0, 0, 0);

    TrajectoryCollection trjCol;
    trjCol.scan(fileTCF);
    // trjCol.print();

    vector<int> state(gln.size(), 0);

    int pos;
    
    do {
        size_t i;
        
        for (i=0; i<trjCol.getNumTrajectories(); i++) {

            if(trjCol.getTrajectory(i).end() !=
               find(trjCol.getTrajectory(i).begin(), 
                    trjCol.getTrajectory(i).end(),
                    state) ) {
                   break;
            }
            
        }
        
        if(i == trjCol.getNumTrajectories()) {
        
            cerr << "ERROR: ";
            
            for(size_t j=0; j<state.size(); j++) { 
                cerr << state[j] << ',';
            }
            cerr << " was not visited in the exhaustive simulation!" << endl;
            exit(EXIT_FAILURE);
            
        } else {
            pos = EnumMethod::nextAssignmentWithRepeat(state, gln.getBases());
        }
        
    } while (pos >= 0);
    
    cout << "OK." << endl;

    remove(fileTCF.c_str());
    
}

void test_comparative_analysis()
{
    string tcf1 = "test_gln1.trj";
    string tcf2 = "test_gln2.trj";

#ifdef ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
#endif
    cout << ">>>> Testing GLN comparative analysis ..." << endl;
    
    GeneralizedLogicalNetwork gln1, gln2;
    
    gln1.generateRandomNetwork(5,3,3,1,3,-1,-2,"",-1);
    gln2.generateRandomNetwork(5,3,3,1,3,-1,-2,"",-1);

    string file_gln1 = "test_gln1.txt", file_gln2 = "test_gln2.txt";
    
    gln1.save(file_gln1);
    gln2.save(file_gln2);
    
    gln1.scan(file_gln1);
    gln2.scan(file_gln2);
    
    gln1.simulate("", tcf1, 40, 10);
    gln2.simulate("", tcf2, 40, 10);

#ifdef ENABLE_MPI
  } 
#endif

    //the following is to test GLN comparison
    GLNConfig gc;
    
    gc.m_alpha = 0.05;
    gc.m_acceptableTopologies = SAME; // DIFFERENT;
    gc.m_fileDOT = "diff.dot";
    
    gc.m_methodCmpxParCmp = BY_EACH_COND; // BY_TOT; // BY_BEST_HOM; // Joe Song 11/3/2011
        
    gc.m_listTrajFiles = vector<string>(2); 
    gc.m_listTrajFiles[0] = tcf1;
    gc.m_listTrajFiles[1] = tcf2;

    TrajectoryCollection trajCol;
    trajCol.scan(tcf1);
    
    gc.m_BP = trajCol.getBeParent();
    gc.m_HP = trajCol.getHaveParent();
    
    
    CompareGLNs(gc);
    
    if(gc.m_methodCmpxParCmp == BY_EACH_COND) {
        BuildTopologyThenCompareGLNs(gc);
    }
    
#ifdef ENABLE_MPI
    if(rank == 0) {
#endif
        remove(tcf1.c_str());
        remove(tcf2.c_str());
        remove("test_gln1.txt");
        remove("test_gln2.txt");
        remove("-Conserved.dot");
        remove("-Diff.dot");
        remove("-WkZn.dot");
        
#ifdef ENABLE_MPI
    }
#endif
    
}

void test_gln() 
{
    int compareTraj(const string & fileTrj1, const string & fileTrj2);

    GeneralizedLogicalNetwork gln1;

#ifdef ENABLE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
#endif

    cout << ">>>> Testing GLN modeling ..." << endl;
    
    //the following code in comment is to test gln simulation and how well we can reconstruct from them
    //i.e. can we recover the GLN from the simulated trjactory?
    //results in testGLNComparison should be like this
	/*
	C	1	0	2	2
	C	1	0	1,5	1,5
	C	1	0	1,5	1,5
	D	1	0.0833333	2,5	2,5
	C	1	0	2,4	2,4
	*/
    
	string testGLNfile = "testGLN.txt";
    
	ofstream out(testGLNfile.c_str(),ios::out);
	out<<"GENERALIZED_LOGICAL_NETWORK_VER2"<<endl;
	out<<"5\n1\nNode1\n2\n0\ni\n1\n2\n-2\n2\n1\n0\n0\t1\n1\t0\n*"<<endl;
	out<<"2\nNode2\n2\n1\ni\n2\n1\t5\n-1\t-1\n2\t2\n0\n1\n1\n1\n1\t0\n0\t1\n0\t1\n0\t1\n*"<<endl;
	out<<"3\nNode3\n3\n2\ni\n2\n1\t5\n-1\t-1\n2\t2\n0\n0\n1\n2\n1\t0\t0\n1\t0\t0\n0\t1\t0\n0\t0\t1\n*"<<endl;
	out<<"4\nNode4\n3\n2\ni\n2\n2\t5\n-2\t-1\n2\t2\n0\n2\n1\n0\n1\t0\t0\n0\t0\t1\n0\t1\t0\n1\t0\t0\n*"<<endl;
	out<<"5\nNode5\n2\n0\ni\n2\n2\t4\n-1\t-1\n2\t3\n0\n1\n1\n1\n1\n1\n1\t0\n0\t1\n0\t1\n0\t1\n0\t1\n0\t1\n*"<<endl;
	out<<"done"<<endl;
	out.close();
	
    GeneralizedLogicalNetwork gln;

	gln.scan(testGLNfile);

	// Simulation 1
    cout <<  endl << "   .... Simulation 1 ... " << endl;
	gln1.scan("testGLN.txt");
	gln1.simulate("","testGLN.trj",80,10);


    // Reconstruction
    cout << endl << "   .... Reconstruction ... " << endl;

#ifdef ENABLE_MPI
  }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    GLNConfig gc;
    gc.m_pValueMode = 3;
    gc.m_alpha = 1;
    gc.m_min_Markov_order = -1; 
	gc.m_max_Markov_order = -2;
	gc.m_startAtNode = 1;
    
	// Input trajectory files
	vector<string> input_files(1, "testGLN.trj");
	gc.m_listTrajFiles = input_files;

	// Output files
	gc.m_fileGLN = "testGLNest.txt";
	gc.m_allowSelfCycle = false;
	gc.m_max_parents_allowed = 2;
	gc.m_min_parents_allowed = 1; 
	vector<bool> BP;
	vector<bool> HP;
	BP.resize(5);
	HP.resize(5);
	for(int i=0;i<5;i++)
	{
		BP[i] = true;
		HP[i] = true;
	}
	gc.m_HP = HP;
	gc.m_BP = BP;
	ReconstructGLN(gc);

#ifdef ENABLE_MPI
  if(rank == 0) {
#endif

    // Simulation 2
    cout << endl << "   ... Simulation 2 ... " << endl;
    GeneralizedLogicalNetwork gln2;
    gln2.scan("testGLNest.txt");
    gln2.simulate("testGLN.trj","testGLNest.trj",80,10);
	
    // Compare trajectories from simulation 1 and 2
    if(! compareTraj("testGLN.trj", "testGLNest.trj")) {
        
        cout << "Simulated trajectory from reconstructed GLN is the same with the original ... OK." << endl;
        
    } else {

        cout << "Simulated trajectory from reconstructed GLN is different from the original ... ERROR." << endl;

    }

    gln1.compareAllTruthTables(gln2, "testGLNComparison.txt");

    // Fisher exact test
    cout << "   .... Reconstruction using fisher exact test (mode 8)  .... " << endl;
    gln1.simulate("","testGLN-short.trj",10,3);

#ifdef ENABLE_MPI
  }
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    GLNConfig gcfisher = gc;
    gcfisher.m_listTrajFiles = vector<string>(1, "testGLN-short.trj");
    gcfisher.m_fileGLN = "testGLN-fisher.txt"; // Output GLN file
      
    gcfisher.m_pValueMode = 8;
    ReconstructGLN(gcfisher);

    cout << "   ....  using fisher exact test (mode 9)  .... " << endl;
    gcfisher.m_pValueMode = 9;
    ReconstructGLN(gcfisher);

    cout << endl << ">>>> Comparative pathway analysis ... " << endl;

    string testPathwayFile = "testPathway.txt";
    
    string testPathwayFile2 = "testPathway2.txt";
    // remove(testPathwayFile.c_str());

#ifdef ENABLE_MPI
    if(rank == 0) {
#endif
    //topology for condition 1
    ofstream outpathway(testPathwayFile.c_str(), ios::out);
    outpathway<<"pathway"<<endl;
    outpathway<<"5"<<endl;
    outpathway<<"Node1*Node2,"<<endl;
    outpathway<<"Node2*Node1,Node5,Node3,"<<endl;
    outpathway<<"Node3*Node1,Node5,"<<endl;
    outpathway<<"Node4*Node2,Node5,"<<endl;
    outpathway<<"Node5*Node2,Node4,"<<endl;
    outpathway.close();
    //topology for condition 2
    outpathway.open(testPathwayFile2.c_str(), ios::out);
    outpathway<<"pathway"<<endl;
    outpathway<<"3"<<endl;
    outpathway<<"Node1*Node2,Node3,Node4,"<<endl;
    outpathway<<"Node2*Node1,Node5,"<<endl;
    outpathway<<"Node3*Node1,Node4,"<<endl;
    outpathway.close();

#ifdef ENABLE_MPI
  }
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    //the following is to test GLN comparison
    gc.m_pValueMode = 3;
    gc.m_alpha = 0.05;
    gc.m_acceptableTopologies = DIFFERENT; // 2;
    gc.m_methodCmpxParCmp = BY_EACH_COND; // MS 6/2/2013 BY_TOT; // BY_BEST_HOM; Joe Song 11/3/2011
    
    gc.m_fileTopology = "./testTopology.txt";
    
    gc.m_listTrajFiles = vector<string>(2); 
    gc.m_listTrajFiles[0] = "./testGLN.trj";
    gc.m_listTrajFiles[1] = "./testGLNest.trj";
    
    TrajectoryCollection trajColFile;
	
    trajColFile.scan(gc.m_listTrajFiles[0]);
	
    vector<string> geneNames = trajColFile.getIntNodeNames();
    
    gc.m_candidateTopologys.resize(1);
    gc.m_candidateTopologys[0].setNodeNames( geneNames );
    gc.m_candidateTopologys[0].setNumConditions( input_files.size() );
    gc.m_candidateTopologys[0].scan(testPathwayFile.c_str());
    
    cout << "   ....  using one topolgy for both conditions  .... " << endl;
    
    CompareGLNs(gc);
    
    //add test case to use two topologies for two conditions

    gc.m_candidateTopologys.clear();
    gc.m_pathways.resize(2);
    
    gc.m_alpha = 1;
    
    gc.m_pathways[0].setNumConditions(2);
    gc.m_pathways[0].setNodeNames(geneNames);
    gc.m_pathways[0].scan(testPathwayFile.c_str());
    
    gc.m_pathways[1].setNumConditions(2);
    gc.m_pathways[1].setNodeNames(geneNames);
    gc.m_pathways[1].scan(testPathwayFile2.c_str());

    cout << "   ....  using separate topolgy files for each condition  .... " << endl;
    
    CompareGLNs(gc);

  
#ifdef ENABLE_MPI
    if(rank == 0) {
#endif    
    remove("testGLN.txt");
    remove("testGLN.trj");
    remove("testGLNest.txt");
    remove("testGLNest.trj");
    remove("testGLNComparison.txt");
    remove("testGLN-short.trj");
    remove("testGLN-fisher.txt");
    remove(testPathwayFile.c_str());
    remove("./testTopology.txt");
    remove("./testGLN.trj");
    remove("./testGLNest.trj");
    remove(testPathwayFile2.c_str());
#ifdef ENABLE_MPI
    }
#endif
}


void test_gln_fictitious() 
{
	cout << "\n>>>> Testing GLN with fictitious parents ..." << endl;

	string testGLNfile = "testGLN.gln";

	ofstream out(testGLNfile.c_str(),ios::out);
	out<<"GENERALIZED_LOGICAL_NETWORK_VER2"<<endl;
	out<<"4\n";
	out<<"1\nNode1\n2\n0\ni\n3\n2\t3\t4\n-1\t-1\t-1\n2\t3\t2\n0\n0\n1\n0\n0\n1\n0\n0\n1\n0\n0\n1\n1\t0\n1\t0\n0\t1\n1\t0\n1\t0\n0\t1\n1\t0\n1\t0\n0\t1\n1\t0\n1\t0\n0\t1\n*"<<endl;
	out<<"2\nNode2\n2\n1\ne\n0\n*"<<endl;
	out<<"3\nNode3\n3\n0\ne\n0\n*"<<endl;
	out<<"4\nNode4\n2\n1\ne\n0\n*"<<endl;
	out<<"done"<<endl;
	out.close();

	GeneralizedLogicalNetwork gln;
	gln.scan(testGLNfile);

	cout <<  endl << "   .... Call void eliminateFictitiousParents();... " << endl;
	gln.eliminateFictitiousParents();

	cout <<  endl << "   .... Save to file result.gln ... " << endl;
	gln.save("result.gln");
    
    remove("result.gln");
    remove(testGLNfile.c_str());
}
