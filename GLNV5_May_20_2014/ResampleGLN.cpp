//
//  ResampleGLN.cpp
//  gln
//
//  Created by Joe Song on 4/23/13.
//
//

#include <cstdio>
using namespace std;

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "GLNConfig.h"
#include "EnvGLNCmp.h"
#include "PathwayStats.h"
#include "GLNGlobals.h"

const EnvGLNCmp &
EnvGLNCmp::
modify(int min_Markov_order, int max_Markov_order,
       const vector<TrajectoryCollection> & trajCols)
{
    setMarkovOrder(min_Markov_order, max_Markov_order);

#ifdef ENABLE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) { // manager
        
        // pack trajCols
        // send trajCols package to each worker node
        
    } else { // worker
        
        // recieve trajCols package from manager
        // unpack trajCols
    }
#endif
    
    setTrajCols(trajCols);
    
    compareChildWorkingZones(); // due to changed trajectories
    
    return *this;
}

EnvGLNCmp CompareGLNsResample(const GLNConfig & gcObs,
                              const string & askedOption)
{
    // Note:
    // option can take one of "BY_FILE"; "BY_MEM"; "BY_BOTH"
    // BY_FILE is designed work for both single processor and MPI, but
    //   is not tested under MPI.
    // BY_MEM is fast, but not implemented completely for MPI mode.
    // BY_BOTH will run both above modes and report an error if the results
    //   are inconsistent

    string option = askedOption;
    
    if (gcObs.m_methodCmpxParCmp == BY_EACH_COND
        && gcObs.m_acceptableTopologies == DIFFERENT) {
        option = "BY_FILE"; // MS. 5/13/2014. Force to take this option as a
                            // temporary solution.
                            // The "BY_MEM" option cannot handle this condition.
    }
    
    EnvGLNCmp envObs;
    
    envObs.initialize(gcObs);
    
    size_t nBootstraps = gcObs.m_numPermutations;

    ResampleType strategy = BOOTSTRAP;
    
    if (gcObs.m_max_Markov_order == 0 &&
        gcObs.m_min_Markov_order == 0) {
        strategy = BOOTSTRAP_PLUS_ONE;
    }
    
    vector<MultiplePathwayStats> stats(nBootstraps);
    
    cerr << "Bootstrapping ";
    
    for(size_t i=0; i<nBootstraps; i++) {

        cerr << ".";
        
        EnvGLNCmp env;
        MultiplePathwayStats pwStatsByFile, pwStatsByMem;

        vector<TrajectoryCollection> trajCols;
        
#ifdef ENABLE_MPI
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) // manager/master node only
#endif
        // Generate bootstrap trajectory collections
        trajCols = resample(envObs.getTrajCols(), strategy);
        // end of MPI manager processing.
    
        if (option == "BY_FILE" || option == "BY_BOTH") { // Share by disk files

            GLNConfig gc = gcObs;
            
            gc.m_verbose = false; // turn off screen output
            
            for (size_t j=0; j<gc.m_listTrajFiles.size(); j++) {

                // Save trajCols to a list of tempory bootstrap TCFs
                // gc.m_listTrajFiles[j] = to_string(i+1) + "-" + to_string(j+1) + ".TCF";
                gc.m_listTrajFiles[j] =  "bootstrap-c" + to_string(j+1) + ".TCF";
#ifdef ENABLE_MPI
                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                if(rank == 0) // manager/master node only
#endif
                trajCols[j].save(gc.m_listTrajFiles[j].c_str());
                // end of MPI manager processing.

            }
            
            if (strategy == BOOTSTRAP_PLUS_ONE) {
                gc.m_max_Markov_order = -1;
                gc.m_min_Markov_order = -1;
            }
            
            if (gc.m_methodCmpxParCmp == BY_EACH_COND
                && gc.m_acceptableTopologies == DIFFERENT) { // 2) {
                
                env = BuildTopologyThenCompareGLNs(gc);
                
            } else {
                env = CompareGLNs(gc);
            }

            // env = CompareGLNs(gc);
            
            pwStatsByFile = env.getPathwayStats();

#ifdef ENABLE_MPI
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if(rank == 0) // manager/master node only
#endif
                for (size_t j=0; j<gc.m_listTrajFiles.size(); j++) {
                    remove(gc.m_listTrajFiles[j].c_str());
                }
            
        }
            
        if (option == "BY_MEM" || option == "BY_BOTH") { // Share by memory
            
            // envTest.initialize(gcTest); // gc
            env = envObs;

            int min_Markov_order, max_Markov_order;
            if (strategy == BOOTSTRAP_PLUS_ONE) {
                min_Markov_order = -1;
                max_Markov_order = -1;
            } else {
                min_Markov_order = gcObs.m_min_Markov_order;
                max_Markov_order = gcObs.m_max_Markov_order;
            }
            
            env.modify(min_Markov_order, max_Markov_order, trajCols);
            
            CompareGLNs(env);

            env.gatherPathwayStats();
            
            pwStatsByMem = env.getPathwayStats();
        }
        
        if (option == "BY_BOTH") {
            if (pwStatsByFile != pwStatsByMem) {
                cerr << "ERROR: inconsistent pathway stats between BY_FILE and BY_MEM!" << endl;
                exit(EXIT_FAILURE);
            }
        }
        
        // Record bootstrap pathway statistics from envBS
        stats[i] = env.getPathwayStats();
    }
    
    cerr << endl;
    
    if (gcObs.m_methodCmpxParCmp == BY_EACH_COND
        && gcObs.m_acceptableTopologies == DIFFERENT) {
        envObs = BuildTopologyThenCompareGLNs(gcObs);
    } else {
        CompareGLNs(envObs);
    }

    // CompareGLNs(envObs);

    envObs.finalize();
    
    // Approximate null distribution and calculate p-value
    MultiplePathwayStats pwStats = envObs.getPathwayStats();
    
    pwStats.approximate(stats);

    envObs.setPathwayStats(pwStats);
    envObs.savePathwayStats();
    envObs.saveAllInteractionToDot();
    
    return envObs;
}

void ReconstructGLNsResample(const GLNConfig &  gcObs)
{
    // ReconstructGLNs(gcObs);
    
    // Loop:
    
        // Boot strap to produce new TrajCollectoins
    
        // ReconstructGLNS
    
        // Accumulate statistics
    
    // Approximate null distribution
    
    // Calculate p-value
}

void test_CompareGLNsResample()
{
    cout << ">>>> Testing Resample ChiNet on comparative pathway analysis..." << endl;
    
    GLNConfig gcObs;
    gcObs.m_alpha = 1;
    gcObs.m_min_Markov_order = 0;
    gcObs.m_max_Markov_order = 0;
    gcObs.m_methodCmpxParCmp = BY_EACH_COND;
    
    gcObs.m_listTrajFiles = vector<string>(2);
    
    //creating TRJ1.trj
    string TCF1 = "TRJ1.trj";
    
	ofstream out(TCF1.c_str(),ios::out);
	out<<"TRAJECTORY_VER2"<<endl;
	out<<"1\t4\t0\n";
    out<<"2\t2\t2\t2\n";
    out<<"Node1\tNode2\tNode3\tNode4\n";
    out<<"10\n";
    out<<"1\t0\t1\t0\t0\t1\t1\t0\t1\n";
    out<<"0\t1\t1\t1\t1\t1\t0\t1\t0\n";
    out<<"1\t1\t0\t1\t0\t0\t0\t1\t1\n";
    out<<"0\t1\t0\t0\t1\t1\t1\t1\t1\n";
    out<<"1\t1\t1\t0\t0\t0\t1\t0\t0\n";
    out<<"1\t0\t0\t1\t1\t1\t1\t0\t0\n";
    out<<"1\t1\t1\t0\t0\t0\t0\t1\t1\n";
    out<<"0\t0\t0\t1\t1\t1\t0\t0\t0\n";
    out<<"1\t1\t1\t1\t0\t0\t0\t0\t0\n";
    out<<"1\t0\t0\t1\t0\t0\t1\t0\t1\n";
	out.close();
    
    //creating TRJ2.trj
    string TCF2 = "TRJ2.trj";
    
	out.open(TCF2.c_str(),ios::out);
	out<<"TRAJECTORY_VER2"<<endl;
	out<<"1\t4\t0\n";
    out<<"2\t2\t2\t2\n";
    out<<"Node1\tNode2\tNode3\tNode4\n";
    out<<"10\n";
    out<<"0\t1\t0\t1\t1\t1\t0\t1\t1\n";
    out<<"1\t1\t1\t1\t1\t0\t1\t0\t1\n";
    out<<"0\t0\t1\t0\t0\t1\t0\t1\t1\n";
    out<<"1\t0\t0\t1\t1\t1\t1\t1\t1\n";
    out<<"1\t1\t0\t0\t1\t1\t0\t1\t0\n";
    out<<"0\t0\t0\t1\t1\t1\t0\t0\t1\n";
    out<<"1\t1\t0\t0\t0\t0\t1\t1\t1\n";
    out<<"0\t0\t1\t1\t1\t0\t0\t0\t0\n";
    out<<"1\t1\t1\t0\t0\t0\t0\t0\t1\n";
    out<<"0\t0\t1\t0\t0\t1\t0\t1\t1\n";
	out.close();
    
    gcObs.m_listTrajFiles[0] = TCF1;
    gcObs.m_listTrajFiles[1] = TCF2;
    
    gcObs.m_numPermutations = 10;
    
    TrajectoryCollection trajCol;
    trajCol.scan(TCF1);
    
    gcObs.m_pathways.resize(1);
    gcObs.m_pathways[0].setNodeNames(trajCol.getIntNodeNames());
    gcObs.m_pathways[0].setNumConditions(2);
    
    //creating pathway
    string pathway = "pwy.txt";
    
	out.open(pathway.c_str(),ios::out);
	out<<"Pathway"<<endl;
	out<<"4\n";
    out<<"Node4*Node3,\n";
    out<<"Node3*Node2,\n";
    out<<"Node2*Node1,\n";
    out<<"Node1*Node4,\n";
	out.close();
    
    
    gcObs.m_pathways[0].scan(pathway);
    gcObs.m_pathways[0].setFile(pathway);
    
    gcObs.m_BP = trajCol.getBeParent();
    gcObs.m_HP = trajCol.getHaveParent();

    EnvGLNCmp envObs = CompareGLNsResample(gcObs, "BY_BOTH");
    
    MultiplePathwayStats pwStats = envObs.getPathwayStats();
    
    cout << endl << endl << "Obtained by independent simulation study -VERSUS- "
    << "Gamma approx with bootstrap mean and variance:" << endl;
    
    cout<<"pD (Simulation): " << 0.192960 << " -VERSUS- pD (Gamma approx): "
    << pwStats.m_het.m_pval << endl;
    cout<<"pC (Simulation):" << 0.369008 << " -VERSUS- pC (Gamma approx): "
    << pwStats.m_hom.m_pval << endl;
    cout<<"pT (Simulation):" << 0.269834 << " -VERSUS- pT (Gamma approx): "
    << pwStats.m_tot.m_pval << endl;
    cout<<"child pZ (Simulation):" << 0.0235999 << " -VERSUS- child pZ (Gamma approx): "
    << pwStats.m_wzChildren.m_pval << endl;
    cout<<"parent pZ (Simulation):" << 0.0345629 << " -VERSUS- parent pZ (Gamma approx): "
    << pwStats.m_wzParents.m_pval << endl <<endl;
    
    remove(TCF1.c_str());
    remove(TCF2.c_str());
    remove(pathway.c_str());
    remove("pwy-Conserved.dot");
    remove("pwy-Diff.dot");
    remove("pwy-WkZn.dot");
    remove("pwy-PWStats.txt");
}

void test_CompareGLNsResample_YZ()
{
    cout << ">>>> Testing Resample ChiNet ..." << endl;
    
    GLNConfig gcObs;
    gcObs.m_alpha = 1;
    gcObs.m_min_Markov_order = 0;
    gcObs.m_max_Markov_order = 0;
    gcObs.m_methodCmpxParCmp = BY_EACH_COND;
    
    gcObs.m_listTrajFiles = vector<string>(2);
    
    //creating TRJ1.trj
    string TCF1 = "TRJ1.trj";
    
	ofstream out(TCF1.c_str(),ios::out);
	out<<"TRAJECTORY_VER2"<<endl;
	out<<"1\t4\t0\n";
    out<<"2\t2\t2\t2\n";
    out<<"Node1\tNode2\tNode3\tNode4\n";
    out<<"10\n";
    out<<"1\t0\t1\t0\t0\t1\t1\t0\t1\n";
    out<<"0\t1\t1\t1\t1\t1\t0\t1\t0\n";
    out<<"1\t1\t0\t1\t0\t0\t0\t1\t1\n";
    out<<"0\t1\t0\t0\t1\t1\t1\t1\t1\n";
    out<<"1\t1\t1\t0\t0\t0\t1\t0\t0\n";
    out<<"1\t0\t0\t1\t1\t1\t1\t0\t0\n";
    out<<"1\t1\t1\t0\t0\t0\t0\t1\t1\n";
    out<<"0\t0\t0\t1\t1\t1\t0\t0\t0\n";
    out<<"1\t1\t1\t1\t0\t0\t0\t0\t0\n";
    out<<"1\t0\t0\t1\t0\t0\t1\t0\t1\n";
	out.close();
    
    //creating TRJ2.trj
    string TCF2 = "TRJ2.trj";
    
	out.open(TCF2.c_str(),ios::out);
	out<<"TRAJECTORY_VER2"<<endl;
	out<<"1\t4\t0\n";
    out<<"2\t2\t2\t2\n";
    out<<"Node1\tNode2\tNode3\tNode4\n";
    out<<"10\n";
    out<<"0\t1\t0\t1\t1\t1\t0\t1\t1\n";
    out<<"1\t1\t1\t1\t1\t0\t1\t0\t1\n";
    out<<"0\t0\t1\t0\t0\t1\t0\t1\t1\n";
    out<<"1\t0\t0\t1\t1\t1\t1\t1\t1\n";
    out<<"1\t1\t0\t0\t1\t1\t0\t1\t0\n";
    out<<"0\t0\t0\t1\t1\t1\t0\t0\t1\n";
    out<<"1\t1\t0\t0\t0\t0\t1\t1\t1\n";
    out<<"0\t0\t1\t1\t1\t0\t0\t0\t0\n";
    out<<"1\t1\t1\t0\t0\t0\t0\t0\t1\n";
    out<<"0\t0\t1\t0\t0\t1\t0\t1\t1\n";
	out.close();
    
    gcObs.m_listTrajFiles[0] = TCF1;
    gcObs.m_listTrajFiles[1] = TCF2;
    
    gcObs.m_numPermutations = 10;
    
    TrajectoryCollection trajCol;
    trajCol.scan(TCF1);
    
    gcObs.m_pathways.resize(1);
    gcObs.m_pathways[0].setNodeNames(trajCol.getIntNodeNames());
    gcObs.m_pathways[0].setNumConditions(2);
    
    //creating pathway
    string pathway = "pwy.txt";
    
	out.open(pathway.c_str(),ios::out);
	out<<"Pathway"<<endl;
	out<<"4\n";
    out<<"Node4*Node3,\n";
    out<<"Node3*Node2,\n";
    out<<"Node2*Node1,\n";
    out<<"Node1*Node4,\n";
	out.close();
    
    
    gcObs.m_pathways[0].scan(pathway);
    gcObs.m_pathways[0].setFile(pathway);
    
    gcObs.m_BP = trajCol.getBeParent();
    gcObs.m_HP = trajCol.getHaveParent();
    
    EnvGLNCmp envObs;
    
    envObs.initialize(gcObs);
    
    size_t nBootstraps = gcObs.m_numPermutations;
    
    ResampleType strategy = BOOTSTRAP_PLUS_ONE;
    
    vector<MultiplePathwayStats> stats(nBootstraps);
    
    for(size_t i=0; i<nBootstraps; i++) {
        
        EnvGLNCmp env;
        
        // Generate bootstrap trajectory collections
        // on master node only for MPI.
        vector<TrajectoryCollection> trajCols
        = resample(envObs.getTrajCols(), strategy);
        
        assert(trajCols.size() > 0);
        
        GLNConfig gc = gcObs;
        
        if (strategy == BOOTSTRAP_PLUS_ONE) {
            gc.m_max_Markov_order = -1;
            gc.m_min_Markov_order = -1;
        }
        
        env.initialize(gc);
        
        env.setTrajCols(trajCols);
        env.compareChildWorkingZones(); // due to changed trajectories
        CompareGLNs(env);
        
        // Record bootstrap pathway statistics from envBS
        stats[i] = env.getPathwayStats();
    }
    
    CompareGLNs(envObs);
    
    // Approximate null distribution and calculate p-value
    MultiplePathwayStats pwStats = envObs.getPathwayStats();
    
    pwStats.approximate(stats);
    
    cout<<endl<<endl<<"Expected pd: " << 0.192960 << " Empirical pd: " << pwStats.m_het.m_pval << endl;
    cout<<"Expected pc:" << 0.369008 << " Empirical pc: " << pwStats.m_hom.m_pval << endl;
    cout<<"Expected pt:" << 0.269834 << " Empirical pt: " << pwStats.m_tot.m_pval << endl;
    cout<<"Expected child pz:" << 0.0235999 << " Empirical child pz: " << pwStats.m_wzChildren.m_pval << endl;
    cout<<"Expected parent pz:" << 0.0345629 << " Empirical parent pz: " << pwStats.m_wzParents.m_pval << endl <<endl;
    
}