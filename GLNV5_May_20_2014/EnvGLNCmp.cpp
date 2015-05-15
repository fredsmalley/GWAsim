// EnvGLNCmp.cpp -- The job environment for transition and enumeration based GLN comparison
//                Multiple trajectories can be compared
//
// Joe Song
// Created: November 30, 2008
// Modified: November 2, 2009. Added working zone change code.
// Updated: April 7, 2011. Joe Song. Updated finalize() to the latest version
//    developed on March 14, 2010, where some redundant code is removed.
// Modified: September 13, 2011.  Introduced new modes of comparing different parent sets in comparative interaction analysis
// Last modified: September 29, 2011.  Added a function to compute parent working zone change.  The major change is that when 
//   parent identities are different, we consider there is a signficiant working zone change using p_z=0, chisq_z=0, and df_z=0. 
//   Before the change, the parent working zone is considered unchanged when the number of parents are different regardless of 
//   parent  identities

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#include "GLNGlobals.h"
#include "EnvGLNCmp.h"
#include "StatDistributions.h"
#include "Topology.h"
#include "PathwayStats.h"
#include "ComparativeChisq.h"

size_t EnvGLNCmp::getMaxNumDiffParents() const 
{    
    size_t n=0;
    switch (m_gc.m_acceptableTopologies) {
        case SAME: // 1:
            n = m_gc.m_max_parents_allowed;
            break;
        case DIFFERENT: //2:
            n = m_gc.m_max_parents_allowed * m_gc.m_listTrajFiles.size();
            break;
        default:
            GLNExit(EXIT_FAILURE, "m_acceptableTopologies is invalid!");
            break;
    }
    // To make sure that each node can have any combination of m_max_parents_allowed 
    return n; 
}

bool EnvGLNCmp::decideInteractionType(int child, const vector< int > & parents,
                                      const vector< int > & delays)
// This member function applies only if the parents are the SAME across all conditions 
{
	bool improved = false;

	// vector<Node> m_node = m_glnPooled.getNodes();  MS Dec 4, 2011

	size_t i;

	for(i = 0; i<parents.size(); i++) {
		// Check if child itself at delay=0 (current time) has been picked as its own parent.
		if(child == parents[i]-1 && delays[i] == 0) {
			break;
		}
	}

	if(i == parents.size() && i>0) { 
		// Only process an enumeration if the parents do not include 
		//   the child itself at delay = 0

		vector<TransitionTable> tts; // (m_trajCols.size());

        /*
        TransitionTable ttHom(child, parents, delays, m_trajCols[0]);

        // Fill in tts using the exact parent set for each condition
        tts = computeIndividualChisq(child,
                                     vector< vector<int> > (m_trajCols.size(), parents),
                                     vector< vector<int> > (m_trajCols.size(), delays),
                                     m_trajCols,
                                     m_gc.m_pValueMode,
                                     ttHom);
        
		TransitionTable ttHet;

        TransitionTable ttTot 
        = computeTotalChisq(child, 
                            vector< vector<int> > (m_trajCols.size(), parents), 
                            vector< vector<int> > (m_trajCols.size(), delays), 
                            m_trajCols, tts, m_gc.m_pValueMode);
        
		// ttHet = applyHeteroChisqTest(tts, ttHom, m_gc.m_pValueMode);
        */
        
        TransitionTable ttHom, ttHet, ttTot;
        
        collectStats(child, vector< vector<int> > (m_trajCols.size(), parents),
                     vector< vector<int> > (m_trajCols.size(), delays),
                     tts, ttTot, ttHom, ttHet);
        
        CmpIntxType type = labelIntxType(ttTot.getpValue(), ttHom.getpValue(),
                                         ttHet.getpValue(), m_gc.m_alpha);
        
		improved = compareParentSets(child, type, tts, ttTot, ttHom, ttHet);
		
		if(improved)
		{
			// Compute work zone change statistics
			double chisq_z, p_z;
			size_t df_z;

			chisq_z = applyKSampleRowChisqTest(tts, df_z, p_z);

			updateEnvironment(child, type, tts, ttTot, ttHom, ttHet, p_z, (int) df_z, chisq_z);
		}
		//save all results
		if(m_gc.m_recordResultFile.length()!=0)
			recordResults(child, ttHet, ttHom, ttTot, string("Differential")); // MS Dec 4, 2011. m_glnPooled.getNodes(),
	}
	return improved;
}

bool EnvGLNCmp::analyze(int child, const vector<vector< int > > & parents, 
                        const vector<vector< int > > & delays)
// Compares interactions regardless of same or different parents across conditions 
{
	bool improved = false;
    
	if( ! isSelfCycleZeroDelay(child, parents, delays) ) {

        TransitionTable ttTot, ttHom, ttHet;
        vector<TransitionTable> tts;
        
        collectStats(child, parents, delays, tts, ttTot, ttHom, ttHet);

        improved = select(child, tts, ttTot, ttHom, ttHet);
        
		if(! m_gc.m_recordResultFile.empty() ) // save all results
			recordResults(child, ttHet, string("Differential")); // MS Dec 4, 2011 m_glnPooled.getNodes(),

	}
	return improved;
}

void EnvGLNCmp::updateEnvironment(int child, int type, 
                                  const vector <TransitionTable> & tts, 
                                  const TransitionTable & ttTot, 
                                  const TransitionTable & ttHom, 
                                  const TransitionTable & ttHet, 
                                  double p_z, int df_z, double chi_z)
{
	m_pooledTransTables[child] = ttHom; 
	m_diffTransTables[child] = ttHet;
	m_TransTables[child] = tts;

	m_chisqs_t[child] = ttTot.getChisq();
	m_dfs_t[child] = ttTot.getDf();
	m_pchisqs_t[child] = ttTot.getpValue();

	m_chisqs_z[child] = chi_z;
	m_dfs_z[child] = df_z;
	m_pchisqs_z[child] = p_z;

	m_types[child] = type;
}


bool EnvGLNCmp::
processOneEnumeration(const int child, const vector< int > & parents,
                      const vector< int > & delays)
{
	bool better = false;
	better = decideInteractionType(child, parents, delays);
	return better;
}

bool EnvGLNCmp::processOneEnumeration(const int child, const vector<vector< int > > & parents, 
                                      const vector<vector< int > > & delays)
{
	bool better = false;
	better = analyze(child, parents, delays); // decideInteractionType(child, parents, delays);
	return better;
}

void EnvGLNCmp::initialize(const GLNConfig & gc)
{
	// Estimation parameter initialization
	m_gc = gc;

	// Load trajectory collection from a list of given files
	m_trajCols.resize( m_gc.m_listTrajFiles.size() );

	for(size_t i=0; i<m_gc.m_listTrajFiles.size(); i++) {

		m_trajCols[i].scan(m_gc.m_listTrajFiles[i]);
        
        if(i > 0) { // validate the compatibility of the trajectory collection files

            ostringstream oss;
            
            if(m_trajCols[i].nNodes() != m_trajCols[i-1].nNodes()) {
                                
                oss << m_gc.m_listTrajFiles[i] << " and " << m_gc.m_listTrajFiles[i]
                    << " do not have the same number of nodes!" << endl;
                GLNExit(EXIT_FAILURE, oss.str().c_str());
                
            } else if(m_trajCols[i].getBases() != m_trajCols[i-1].getBases()) {
                
                oss << m_gc.m_listTrajFiles[i] << " and " << m_gc.m_listTrajFiles[i]
                << " do not have the same radix for each node!" << endl;
                GLNExit(EXIT_FAILURE, oss.str().c_str());

            } else if(m_trajCols[i].getExtNodes() != m_trajCols[i-1].getExtNodes()) {
                
                oss << m_gc.m_listTrajFiles[i] << " and " << m_gc.m_listTrajFiles[i]
                << " do not have the same external nodes!" << endl;
                GLNExit(EXIT_FAILURE, oss.str().c_str());
                
            } else if(m_trajCols[i].getExtNodeNames() != m_trajCols[i-1].getExtNodeNames()) {
                
                oss << m_gc.m_listTrajFiles[i] << " and " << m_gc.m_listTrajFiles[i]
                << " do not have the same external node names!" << endl;
                GLNExit(EXIT_FAILURE, oss.str().c_str());

            } else if(m_trajCols[i].getIntNodeNames() != m_trajCols[i-1].getIntNodeNames()) {
                
                oss << m_gc.m_listTrajFiles[i] << " and " << m_gc.m_listTrajFiles[i]
                << " do not have the same internal node names!" << endl;
                GLNExit(EXIT_FAILURE, oss.str().c_str());
                
            }
            
        }

	}

	// Initialize the GLN
	//m_glnPooled.initialize( m_trajCols[0].nNodes(), m_trajCols[0].getExtNodes(), m_trajCols[0].getExtNodeNames(), m_trajCols[0].getIntNodeNames() );
	m_glnPooled.initialize( m_trajCols[0].nNodes(), 
		m_trajCols[0].getExtNodes(), 
		m_trajCols[0].getExtNodeNames(), 
		m_trajCols[0].getIntNodeNames(),
		m_trajCols[0].getBeParent(), 
		m_trajCols[0].getHaveParent());

	m_glnDiff = m_glnPooled;

    setMarkovOrder(gc.m_min_Markov_order, gc.m_max_Markov_order);
    
	if(m_gc.m_max_parents_allowed == -1) {
		 m_gc.m_max_parents_allowed = (int) m_glnPooled.size();
	 }

	// Initialize final result variables used in estimation
	m_diffTransTables.resize( m_trajCols[0].nNodes() );
	m_pooledTransTables.resize( m_trajCols[0].nNodes() );

	for(int i=0; i < (int) m_trajCols[0].nNodes(); i++)
	{
		m_diffTransTables[i].setChild( i, m_trajCols[0].base( i ), m_glnDiff.getNodeName(i) );
		m_pooledTransTables[i].setChild( i, m_trajCols[0].base( i ), m_glnPooled.getNodeName(i) );
	}

	m_recordresult.resize((int) m_glnPooled.size());
	m_recordcounts.resize((int) m_glnPooled.size());
	
	//added by Yang Zhang 5/9/2009
	m_TransTables.resize( m_trajCols[0].nNodes() );

	m_chisqs_t.resize( m_trajCols[0].nNodes() );
	m_dfs_t.resize( m_trajCols[0].nNodes() );
	m_pchisqs_t.resize( m_trajCols[0].nNodes() );

	m_types.resize( m_trajCols[0].nNodes() ); // interaction type vector

	// Parent working zone statistics:
	m_chisqs_z.resize( m_trajCols[0].nNodes() );
	m_dfs_z.resize( m_trajCols[0].nNodes() );
	m_pchisqs_z.resize( m_trajCols[0].nNodes() );
	m_workingZones.resize( m_trajCols[0].nNodes() );
	
	for(int i=0; i < (int) m_trajCols[0].nNodes(); i++)
	{
		m_types[i] = NULL_INTX;  //0;
		m_workingZones[i] = UNCHANGED; 

		m_chisqs_t[i] = 0;
		m_dfs_t[i] = 0;
		m_pchisqs_t[i] = 1;

		m_chisqs_z[i] = 0;
		m_dfs_z[i] = 0;
		m_pchisqs_z[i] = 1;
	}

	// Child working zone statistics:
	//added by Yang Zhang 11/08/09
	m_childchisqs_z.resize( m_trajCols[0].nNodes() );
	m_childdfs_z.resize( m_trajCols[0].nNodes() );
	m_childpchisqs_z.resize( m_trajCols[0].nNodes() );
	m_childWorkingZones.resize( m_trajCols[0].nNodes() );

	compareChildWorkingZones();
}

void EnvGLNCmp::compareChildWorkingZones()
{
	for(int i=0; i < (int) m_trajCols[0].nNodes(); i++)
	{
		m_TransTables[i].resize(m_trajCols.size());
		for(size_t j=0; j<m_trajCols.size();j++)
		{
			m_TransTables[i][j].setChild(i, m_trajCols[0].base( i ), m_glnPooled.getNodeName(i));
		}

		// Differential gene expression analysis: compute child working zone change

		vector<TransitionTable> tts(m_trajCols.size());
		
        vector<int> delays (1, 0);
		vector<int> parents (1, i+1);

		for(size_t j=0; j<m_trajCols.size(); j++) {
			
			tts[j].initialize(i, parents, delays, m_trajCols[j]);
			tts[j].fill(m_trajCols[j]);

		}

		double chisq_z, p_z;
		size_t df_z;

		// chisq_z = applyKSampleChisqTest(tts, df_z, p_z);
        computeChildWorkingZoneChange(tts, chisq_z, df_z, p_z);
        
		m_childchisqs_z[i] = chisq_z;

		m_childdfs_z[i] = (int)df_z;

		m_childpchisqs_z[i] = p_z;

		m_childWorkingZones[i] = p_z <= m_gc.m_alpha ? CHANGED : UNCHANGED;

	}

}

void EnvGLNCmp::gatherPathwayStats()
{
    Topology topo;
    
    if(m_gc.m_pathways.size()>=1)
    {
        topo = m_gc.m_pathways[0];
    }else if(m_gc.m_pathways.size() == 0 && m_gc.m_candidateTopologys.size()>=1) {
        topo = m_gc.m_candidateTopologys[0];
    }
	
    if (topo.size() == 0) {
        return;
    }
    
    m_pathwayStats.compare(topo, m_diffTransTables, m_pooledTransTables,
                           m_chisqs_t, m_dfs_t, m_chisqs_z, m_dfs_z,
                           m_childchisqs_z, m_childdfs_z);

    // pwStats.scale(m_trajCols, topo.getAllNodes());
        
}

void EnvGLNCmp::savePathwayStats() const
{
    string file;
    
    if(m_gc.m_pathways.size()==1)
    {
     file = m_gc.m_pathways[0].getFile();
    }

    if(m_gc.m_candidateTopologys.size()==1)
    {
     file = m_gc.m_candidateTopologys[0].getFile();
    }

    if(m_gc.m_pathways.size()>1)
	{
        for(size_t i=0; i<m_gc.m_pathways.size(); i++)
        {
            if((i<m_gc.m_pathways.size()-1)&&m_gc.m_pathways[i].getFile().length()>0) //no dash after last name
            {
                file = file + m_gc.m_pathways[i].getFile().substr(0, m_gc.m_pathways[i].getFile().find(".txt")) + "-";
            }else
            {
                file = file + m_gc.m_pathways[i].getFile();
            }
        }
	}else if(m_gc.m_candidateTopologys.size()>1){
        
        for(size_t i=0; i<m_gc.m_candidateTopologys.size(); i++)
        {
            if((i<m_gc.m_candidateTopologys.size()-1)&&m_gc.m_candidateTopologys[i].getFile().length()>0)  //no dash after last name
            {
                file = file + m_gc.m_candidateTopologys[i].getFile().substr(0, m_gc.m_candidateTopologys[i].getFile().find(".txt")) + "-";
            }else
            {
                file = file + m_gc.m_candidateTopologys[i].getFile();
            }
        }
	}
    
    if(file.empty()) {
        return;
    }

    string base = file.substr(0, file.find(".txt"));

    string filename = base + "-PWStats.txt";

    m_pathwayStats.save(filename);
}

void EnvGLNCmp::comparePathway()
{

	//added by Yang Zhang 2010.1.18

	if(m_gc.m_pathways.size()==1
	|| m_gc.m_candidateTopologys.size()==1)
	{
        /*
		string filename = m_gc.m_pathway.getFile();
		if(m_gc.m_candidateTopology.size()>0)
		{
			filename = m_gc.m_candidateTopology.getFile();
		}
		filename = filename.substr(0,filename.find(".txt"));
		filename = filename + "result.txt";
        */
		string base;
        if(m_gc.m_pathways.size()==1)
		{
			base = m_gc.m_pathways[0].getFile().substr(0, m_gc.m_pathways[0].getFile().find(".txt"));
		}else{
			base = m_gc.m_candidateTopologys[0].getFile().substr(0, m_gc.m_candidateTopologys[0].getFile().find(".txt"));
		}

        string filename = base + "-PWStats.txt";
                
        MultiplePathwayStats mps;

		mps.m_het.m_df = 0; // heterogeneity degree of freedom
		mps.m_het.m_x = 0.0;// heterogeity chisq in pathway
		for(int i=0;i < (int) m_diffTransTables.size();i++)
		{
			mps.m_het.m_df += m_diffTransTables[i].getDf();
			mps.m_het.m_x += m_diffTransTables[i].getChisq();
		}
        mps.m_het.m_pval = ChisqPvalue(mps.m_het.m_x, mps.m_het.m_df);
        
		mps.m_hom.m_df = 0; // homogeneity degree of freedom
		mps.m_hom.m_x = 0.0;// homogeneity chisq in pathway
		for(int i=0; i < (int) m_pooledTransTables.size(); i++)
		{
			mps.m_hom.m_df += m_pooledTransTables[i].getDf();
			mps.m_hom.m_x += m_pooledTransTables[i].getChisq();
		}
        mps.m_het.m_pval = ChisqPvalue(mps.m_hom.m_x, mps.m_hom.m_df);
        
		mps.m_tot.m_df = 0; // total degree of freedom
		mps.m_tot.m_x = 0.0;// total chisq in pathway
		for(int i=0; i < (int) m_chisqs_t.size(); i++)
		{
			mps.m_tot.m_df += m_dfs_t[i];
			mps.m_tot.m_x += m_chisqs_t[i];
		}
		mps.m_tot.m_pval = ChisqPvalue(mps.m_tot.m_x, mps.m_tot.m_df);
        
		mps.m_wzParents.m_df = 0; // pathway parent working zone degrees of freedom
		mps.m_wzParents.m_x = 0.0;// pathway parent working zone chisq in pathway
		for(int i=0; i < (int) m_chisqs_z.size(); i++)
		{
			mps.m_wzParents.m_df += m_dfs_z[i];
			mps.m_wzParents.m_x += m_chisqs_z[i];
		}
		mps.m_wzParents.m_pval = ChisqPvalue(mps.m_wzParents.m_x,
                                             mps.m_wzParents.m_df);
        
		mps.m_wzChildren.m_df = 0; // pathway child working zone degrees of freedom
		mps.m_wzChildren.m_x = 0.0; // pathway child working zone chisq
        
		//only calculate those node on pathway, but including nodes have parents and nodes have not
		//if(m_dfs_z[i] != 0) this maynot be true, when ith node have different # of parents, m_dfs_z[i] will be 0
		if(m_gc.m_pathways.size()==1)
		{
            //	if(m_gc.m_pathway.getPathway()[i].size() != 0)
			vector<int> NodesIdOnPathway = m_gc.m_pathways[0].getAllChildNodes();
			for(size_t index =0; index < NodesIdOnPathway.size(); index ++ )
			{
				mps.m_wzChildren.m_df += m_childdfs_z[NodesIdOnPathway[index]-1];
				mps.m_wzChildren.m_x += m_childchisqs_z[NodesIdOnPathway[index]-1];
			}
		}else{
			vector<int> NodesIdOnCandidateFile = m_gc.m_candidateTopologys[0].getAllChildNodes();
			for(size_t index =0; index < NodesIdOnCandidateFile.size(); index ++ )
			{
				mps.m_wzChildren.m_df += m_childdfs_z[NodesIdOnCandidateFile[index]-1];
				mps.m_wzChildren.m_x += m_childchisqs_z[NodesIdOnCandidateFile[index]-1];
			}
		}
        mps.m_wzChildren.m_pval = ChisqPvalue(mps.m_wzChildren.m_x, mps.m_wzChildren.m_df);
        
        mps.save(filename);
                 
		if(m_gc.m_numPermutations==0)
		{
			saveAllInteractionToDot();
		}

	}

}


void EnvGLNCmp::finalize()
{
    
#ifdef ENABLE_MPI
  //only write result file on master node, which is the head node with rank = 0
  int rank =0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank == 0) {
#endif

	for(int i=0; i < (int) m_glnPooled.size(); i++){
        
        // Create generalized truth tables
        
		GeneralizedTruthTable gtt;

		gtt.convert(m_pooledTransTables[i], m_gc.m_alpha);

		m_glnPooled.setGTT(i, gtt);

		m_glnPooled.initializeNode(i, m_trajCols[0].get(0, 0, i), m_trajCols[0].base(i), i+1);		//initialize the nodes in the network

		gtt.convert(m_diffTransTables[i], m_gc.m_alpha);

		m_glnDiff.setGTT(i, gtt);

		m_glnDiff.initializeNode(i, m_trajCols[0].get(0, 0, i), m_trajCols[0].base(i), i+1);		//initialize the nodes in the network

		m_workingZones[i] = (m_pchisqs_z[i] <= m_gc.m_alpha) ? CHANGED : UNCHANGED; 
	}

	// calculateOverallPval(m_gc.m_pValueMode);
	
	m_glnDiff.save(m_gc.m_fileGLN);
	
    switch (m_gc.m_acceptableTopologies) {
            
        case SAME: // 1:
            
            if(! m_gc.m_fileDOT.empty()) {
                saveDotFile(m_gc.m_fileDOT);
            }
            if( m_gc.m_verbose ) {
                print();
            }
            
            break;

        case DIFFERENT: // 2:

            if(! m_gc.m_fileDOT.empty()) {
                saveDOTFilebyCondition(m_gc.m_fileDOT); //saveDotFilek2(m_gc.m_fileDOT);
            }
            
            if( m_gc.m_verbose ) {
                printk2();
            }
            
            break;
            
        default:
            break;
    }

	save(m_gc.m_fileIntxCmp);
	
    gatherPathwayStats(); // comparePathway();

    savePathwayStats();
    
    if(m_gc.m_numPermutations==0) {
        saveAllInteractionToDot();
    }
    
    // Saving non-null and child and parent working zone changed interactions
    saveTopology(m_gc.m_fileTopology);

#ifdef ENABLE_MPI
  }
#endif
}

//modified by YangZhang 4/13/2009, 2/15/2009
//added by YangZhang 2/7/2009
//dflag denotes different flag
int EnvGLNCmp::checkdiff(const vector<TransitionTable> & diffTransTables,
                         const vector<TransitionTable> & pooledTransTables,
                         double alpha) const
{
	/*
	let dflag to denote whether or not two networks are different,
	0 means no common interactions,
	1 means share common interactions,
	2 means have different interactions
	*/

	int dflag = 0;
	for(size_t i = 0; i <diffTransTables.size(); i++)
	{
		double tdiffpvalue = diffTransTables[i].getpValue();//temp pvalue
		double tpoolpvalue = pooledTransTables[i].getpValue();//temp pvalue
		// double tdiffchisq = diffTransTables[i].getChisq();//temp chisq
		// double tpoolchisq = pooledTransTables[i].getChisq();//temp chisq
		// int tdiffdf = diffTransTables[i].getDf();//temp degree of freedom
		// int tpooldf = pooledTransTables[i].getDf();//temp degree of freedom
		
		if (tdiffpvalue<=alpha)
		{
			dflag = 2;
		}
		else
		{
			if(tpoolpvalue<=alpha)
			{
				dflag = 1;
			}
			else
			{
				dflag = 0;
			}
		}
		if(dflag==2) break;
	}
	return dflag;
}
