// EnvGLNRec.cpp -- The job environment for transition and enumeration based GLN reconstruction
//
// Joe Song
// Created: November 28, 2008

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>

using std::ios;
using std::cerr;
using std::endl;
using std::cout;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;

#include "GLNGlobals.h"
#include "EnvGLNRec.h"
#include "Topology.h"
#include "InteractionEvaluationController.h"
#include "PathwayStats.h"
#include "StatDistributions.h"

void adjusted_pvalue(GeneralizedLogicalNetwork & net, vector<TransitionTable> & tts, int numParents);

bool containedInCandidates(int childid, const vector<int> &parents);

void EnvGLNRec::initialize(const GLNConfig & gc)
{
	// Estimation parameter initialization
	m_gc = gc;
    
	// Load trajectory collection from a list of given files
	m_trajCol.scan(m_gc.m_listTrajFiles);
    
	// Initialize the GLN
    //m_gln.initialize( m_trajCol.nNodes(), m_trajCol.getExtNodes(), m_trajCol.getExtNodeNames(), m_trajCol.getIntNodeNames() );
	//modified by Yang Zhang 4/24/2009
	m_gln.initialize( m_trajCol.nNodes(), m_trajCol.getExtNodes(), m_trajCol.getExtNodeNames(),
                     m_trajCol.getIntNodeNames(),m_trajCol.getBeParent(), m_trajCol.getHaveParent());
	m_gln.setMaxMarkovOrder(m_gc.m_max_Markov_order);
    
    if(m_gc.m_max_parents_allowed == -1) {
        m_gc.m_max_parents_allowed = (int) m_gln.size();
    }

    if( (int) m_gln.size() < m_gc.m_min_parents_allowed )
        GLNExit( EXIT_FAILURE, "Error: minimum number of parents required is larger than the size of whole network.\n Condiser revise '-g' option.\n");
    if( m_gc.m_max_parents_allowed < m_gc.m_min_parents_allowed )
        GLNExit( EXIT_FAILURE, "Error: minimum number of parents required is larger than maximum number of parents allowed.\n Condiser revise '-p' and '-g' option.\n");

	// Initialize final result variables used in estimation
	m_transTables.resize( m_trajCol.nNodes() );
    
	for(int i=0; i < (int) m_transTables.size(); ++i)
	{
		m_transTables[i].setChild( i, m_trajCol.base( i ), m_gln.getNodeName(i) );
	}
	
	//resize recordresult
	m_recordresult.resize((int) m_gln.size());
	m_recordcounts.resize((int) m_gln.size());
}

//child starts from 0, parents start from 1
bool EnvGLNRec::processOneEnumeration(const int child, const vector< int > & parents, const vector< int > & delays)
{
	//added by Yang Zhang 9/12/2009 if parents not in candidate parents list then return
	//no need of this checking after changing parent enumeration strategy in EnumMethod
	/*if(m_gc.m_candidateTopology.getCandidateParents().size()>0)
     {
     CandidateParents candidates;
     candidates.setCandidateParents(m_gc.m_candidateTopology.getCandidateParents());
     if(!candidates.containedInCandidates(child,parents))
     {
     return false;
     }
     }*/
    
	bool better = false;
	//parent index starts from 1,child index starts from 0
	//if this node can not have parent,just return added by Yang Zhang
	// vector<Node> m_node = m_gln.getNodes();  MS Dec 4, 2011
    
	size_t i;
	for(i=0; i<parents.size(); ++i) {
		// Check if child itself at delay=0 (current time) has been picked as its own parent.
		if(child == parents[i]-1 && delays[i] == 0) {
			break;
		}
	}
    
	if( i==parents.size() && i>0 ) {
		// Only process an enumeration if the parents do not include
		//   the child itself at delay = 0
        
		TransitionTable tt(child, parents, delays, m_trajCol);
        
        if (m_gc.m_pValueMode > 100)
            InteractionEvaluationController::getController().evaluate(tt, m_trajCol);
        else{
            tt.fill(m_trajCol);
            
            switch(m_gc.m_pValueMode) {
                case 1:
                case 2:
                case 3:
                case 5:
                case 6:
                case 7:
				case 8:
                case 9:
                case 10:
                case 15:
                    tt.applyChisquareTest(m_gc.m_pValueMode);
                    break;
                case 4:
                    //forMultinomial(trans_table, m_trajCol.base(child), m_gc.m_pValueMode, chisq);
                    //p_val = (m_gc.m_permStatTable).compute_pValue(child, (int) k, m_gc.m_max_Markov_order,
                    //	m_gc.m_min_Markov_order, chisq);
                    tt.applyChisqPermTest(m_gc.m_pValueMode,m_gc.m_permStatTable,m_gc.m_max_Markov_order,m_gc.m_min_Markov_order);
                    break;
                default:
                    cerr << "ERROR: p-value mode must be 1 to 9!" << endl;
                    GLNExit(EXIT_FAILURE);
            }
        }
        
		if(tt.getpValue() <= m_gc.m_alpha) {
            
			better = compareParentSets(tt, m_transTables[child]);
			
			//added by YangZhang 2/22/2009 to record all parents with Chi-square test p-value less than alpha(say 0.05)
            if(m_gc.m_recordResultFile.length()!=0)
            {
                recordResults(child, tt);
			}
            
            if(better) {
				m_transTables[child] = tt;
			}
		}
	}
    
	return better;
}


//this function is added to be compatible with comparative GLN
//currently it is not utilized
bool EnvGLNRec::processOneEnumeration(const int child, const vector<vector< int > > & parents,
                                      const vector<vector< int > > & delays)
{
	//added by Yang Zhang 9/12/2009 if parents not in candidate parents list then return
	//no need of this checking after changing parent enumeration strategy in EnumMethod
	/*if(m_gc.m_candidateTopology.getCandidateParents().size()>0)
     {
     CandidateParents candidates;
     candidates.setCandidateParents(m_gc.m_candidateTopology.getCandidateParents());
     if(!candidates.containedInCandidates(child,parents[0]))
     {
     return false;
     }
     }*/
    
	bool better = false;
    
    // vector<Node> * p_nodes = & m_gln.getNodes();  // MS Dec 4, 2011.
    
	/*if this node can not have parent,just return added by Yang Zhang
     if (!m_node[child].getHP())
     {
     return better;
     }*/
    
	vector<int> canbeparents; //added by Yang Zhang 4/24/2009
	vector<int> canbedelays;	//added by Yang Zhang 4/27/2009
    
	for(size_t t=0; t <parents[0].size(); t++)
	{
		if(m_gln.getNode(parents[0][t]-1).getBP())
		{     //parent index starts from 1,child index starts from 0
			canbeparents.push_back(parents[0][t]);
			canbedelays.push_back(delays[0][t]);
		}
	}
    
	size_t i;
    
	for(i = 0; i<canbeparents.size(); i++) {
		// Check if child itself at delay=0 (current time) has been picked as its own parent.
		if(child == canbeparents[i]-1 && canbedelays[i] == 0) {
			break;
		}
	}
    
	if(i == canbeparents.size() && i>0 ) {
		// Only process an enumeration if the parents do not include
		//   the child itself at delay = 0
        
		TransitionTable tt(child, canbeparents, canbedelays, m_trajCol);
		tt.fill(m_trajCol);
        
		switch(m_gc.m_pValueMode) {
			case 1:
			case 2:
			case 3:
            case 5: // case 5 is added by Joe Song 9/11/2011
            case 10: // case 10 is added by Joe Song 2/28/2015
				tt.applyChisquareTest(m_gc.m_pValueMode);
				break;
			case 4:
				//forMultinomial(trans_table, m_trajCol.base(child), m_gc.m_pValueMode, chisq);
				//p_val = (m_gc.m_permStatTable).compute_pValue(child, (int) k, m_gc.m_max_Markov_order,
				//	m_gc.m_min_Markov_order, chisq);
				tt.applyChisqPermTest(m_gc.m_pValueMode,m_gc.m_permStatTable,m_gc.m_max_Markov_order,m_gc.m_min_Markov_order);
				break;
		}
        
		if(tt.getpValue() <= m_gc.m_alpha) {
            
			better = compareParentSets(tt, m_transTables[child]);
            
			if(m_gc.m_recordResultFile.length()!=0)
			{
                //added by YangZhang 2/22/2009 to record all parents with
                //  Chi-square test p-value less than alpha(say 0.05)
                
                recordResults(child, tt);
			}
            
			if(better) {
				m_transTables[child] = tt;
			}
		}
	}
    
	return better;
}


void EnvGLNRec::calculateOverallPval(int pValMode)
{
	switch(pValMode) {
            
		case 1:
		case 2:
		case 4:
		case 8:
        case 9:
        case 15:
			m_p_val=1;
			//Calculate p vavector <int> parentoffsetslue from the array of p values
			for(size_t i=0; i < m_transTables.size(); i++){
				if( m_transTables[i].getParents().size() > 0 ) {  // Joe Song 08/05/06
					m_p_val = m_p_val * (1.0 - m_transTables[i].getpValue());
				}
			}//end for
			m_p_val = 1 - m_p_val;
            
			break;
            
		case 3:
		case 5: // Joe Song. Added October 1, 2010
        case 10: // case 10 is added by Joe Song 2/28/2015
		case 6: // Tyler Hunt added September 22, 2011
        {
            m_chisq = 0;
            int df = 0;
            for(size_t i=0; i < m_transTables.size(); i++){
                // accumulate chisq
                if(m_transTables[i].getParents().size() > 0) {
                    m_chisq += m_transTables[i].getChisq();
                    df += m_transTables[i].getDf();
                }
            }
            //m_p_val = qchisq(df, m_chisq);
            m_p_val = ChisqPvalue(m_chisq, df);
        }
			break;
        case 101:
        case 102:
        case 103:
        case 104:
        case 105:
        case 106:
        case 107:
		case 7:
			m_p_val = .5; // not sure how to handel this ill set it to .5 for now
			break;
		default:
			GLNExit(EXIT_FAILURE, "ERROR: Invalid p-value mode!");
	}
}

void EnvGLNRec::calculatePathwayChisqStats() const
{
	//added by YangZhang 06/12/2009
	if(m_gc.m_pathways.size()==1)
	{
		string filename = m_gc.m_pathways[0].getFile();
		filename = filename.substr(0,filename.find(".txt"));
		filename = filename + "result.txt";
        
        /* MS 4/22/2013
         ofstream out(filename.c_str());
         //out.precision(2);
         out<<"";
         out.close();
         */
        
		int df = 0;//total degree of freedom
		double chisquare = 0;//total chisq in pathway
		for(int i=0;i<(int)m_transTables.size();i++)
		{
			df += m_transTables[i].getDf();
			chisquare += m_transTables[i].getChisq();
		}
        
        // MS 4/22/2013
        SinglePathwayStats sps;
        
        sps.m_tot = Chisq(chisquare, df, ChisqPvalue(chisquare, df));
        sps.save(filename);
        
	}
    
}

void EnvGLNRec::finalize()
{
	calculateOverallPval(m_gc.m_pValueMode);
    
	for(int i=0; i < (int) m_gln.size(); i++) {
        //initialize the nodes in the network
        
		GeneralizedTruthTable gtt;
		
		gtt.convert(m_transTables[i], m_gc.m_alpha);
        
		m_gln.setGTT(i, gtt);
        
		m_gln.initializeNode(i, m_trajCol.get(0, 0, i), m_trajCol.base(i), i+1);
        
	}//end for
    
    if(! m_gc.m_fileGLN.empty()) {
        m_gln.save(m_gc.m_fileGLN);
    }
    
	if(! m_gc.m_fileDOT.empty() ) {
		m_gln.saveDot(m_gc.m_fileDOT);
	}
    
	adjusted_pvalue(m_gln, m_transTables, m_gc.m_max_parents_allowed);
	
	if(! m_gc.m_adjustedpValFile.empty() ) {
		m_gln.save(m_gc.m_adjustedpValFile);
	}
    
	if(! m_gc.m_adjustedDOTFile.empty() ) {
		m_gln.saveDot(m_gc.m_adjustedDOTFile);
	}
    
    calculatePathwayChisqStats();
    
    if(m_gc.m_verbose) {
        print();
    }
    
    // Saving non-null and child and parent working zone changed interactions
    saveTopology(m_gc.m_fileTopology);

}
