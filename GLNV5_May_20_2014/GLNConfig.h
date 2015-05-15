// GLNConfig.h
//
// Joe Song
// Created: September 27, 2008
// Modified:
//   October 1, 2011. Joe Song. Changed file name from
//      "GLN_EstParameters.h" to "GLNConfig.h"
//   June 10, 2013. MS.

#pragma once

#include <iostream>
using std::cerr;
using std::endl;

#include <string>
using std::string;

#include "PermStatTable.h"
#include "Topology.h"
#include "ComparativeChisq.h"

//added by Yang Zhang 2011.10.6 only compare per degree improvement for subset parents
//otherwise still use pvalue for selection
enum REC_CHISQ_CMP_METHOD {BY_PVAL, BY_PVAL_IMP, BY_PVAL_F_OVER_CHISQ,
    BY_IMP_PER_DF, BY_IMP_AND_PVAL};

enum CMPX_CHISQ_CMP_METHOD {BY_INTX_TYPE, BY_TOT, BY_HET, BY_HOM,
    BY_T_IMP_PER_DF, BY_T_C_IMP, BY_EACH_COND, BY_BEST_HOM};

enum CMPX_TOPOLOGY {SAME=1, DIFFERENT=2};

class GLNConfig
{
public:
    
    GLNConfig(void):m_mustIncludeParents()
    {        
        // int pval_mode=3;	 // method to calculate p-value, default = mode 3
        m_pValueMode = 3;

        // double alpha= .05; // max false positive rate, default = .05 
        m_alpha = 0.05;

        // int max_parents_allowed = -1; //Max num of parents for estimation or generation
        m_max_parents_allowed = -1; // 3;
        
        // int min_parents_allowed = 1;	// Joe Song, September 5, 2008        
        m_min_parents_allowed = 0; 
            // default to zero necessary for the comparison mode 
            // 0 or 1 is the same for the estimation mode
 
        // int compareMethod=1;	// using 2 in comparative program we can compare
                        //   different parents, added by Yang Zhang 2/19/2010
        m_acceptableTopologies = SAME;

        // int min_Markov_order = -1, // max parent time offset initialized to -1
        m_max_Markov_order = -1;

        // max_Markov_order = -1;
        m_min_Markov_order = -1;

        // int startAtNode = 1;
        m_startAtNode = 1;
        
        // bool allowSelfCycle = false;
        m_allowSelfCycle = false;
                
        // int numberOfPermutations = 0; //to specify number of permutations, added by Yang Zhang 7/26/2010
        m_numPermutations = 0;
        
        m_methodRecParCmp = BY_PVAL; // BY_IMP_PER_DF
        m_methodCmpxParCmp = BY_TOT; // BY_EACH_COND;
        m_methodMarginal = INDEP_SUB; // INDEP
        m_verbose = true;
    }
    
	bool allowSelfCycle () const
	{	
		return m_allowSelfCycle;
	}

    void setMethodCmpParSets(const string & method, const string & mode) 
    {
        if(mode == "estimation") {
            
            if ( method.empty() ) {
                
                // use default value declared in the constructor
                
            } else if(method == "BY_PVAL") {
                
                m_methodRecParCmp = BY_PVAL;
                
            } else if(method == "BY_PVAL_IMP") {
                
                m_methodRecParCmp = BY_PVAL_IMP;

            } else if(method == "BY_PVAL_F_OVER_CHISQ") {
            
                m_methodRecParCmp = BY_PVAL_F_OVER_CHISQ;
            
            } else if(method == "BY_IMP_PER_DF") {
            
                m_methodRecParCmp = BY_IMP_PER_DF;
            
            } else if(method == "BY_IMP_AND_PVAL") {

				m_methodRecParCmp = BY_IMP_AND_PVAL;

			} else {
            
                cerr << "ERROR: Invalid parent set comparison method \"" << method << "\" under \"" << mode << "\" mode!" << endl;
            }
            
        } else if(mode == "comparison") {
            
            if ( method.empty() ) {
                
                // use default value declared in the constructor
                
            } else if(method == "BY_INTX_TYPE") {
                
                m_methodCmpxParCmp = BY_INTX_TYPE;
                
            } else if(method == "BY_TOT") { 

                m_methodCmpxParCmp = BY_TOT;

            } else if(method == "BY_HET") {
                
                m_methodCmpxParCmp = BY_HET;
                
            } else if(method == "BY_HOM") {
                
                m_methodCmpxParCmp = BY_HOM;
                
            } else if(method == "BY_T_IMP_PER_DF") {
                
                m_methodCmpxParCmp = BY_T_IMP_PER_DF;
                
            } else if(method == "BY_T_C_IMP") {
                
                m_methodCmpxParCmp = BY_T_C_IMP;
                
            } else if(method == "BY_EACH_COND") {
                
                m_methodCmpxParCmp = BY_EACH_COND;

            } else if(method == "BY_BEST_HOM") {
                
                m_methodCmpxParCmp = BY_BEST_HOM;
                
            } else {
                cerr << "ERROR: Invalid parent set comparison method \""
                    << method << "\" under \"" << mode << "\" mode!" << endl;
            }
            
        } else { // other modes
            
            // this parameter is irrelevant
        }
    }
    
    void setMarginalMethod(const string & method)
    {
        if (method == "INDEP" || method == "INDEP_SUB") {
            m_methodMarginal = INDEP_SUB;
        } else if (method == "INDEP_SUP") {
            m_methodMarginal = INDEP_SUP;
        } else if (method == "POOLED") {
            m_methodMarginal = POOLED;
        } else if (method == "POOLED_NO_ROWCOL") {
            m_methodMarginal = POOLED_NO_ROWCOL;
        } else if(method == "STACKED") {
            m_methodMarginal = STACKED;
        } else if (method == "DIFF_COR_SUP") {
            m_methodMarginal = DIFF_COR_SUP;
		} else if (method == "CPFUNX2") {
			m_methodMarginal = CPFUNX2;
		}
		else {
            cerr << "ERROR: Invalid parent set comparison method \"" << method
            << "!" << endl;
        }
    }
    
public:

	int m_pValueMode;

	int m_max_parents_allowed;
	int m_min_parents_allowed; // Joe Song September 5, 2008

	int	m_min_Markov_order;
	int m_max_Markov_order;

	double m_alpha;

	int m_startAtNode; // 1-based

	bool m_allowSelfCycle; // Whether self cycle is allowed

	string m_dataType; // Added. Joe Song. September 5, 2008

	CPermStatTable m_permStatTable;
	
	vector<string> m_listTrajFiles; // Input files
	
	vector<Topology> m_pathways; //added by yangzhang 1.25.2009
    
	vector<Topology> m_candidateTopologys; //added by Yang Zhang 9/12/2009

	//added by Yang Zhang 9/15/2010
	vector<bool> m_HP;  //nodes can have parents or not
	vector<bool> m_BP;	//nodes can be parents or not used in JobGLN.cpp

	// Output files
	string m_fileGLN;
    
	string m_fileDOT; // string dotNetworkFile;
	string m_adjustedDOTFile; // string adjusteddotfile

	//added by YangZhang 2.23.2009
	string m_recordResultFile; 	// string recordResultFile; //added by yangzhang 2/23/2009 to record all chi-square test with p-value less than threshhold

	string m_adjustedpValFile; 	// string adjustedpvaluefile;

    // Interaction comparison file 
	string m_fileIntxCmp;    // string fileIntxCmp; // Joe Song. A file of a table format to save the interaction comparison results

    // Comparative modeling topology file
	string m_fileTopology; // string fileTopology; // Joe Song.  Dec 8, 2009. A file of the pathway format for saving the topology of comparative modeling

	// added by Yang Zhang 9.4.2010
    int m_numPermutations;
    
    // -----------------------------------------------------------------
    // Reconstruction mode: method to compare parent sets    
    REC_CHISQ_CMP_METHOD m_methodRecParCmp;
    
    // -----------------------------------------------------------------
    // Comparative mode: method to compare parent sets
    CMPX_CHISQ_CMP_METHOD m_methodCmpxParCmp;  

    CMPX_TOPOLOGY m_acceptableTopologies; //corresponds to different stretegy,
      // 2 will consider different parents when comparing while 1 (default) can not

    CMPX_MARGINAL m_methodMarginal;
    
    // -----------------------------------------------------------------
    //Added by Haizhou Wang, Jan 29, 2013
    Topology m_mustIncludeParents;
    //End of Haizhou's Code
    
    string m_filePSSN; // Proxy signature subnetwork classifier data file
    bool m_verbose; // a switch to turn on/off verbose output
};
