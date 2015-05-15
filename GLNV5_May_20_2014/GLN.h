// GeneralizedLogicalNetwork.h
//
// Eric Lance, Chris Lewis, Curtis Luce, Joe Song
// Created: 2006
// Modified: 
//   September 7, 2008
//   Dec 4, 2011 (MS).  
//     - Deleted the member function getNodes() to hide data structure on how 
//       nodes are stored within the GeneralizedLogicalNetwork class. 
//     - Added several get functions to access information related to a node:
//       getNode(), getNodeName(), getNodeType()
//     - Added an empty print() function for future coding

#pragma once

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <utility>

#include "Node.h"
#include "PermStatTable.h"
#include "TrajectoryCollection.h"
#include "GeneralizedTruthTable.h"
#include "Topology.h"

class GeneralizedLogicalNetwork
{
protected:

	vector<Node> m_nodes;

	int m_min_Markov_order;
	int m_max_Markov_order;

	int m_current_time;

	string m_version;

	vector<GeneralizedTruthTable> m_gtts;  // generalized truth tables for each node

	bool m_ready;			//Has the network been initialized?

public:

	GeneralizedLogicalNetwork(void) { }

	// ~GeneralizedLogicalNetwork(void) {}

	int detectCycle(const vector< vector<int> > & trajectory) const;

	/* Initialization */
	void initialize(int networksize);
	void initialize(int networksize, const vector<int> & extNodes,
		const vector<string> & extNodeNames, const vector<string> & intNodeNames);
	//added by YangZhang 4/24/2009
	void initialize(int networksize, const vector<int> & extNodes,const vector<string> & extNodeNames, 
		const vector<string> & intNodeNames, const vector<bool> & beParent, const vector<bool> & haveParent);
	/*
	void initializeNode(int nodenum, unsigned short initialValue,
		unsigned short radix, int num_parents, vector<int> parentindexes,
		Predictor pred, double p_val = -1.0, string node_name = "!!", 
		unsigned char node_type = 'n');

	void initializeNode(int nodenum, unsigned short initialValue,
		unsigned short radix, int num_parents, vector<int> parentindexes,
		Predictor pred, vector<int> parentoffsets, double p_val = -1.0, 
		int node_id = -1, string node_name ="!!", unsigned char node_type = 'n');
	*/

	void initializeNode(int nodenum, unsigned short initialValue,
		unsigned short radix, string node_name = "!!", unsigned char node_type = 'n');

	void initializeNode(int nodenum, unsigned short initialValue,
		unsigned short radix, int node_id = -1, string node_name ="!!", unsigned char node_type = 'n');

	/* Random generation */
	void initialStates(string infile, string outfile, int numTraj, int min_Markov_order, int max_Markov_order);

	void generateRandomNetwork(int networkSize, int minBase, int maxBase, int minParents, int maxParents,
		int min_Markov_order, int max_Markov_order, string init_traj_file, int numb_of_init);
    
	/* I/O functions */
	void save(const string & filename);
	void scan(const string & filename);

	void saveEdgesDotFormat(const string & filename);
	void saveDot(const string & filename);

	void printGTTs() const;
    void print() const;

	/* Simulation */
	void generateRandInitTrajCol(TrajectoryCollection & trjcollection, int numberoftrj, int timesteps) const;

	//added by yangzhang 11/15/2008,simulate with external nodes enabled
	//modified by Yang Zhang 7/22/2012 to add new house noise model
	void simulate(const string & inittrjfile, const string & outputfile, int numberoftrj, int timesteps, double levelNoise=0, const string & houseModel="default") const;
	void simulate(TrajectoryCollection & trjCol) const;
	vector<int> simulate(const vector< vector<int> > & traj, int t) const; 

    //Haizhou Wang, Oct 11,2012
    /*
    Need the result of simulate() function within the program rather than save the result out to a file
    In the future, when rewriting this part, one should consdier let simulate() return a TraCol.
    Then, use "simulate(....).save(OUT_PUT_FILE_NAME) to save the result to an external file.

    Right now, I just simplely overlod simulate(), duplicaing its function body. The only difference is that this overloaded simulate() will return a TrajectoryCollection object.
    */
    const TrajectoryCollection simulate(const string & inittrjfile, int numberoftrj, int timesteps, double levelNoise=0, const string & houseModel="default") const;



    void simulateExhaustively(TrajectoryCollection & trjCol) const;
    
	list< vector<int> > computeStateTrans(); /** have it return a list**/

	/* Evaluation */
	void evaluate(int mode, string trajectoryfile, string outputfile, int offset, string differecefile, string summaryfile);

	/* set/get functions */

	void setMaxMarkovOrder(int order);
	int getMaxMarkovOrder() const;

    //Added by Haizhou Wang, Oct 25, 2012
    /*
    *NOTE: Neither m_min_Markov_order nor m_max_Markov_order are initialized in constructor or initializer.
    */
    void setMinMarkovOrder( int order){ m_min_Markov_order = order; }
    int getMinMarkovOrder() const { return m_min_Markov_order; }
    //End of Haizhou's Code


	void setVersion(const string & version) { m_version = version;}
	const string & getVersion() const { return m_version; }

	// const vector<Node> & getNodes() const { return m_nodes;} Removed MS Dec 4, 2011 

    const Node & getNode(int i) const { return m_nodes[i];} // Added MS Dec 4, 2011
    
    const vector<string> getNodeNames(char type) const;
    const vector<string> getNodeNames(const vector<int> & ids) const;
    
    const string & getNodeName(int i) const { return m_nodes[i].getName(); }
    
	//added by yangzhang 11.15.2008,return number of internal nodes to decide how to generate random initial values
	int getNumofInternalNode() const;
    
    void setNodeType(int i, char type) { m_nodes[i].setType(type); }
    char getNodeType(int i) const { return m_nodes[i].getType(); }
    
	// void setNodePvalue(int i, double pValue) { m_nodes[i].setPvalue(pValue); }

	//added by yangzhang 11/16/2008
	vector<int> getBases() const;

	void setGTTs(const vector<GeneralizedTruthTable> & gtts) { m_gtts = gtts; }
	void setGTT(int i, const GeneralizedTruthTable & gtt) { m_gtts[i] = gtt; }

	const vector<GeneralizedTruthTable> & getGTTs(void) const { return m_gtts; }
	const GeneralizedTruthTable & getGTT(int i) const { return m_gtts[i]; }

    size_t size() const { return m_nodes.size(); }

    //added by Yang Zhang 4.3.2009 to compare truth tables and parents in two gln files
	bool compareTruthTables(const GeneralizedLogicalNetwork & gln);
    
	void compareAllTruthTables(const GeneralizedLogicalNetwork & gln, 
                               const string & fileIntxCmp);

    // #########################################################################
	// ### SPECIFIC FOR THE ALCOHOL PROJECT WITH ELISSA CHESLER AND SUSAN 
    // BERGESON USED IN ERUASIP 2009 PAPER
	void simulateAcohol(const char * fileAcoNet, const char * fileDOT);
	void detectCycleAcohol(vector< vector<int> > & trajectory, const vector<int> & bases, list< vector< int > > & edgeList, vector<int> & vecStates, int prev_li);

	//Added by Haizhou Wang, March 23, 2012
	void eliminateFictitiousParents(const Topology & mustIncludeParents = Topology() );

    //Added by Haizhou Wang, Oct 25, 2012
    bool operator==( const GeneralizedLogicalNetwork & gln ) const;
    //End of Haizhou's code

    //Added by Haizhou Wang, March 18, 2013
    vector< std::pair<int,int> > getEdgeList() const;
    void saveEdgeList(const string & outputFileName) const;

};
