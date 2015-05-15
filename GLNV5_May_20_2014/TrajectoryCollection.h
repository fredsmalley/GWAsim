// TrajectoryCollection.h
//
// Joe Song
// Created: September 27, 2008
// Last modified: November 9, 2008.  Joe Song.  Added missing value handling.

#pragma once

#include <string>
using std::string;

#include <list>
using std::list;

#include <vector>
using std::vector;

#include <fstream>
using std::ifstream;

class TrajectoryCollection
{

public:
    
    TrajectoryCollection(void) : m_hasMissingValues(true) {}

	TrajectoryCollection(const vector<string> & listTrajFiles);

	void allocate(int nTrajs, int lenEachTraj, int nVariables, int initValue);

    // I/O member functions
	size_t scan(const vector<string> & listTrajFiles);
	bool scan(const string & fileTraj);

	void print() const;
	
	void save(const char* path) const; //added by yangzhang 11.12.2008, save trajectory collection file

	void setVersion(const string & version) { m_version = version; }
	const string & getVersion() const { return m_version; }

	int get(int trajID, int t, int i) const { return m_trajTables[trajID][t][i]; }
    
	void set(int trajID, int t, int i, int value) { 
        m_trajTables[trajID][t][i] = value; 
        if(m_missingValue == value) {
            m_hasMissingValues = true;   
        }    
    }

	bool missing(int trajID, int t, int nodeID) const { return m_missingValue == m_trajTables[trajID][t][nodeID]; }

    bool hasMissingValues() const { return m_hasMissingValues; }
    
    bool checkMissingValues();
    
	size_t getNumTrajectories() const { return m_trajTables.size(); }
	
	// int getTrajectoryLength(int trajID) const { return m_trajLengths[ trajID ]; }
    int getTrajectoryLength(int trajID) const {
        int length = 0;
        if(trajID < (int) m_trajTables.size()) {
            length = m_trajTables[ trajID ].size();
        }
        return length;
    }
    
	// void setTrajectoryLength(vector<int> trjlength) { m_trajLengths = trjlength;}

	// Modified from const to non-const by yangzhang 11.16.2008
	// Modified back to const again. Joe Song. November 28, 2008 
	const vector< vector< vector<int> > > & getTrajTables() const { return m_trajTables; } 

	//added by yangzhang 11.15.2008
	void setTrajTables(const vector< vector< vector<int> > > & trjTables);

	const vector< vector< int > > & getTrajectory(int i) const { return m_trajTables[i]; } 
	void setTrajectory(int i, const vector< vector< int > > & trj) { m_trajTables[i] = trj; }

	void setState(int i, int t, const vector<int> & state) { m_trajTables[i][t] = state; }
	const vector<int> & getState(int i, int t) const { return m_trajTables[i][t]; }

	int base(int i) const { return m_bases[i]; }
	
	void setBases(const vector<int> & bases) { m_bases = bases; }
	const vector<int> & getBases() const { return m_bases; }

	int nNodes() const { return m_nNodes; }
	//added by yangzhang 11/16/2008
	void setn_nodes(int n_nodes) { m_nNodes = n_nodes;}	

	const vector<int> & getExtNodes() const { return m_extNodes; }
	const vector<string> & getIntNodeNames() const { return m_intNodeNames; }
	const vector<string> & getExtNodeNames() const { return m_extNodeNames; }

    const vector<string> getNodeNames(const vector<int> & ids) const;
    
	//added by yangzhang 11/16/2008
	void setIntNodeNames(const vector<string> & intnames) { m_intNodeNames = intnames;}
	void setExtNodeNames(const vector<string> & extnames) { m_extNodeNames = extnames;}

    
	//added by Yang Zhang 4/24/2009
	const vector<bool> & getBeParent() const {return m_beParent; }
	const vector<bool> & getHaveParent() const {return m_haveParent; }

	void setBeParent(const vector<bool> & beParent) { m_beParent =  beParent; }
	void setHaveParent(const vector<bool> & haveParent) { m_haveParent =  haveParent; }

	//added by yangzhang 2/8/2009 add noise to one node
	//modified to house noise model 7/21/2012
	void addNoiseToNode(int base,double noise,int nodeid,const string & HouseNoiseModel);
	
	//add noise to whole trj
	//modified by Yang Zhang 7/22/2012 to add option to choose new house noise model
	void addNoise(double noiselevel, const string & housemodel);

	//added by YangZhang 8/21/2009
	const vector<vector<int> > & getexcludeTrajectory() const {return m_excludeTrajectory; }
	void setexcludeTrajectory(const vector<vector<int> > & excludeTrajectory) 
	{ m_excludeTrajectory = excludeTrajectory; }

	//added by Yang Zhang 7.30.2010
	void permuteTrjtableSwapingColumn();
	void permuteTrjtableOnTime();
	void permuteTrjtableOnColumn();


    //Added by Haizhou Wang, July 24, 2012
    bool CheckTraDeterministic() const;

private:

	static const int m_missingValue = -1;

    bool m_hasMissingValues;
    
	string m_version;

	vector<string> m_files;

	/*added by YangZhang 8/21/2009, a matrix indicate 
	which trajectory can be used to accumulate transition table 
	for parent selecting.
	The idea is coming from the need for Dream 4 challenge,
	in a trajectory which one node is knock-out/knock-down or inhibited,
	this trajectory can not be used to select parents for this node anymore.
	*/
	vector<vector<int> > m_excludeTrajectory;

	vector< vector< vector<int> > > m_trajTables;		//The list holds trajectory tables
	// vector<int> m_trajLengths;	// Removed MS 4/28/2013	//Holds the number of rows for each trajectory table

	int m_nNodes;		// Number of nodes in the network

	vector<int> m_bases;		//Will be pointing to array that holds the bases (radix)

	vector<int> m_extNodes;
	vector<string> m_intNodeNames;
	vector<string> m_extNodeNames;

	//added by Yang Zhang
	vector<bool> m_beParent; 
	vector<bool> m_haveParent; 

	//added by Yang Zhang
	vector< vector<int> > readTraj(ifstream &in, int iteration,int &columns, int &rows, int &num_of_traj, vector <bool> &beParent,vector <bool> &haveParent);

	vector< vector<int> > readTraj(ifstream &in, int iteration,int &columns, int &rows, int &num_of_traj);
};
