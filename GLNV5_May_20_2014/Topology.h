//
//  Topology.h
//  gln
//
//  Created by Joe Song on 10/6/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//

#ifndef gln_Topology_h
#define gln_Topology_h

#include <string>
using std::string;

#include <vector>
using std::vector;

class Topology {
    
public:
    
    Topology(): m_size(0), m_nConditions(0) {}
    
    Topology(const vector<string> & nodeNames, int nConditions) 
    {
        setNodeNames(nodeNames);
        m_nConditions = nConditions;
    }
    
    void set(int i, const vector<int> & parents) { m_topology[i][0] = parents; }
    
    void setNodeNames(const vector<string> & nodeNames) {
        m_nodeNames = nodeNames;
        m_size = nodeNames.size();
        m_topology = vector< vector< vector<int> > >( m_size, vector< vector<int> >(2) );                
    }
    
    void setNumConditions(int nConditions) { m_nConditions = nConditions; }
    
    void set(const vector<vector< vector<int> > > & topology)
    {
        m_topology = topology;
    }
    
    void setFile(const string & file) { m_file = file; }
	const string & getFile() const { return m_file; }

    const vector<vector< vector <int> > > & get() const { return m_topology; }
    
    vector<int> getAllNodes() const;
    vector<int> getAllChildNodes() const;
    
    size_t size() const { return m_size; }
    
    void scan(const string & file);
    void save(const string & file) const;

    void randomGenerate(vector< vector<int> > ParentsVector,
                        int NumberOfNodes,const char * fileTopology);

protected:
    
    // Changed from two dimension to three dimension, to include time delay 
    // for each parent by Yang Zhang 2.28.2011

    //childID: 0-based
    //m_topology[childID][0]: parent ID list, 1-based
    //m_topology[childID][1]: parent delay list.
    vector < vector < vector<int> > > m_topology;
    
    size_t m_size;
    
    vector<string> m_nodeNames;
    int m_nConditions;
    
    string m_file;

};

#endif
