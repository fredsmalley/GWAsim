//
//  PathwayStats.h
//  gln
//
//  Created by Joe Song on 3/11/13.
//
//

#ifndef __gln__PathwayStats__
#define __gln__PathwayStats__

#include <iostream>

#include "ChisqTests.h"
#include "TransitionTable.h"
#include "Topology.h"

class SinglePathwayStats {
public:
    Chisq   m_tot;
    void append(const string & file) const;
    void save(const string & file) const;
};

class MultiplePathwayStats {

public:

    friend bool operator == (const MultiplePathwayStats & s1, const MultiplePathwayStats & s2);
    friend bool operator != (const MultiplePathwayStats & s1, const MultiplePathwayStats & s2);
    
    MultiplePathwayStats() {}
    
    const MultiplePathwayStats & approximate(const vector<MultiplePathwayStats> & statsBS);
    void compare(const Topology & topo,
                 const vector<TransitionTable> & tts,
                 const vector<TransitionTable> & ttsPooled,
                 const vector<double> & chisqs_t,
                 const vector<int> & dfs_t,
                 const vector<double> & chisqs_z,
                 const vector<int> & dfs_z,
                 const vector<double> & child_chisqs_z,
                 const vector<int> & child_dfs_z);
    
    void scale(const vector<TrajectoryCollection> & trjCols,
               const vector<int> & nodeIDsOnPathway);
    
    void save(const string & file) const;
    
public:
    Chisq m_het;    // pathway interaction heterogeneity
    Chisq m_hom;    // pathway interaction homogeneity
    Chisq m_tot;    // pathway total interaction strengh
    Chisq m_wzChildren;  // pathway children working zone
    Chisq m_wzParents;   // pathway parent working zone
};


#endif /* defined(__gln__PathwayStats__) */
