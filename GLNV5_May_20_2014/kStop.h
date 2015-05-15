#pragma once

#include <cstddef>
#include <ctime>
#include <fstream>
#include <string>
#include <set>
#include <stack>
#include <vector>

#include "DS_LocalConstraint.h"
#include "EffCGLN_GlobalConstraint.h"
#include "GLNConfig.h"
#include "kStop_Node.h"
#include "kStop_SegmentPath.h"
#include "TrajectoryCollection.h"


//Random number generator to be used in random_shuffle()
inline ptrdiff_t RNG( ptrdiff_t i){return rand()%i;}



class kStop
{
public:
    kStop(const std::vector<kStop_Node> & graph, const std::vector<int> & mustInclude):m_graph(graph),
        m_MustInclude(mustInclude),
        m_FinalPath(),
        m_segChangeable(m_MustInclude.size()-1, true),
        m_lc(),m_gc(),
        m_SearchDirection("forward"),
        m_SearchDepth(5)
    {}


    kStop(const DS_LocalConstraint & lc, const EffCGLN_GlobalConstraint & gc, const std::string SearchDirection, const int depth);

    //As of Oct 22, 2012, delegating consturctor is not widely supported.
    //kStop(const std::string LocalConstraintFileName, const std::string GlobalConstraintFileName, const std::string SearchDirection, const int depth): kStop( DS_LocalConstraint( LocalConstraintFileName), EffCGLN_GlobalConstraint( GlobalConstraintFileName), SearchDirection, depth) {}
    
    kStop(const std::string LocalConstraintFileName, const std::string GlobalConstraintFileName, const std::string SearchDirection, const int depth);

    kStop(){}
    ~kStop(){}

    void showFinalPath() const;
    void saveTra( const std::string outputfile ) const;
    bool checkMust_Include_Set_Pure() const;
    bool FindPath();
    const GLNConfig getGLNConfig(const string &, string) const;
    const TrajectoryCollection getTrajectoryCollection() const;



private:
    std::vector<kStop_Node> m_graph;
    std::vector<int> m_MustInclude; 		//M is the set of k stops of vertices/nodes
    std::vector<SegmentPath> m_FinalPath; 	//Simple path for each segment
    std::vector<bool> m_segChangeable;

    DS_LocalConstraint m_lc;
    EffCGLN_GlobalConstraint m_gc;

    std::string m_SearchDirection;
    int m_SearchDepth;

private:
    bool initialPathFinding();
    bool detour();
    bool detourSegment( unsigned segID);  //Given a segment ID, find a detour

    const SegmentPath SIDDFS(const std::set<int> & excludeSet, const int segmentID, const int startID, const std::vector<int> & targetPattern, const std::string & searchDirection, const int depth);

    const SegmentPath DFS_Visit(const std::set<int> & excludeSet, int segmentID, int startID, const std::vector<int> & targetPattern, const std::string & searchDirection, const int depth);

    const bool matchPattern( const int state, const std::vector<int> & pattern);

    void reset();

    //Copy constructor and assignment operator are disabled.
    kStop(const kStop &);
    kStop & operator=(const kStop &);
};
