#pragma once

#include <deque>
#include <set>
#include <vector>

#include "DS_LocalConstraint.h"

enum Color{WHITE, BLACK, GREY};


const std::vector<int> getBits(int number, const std::vector<int> & bitRadixs);
const vector<int> nextVertexContainsPattern(const std::vector<int> & Pattern, const vector<int> & unfixedRadix, int k);
const std::deque<int> getVertexContainPattern(const std::vector<int> & Pattern, const std::vector<int> & radixs);

template<typename T>
const bool sameSign( const vector<T> & v );


class kStop_Node
{
public:
    kStop_Node():m_ID(-1),m_color(WHITE),m_mark(-1),m_ParentID(-1),m_depthLooked(0), m_adjList(){}

    kStop_Node(const kStop_Node & rhs);
    kStop_Node & operator=(const kStop_Node & rhs);

    ~kStop_Node(){};

    void buildAdjList( const DS_LocalConstraint & lc );
    void clearAdjList(){ m_adjList.clear(); }
    void showAdj() const;




/*
In the view of super state transition diagram,
a state in the diagram is a list of variables with each variable corresponds to a node in the topology graph.

Thus, the ID of the state is the list of variables in decimal format.

For example, if a system has 3 binary nodes and a state from the super state transition diagram is (0,1,0)
Then the ID for this state is 2.

If the state is (1,0,1), the ID of this state is 5.
*/
public:
    int m_ID;    //ID of state in super state transition diagram
    Color m_color;
    int m_mark;  //mark represents this node is used by which segment
    int m_ParentID;
    int m_depthLooked;

    std::deque<int> m_adjList; //Records the ID of adjacent nodes
};
