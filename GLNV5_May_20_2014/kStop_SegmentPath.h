#pragma once

#include <vector>

#include "kStop_Node.h"



class SegmentPath
{
public:
    SegmentPath():m_SegmentPath(){}
    SegmentPath(int size): m_SegmentPath(size){}

    ~SegmentPath(){}

public:
    bool empty(){ return m_SegmentPath.empty(); }
    bool pure( const std::vector<kStop_Node> & G, int segID ) const;
    void showSegmentpath() const;
    void push_back(int ID){ m_SegmentPath.push_back(ID); }
    void clear(){ m_SegmentPath.clear(); }

public:
    std::vector<int> m_SegmentPath;  //A sequence of IDs of states(vertices) in super state transition diagram 
};
