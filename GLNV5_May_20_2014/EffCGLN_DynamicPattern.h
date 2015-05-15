#pragma once

#include <vector>

class EffCGLN_DynamicPattern
{
public:
    EffCGLN_DynamicPattern();
    ~EffCGLN_DynamicPattern();

public:
    //It's better to use:
    //std::vector< vector< pair<NodeID, NodeValue> > >
    std::vector< std::vector<int> > m_DynamicPattern;
};