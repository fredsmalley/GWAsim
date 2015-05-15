#pragma once

#include <string>
#include <set>
#include <utility>
#include <vector>

#include "GLN.h"
#include "EffCGLN_GlobalConstraint.h"
#include "DS_LocalConstraint.h"
#include "TrajectoryCollection.h"

class kStop_TestCase
{
public:

    static const DS_LocalConstraint       extractLocalConstraints ( const GeneralizedLogicalNetwork & gln, const double percentage=0.5);
    static const EffCGLN_GlobalConstraint extractGlobalConstraints( const TrajectoryCollection & traCol, const int numOfGlobalConstraints=5, const bool incompleteDynamicPattern=false, const bool sequential = false);

    static const std::vector< std::pair< std::pair<int,int>, unsigned> > RankEdges( const std::vector<GeneralizedLogicalNetwork> & GLNList );

    static const int findEdgeCount( const std::vector< std::pair< std::pair<int,int>, unsigned> > & edgeCounts, const std::pair<int,int> & e);

    static const double getGLNRankByEdgeMedian(const std::vector< std::pair< std::pair<int,int>, unsigned> > & edgeCounts, const GeneralizedLogicalNetwork & gln);

    static void outputGLNsRankedByEdgeMedian( const std::vector< std::pair< std::pair<int,int>, unsigned> > & edgeCounts, const std::vector<GeneralizedLogicalNetwork> & GLNList, const std::vector<std::string> & GLNNames);


    static const bool comp( const std::pair< std::pair<int,int>, unsigned> & op1, const std::pair< std::pair<int,int>, unsigned> & op2){return op1.second > op2.second;}

    static const bool comp2( const std::pair< int, double> & op1, const std::pair< int, double> & op2){return op1.second > op2.second;}

    static void test(const string & inputGLNFileName, const bool incompleteDynamicPattern=false, const double percentage=0.5, const int numOfGlobalConstraints=5, const std::string searchDirection="random", const int depth=5 );

    static void output_PR_ROC( const unsigned networkSize, const std::set< std::pair<int,int> > & originalEdgeList, const std::vector< std::pair<int, int> > & predictedEdgeList, const string & outputFileHead="");
};
