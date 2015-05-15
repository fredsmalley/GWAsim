#include "kStop.h"

#include <cassert>
#include <cstdlib>

#include <algorithm>
using std::random_shuffle;
using std::equal;
using std::make_pair;

#include <ctime>

#include <deque>
using std::deque;

#include <numeric>    //In order to use std::accumulate()
using std::accumulate;

#include <functional> //In order to use std::multiplies<>()

#include <iostream>
using std::cerr;
using std::endl;
using std::cout;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <stack>
using std::stack;

#include <set>
using std::set;



#include "Transition.h"   //In order to use inline function digitToNumber()


//As of Oct 22, 2012, delegating consturctor is not widely supported.
//The following two constructor have duplicate codes.
kStop::kStop(const string LocalConstraintFileName, const string GlobalConstraintFileName, const string SearchDirection, const int depth):
         m_graph(), m_MustInclude(), m_FinalPath(), m_segChangeable(),
         m_lc(LocalConstraintFileName), m_gc(GlobalConstraintFileName),
         m_SearchDirection(SearchDirection), m_SearchDepth(depth)
{

    if(  m_lc.m_nodes!= m_gc.m_nodes  ||  m_lc.m_radixs != m_gc.m_radixs   )
    {
        cerr<<"kStop::kStop(): Error!\nInformation in local constraint file does not match global constraint file.\n"<<endl;
        exit(-1);
    }

    for( const auto & DP: m_gc.m_DynamicPatternList)
    {
        if( DP.m_DynamicPattern.size() == 1 )
        {
            cerr<<"In one of the global constraints, one may have only one dynamic pattern. Applications on only one dynamic pattern is meaningless\n"<<endl;
            exit(-1);
        }
    }

    //Calculate how many states are there
    int numberOfVertices = accumulate( m_lc.m_radixs.begin(), m_lc.m_radixs.end(), 1, std::multiplies<int>() );

    //Initialize node IDs
    m_graph.resize( numberOfVertices );
    for(unsigned t=0; t<m_graph.size(); ++t)
    {
        m_graph[t].m_ID = t;
        m_graph[t].m_color = WHITE;
        m_graph[t].m_mark = -1;
        m_graph[t].m_ParentID = -1;
    }

    m_segChangeable = vector<bool>(m_gc.m_DynamicPatternList[0].m_DynamicPattern.size()-1, true);
    assert( m_segChangeable.size() >= 1);
}



kStop::kStop(const DS_LocalConstraint & lc, const EffCGLN_GlobalConstraint & gc, const std::string SearchDirection, const int depth)
        :m_graph(), m_MustInclude(), m_FinalPath(), m_segChangeable(),
         m_lc(lc), m_gc(gc),
         m_SearchDirection(SearchDirection), m_SearchDepth(depth)

{

    if(  m_lc.m_nodes!=m_gc.m_nodes  ||  m_lc.m_radixs!=m_gc.m_radixs   )
    {
        cerr<<"kStop::kStop(): Error!\nInformation in local constraint file does not match global constraint file.\n"<<endl;
        exit(-1);
    }

    for( const auto & DP: m_gc.m_DynamicPatternList)
    {
        if( DP.m_DynamicPattern.size() == 1 )
        {
            cerr<<"In one of the global constraints, one may have only one dynamic pattern. Applications on only one dynamic pattern is meaningless\n"<<endl;
            exit(-1);
        }
    }

    //Calculate how many states are there
    int numberOfVertices = accumulate( m_lc.m_radixs.begin(), m_lc.m_radixs.end(), 1, std::multiplies<int>() );

    m_graph.resize( numberOfVertices );
    for(unsigned t=0; t<m_graph.size(); ++t)
    {
        m_graph[t].m_ID = t;
        m_graph[t].m_color = WHITE;
        m_graph[t].m_mark = -1;
        m_graph[t].m_ParentID = -1;
    }

    m_segChangeable = vector<bool>(m_gc.m_DynamicPatternList[0].m_DynamicPattern.size()-1, true);
    assert( m_segChangeable.size() >= 1);
}



void kStop::reset()
{
    m_MustInclude.clear();
    m_FinalPath.clear();
    m_segChangeable = vector<bool>(m_gc.m_DynamicPatternList[0].m_DynamicPattern.size()-1, true);
    for( auto & node: m_graph)
    {
        node.m_color = WHITE;
        node.m_mark = -1;
        node.m_ParentID = -1;
        node.m_depthLooked = 0;
    }
}


//Need to make sure all member data are initialized
//For example, FinalPath must be empty everytime FindPath() is called, otherwise the reuslt is  mysterious.
bool kStop::FindPath()
{
    //Get a set of starting states if the first dynamic pattern is incomplete.
    deque<int> startingSet( getVertexContainPattern( m_gc.m_DynamicPatternList[0].m_DynamicPattern.front(), m_gc.m_radixs) );
    random_shuffle( startingSet.begin(), startingSet.end(), RNG);

    //Try enumerations of the first dynamic pattern.
    //Once a solution is found, stop.
    for( const auto & t: startingSet)
    {
        m_MustInclude.push_back( t );

        if( initialPathFinding() && detour() )
            return true;

        reset();  //Reset kStop object
    }

    //After tried all possible enumeration of the fisrt dynamic pattern, still no solution found.
    return false;
}



bool kStop::initialPathFinding()
{
    size_t seg=0;
    while( seg < m_segChangeable.size() )
    {
        set<int> forbidSet( m_MustInclude.begin(), m_MustInclude.end() ); 
        forbidSet.erase( m_MustInclude[seg] );

        set<int> usedVertexSet;
        for(const auto & it: m_FinalPath )
        {
            usedVertexSet.insert( (it).m_SegmentPath.begin(), (it).m_SegmentPath.end() );
        }

        set<int> avoidSet(forbidSet);
        avoidSet.insert( usedVertexSet.begin(), usedVertexSet.end() );
        if( avoidSet.find( m_MustInclude[seg] ) != avoidSet.end() )
            avoidSet.erase( m_MustInclude[seg] );

        SegmentPath segPath( SIDDFS(avoidSet, seg, m_MustInclude[seg], m_gc.m_DynamicPatternList[0].m_DynamicPattern[seg+1], m_SearchDirection, m_SearchDepth)    ) ;

        if( segPath.empty() )
        {
            if( forbidSet != avoidSet)//if forbidSet == avoidSet, no need to apply SIDDFS again.
            {
                //Use more loose constraints, search again
                segPath = SIDDFS(forbidSet, seg, m_MustInclude[seg], m_gc.m_DynamicPatternList[0].m_DynamicPattern[seg+1], m_SearchDirection, m_SearchDepth);

                if( segPath.empty() )
                {   
                    return false;
                }
                else
                {
                    m_FinalPath.push_back(segPath);
                    m_MustInclude.push_back(segPath.m_SegmentPath.front());
                }
            }
            else
            {
                return false;
            }
        }
        else
        {//This segment path may be changed in detour process
            m_FinalPath.push_back(segPath);
            m_MustInclude.push_back(segPath.m_SegmentPath.front());
        }

        ++seg; //Move to find simple path for next segment
    }//end of while()


    assert( m_MustInclude.size() == m_gc.m_DynamicPatternList[0].m_DynamicPattern.size() );
    return true;
}



bool kStop::detour()
{
    for(int t=0; t<(int)m_FinalPath.size(); ++t)
    {
        if( !m_FinalPath[t].pure(m_graph,t) )
        {
            if( !detourSegment(t) )
                return false;
            t = -1; 
            //After detour, previous segment may also being changed. Thus, need to check from the beginning again. This might be very time consuming.
        }
    }
    return true;
}




/*
* Given a segment to detour
* 1. Try to find a detour without using any vertex being used or any vertex in the must-include set
* 2. If no solution, then try to find a detour by considering using vertex on other segment paths; but, vertex on UNCHANGEABLE segment path or in must-include set can NOT be used.
* 3. If still no solution, stop and report failure.
*/
bool kStop::detourSegment(unsigned segID)
{
    //Detour will take the opposite order to traverse adjacent list
    string searchDirection;
    if( m_SearchDirection == "forward")
        searchDirection = "backward";
    else if( m_SearchDirection == "backward")
        searchDirection = "forward";
    else if( m_SearchDirection == "random")
        searchDirection = m_SearchDirection;
    else
    {
        cerr<<"Error: Unknown search order\n"<<endl;
        exit(-1);
    }

    set<int> forbidSet(m_MustInclude.begin(), m_MustInclude.end() );
    forbidSet.erase( m_MustInclude[segID] );
    forbidSet.erase( m_MustInclude[segID+1] );
    
    set<int> usedVertexSet;
    for(size_t t=0; t<m_FinalPath.size(); ++t)
    {
        if( t == segID )
        {
            for( const auto & v: m_FinalPath[t].m_SegmentPath )
            {
                if( m_graph[ v ].m_mark != (int)t )
                    usedVertexSet.insert( v );
            }
        }
        else
            usedVertexSet.insert( m_FinalPath[t].m_SegmentPath.begin(), m_FinalPath[t].m_SegmentPath.end() );
    }
    if( usedVertexSet.find(m_MustInclude[segID]) != usedVertexSet.end() )
        usedVertexSet.erase( m_MustInclude[segID] );
    if( usedVertexSet.find(m_MustInclude[segID+1]) != usedVertexSet.end() )
        usedVertexSet.erase( m_MustInclude[segID+1]);
    

    set<int> avoidSet(forbidSet);
    avoidSet.insert(usedVertexSet.begin(), usedVertexSet.end());
    assert( avoidSet.find( m_MustInclude[segID] ) == avoidSet.end() );
    assert( avoidSet.find( m_MustInclude[segID+1] ) == avoidSet.end() );

    SegmentPath segPath( SIDDFS(avoidSet, segID, m_MustInclude[segID], getBits(m_MustInclude[segID+1],m_gc.m_radixs), searchDirection, m_SearchDepth) ) ;

    if( segPath.empty() ) //Can not find a detour with tight constraint. Loose the constraint, search again.
    {
        //Vertices on unchangeable path can NOT be used.
        set<int> mustAvoidVertexSet;
        for(size_t t=0; t<m_segChangeable.size(); ++t)
        {
            if( !m_segChangeable[t] )
            {
                assert( m_FinalPath[t].pure(m_graph,t) );
                mustAvoidVertexSet.insert( m_FinalPath[t].m_SegmentPath.begin(), m_FinalPath[t].m_SegmentPath.end() );
            }
        }
        if( mustAvoidVertexSet.find(m_MustInclude[segID]) != mustAvoidVertexSet.end() )
		    mustAvoidVertexSet.erase( m_MustInclude[segID] );
        if( mustAvoidVertexSet.find(m_MustInclude[segID+1]) != mustAvoidVertexSet.end() )
		    mustAvoidVertexSet.erase( m_MustInclude[segID+1] );
		
        avoidSet.clear();
        avoidSet = forbidSet;
        avoidSet.insert( mustAvoidVertexSet.begin(), mustAvoidVertexSet.end() );
        assert( avoidSet.find( m_MustInclude[segID] ) == avoidSet.end() );
        assert( avoidSet.find( m_MustInclude[segID+1] ) == avoidSet.end() );


        segPath = SIDDFS(avoidSet, segID, m_MustInclude[segID], getBits(m_MustInclude[segID+1],m_gc.m_radixs), searchDirection, m_SearchDepth);
        
        if( segPath.empty() ) //Can NOT find a detour even with loosed constraint. Report failure, quit program.
        {
            return false;
        }
        else
        {
            m_FinalPath[segID] = segPath;
            m_segChangeable[segID] = false;
        }
    }
    else
    {
        m_FinalPath[segID] = segPath;
    }

    return true;
}





/*
* A state goes back to itself is not allowed.
* Assume a state could go back to itself, this may lead to a state matching two dynamic patterns.
*/
const SegmentPath kStop::DFS_Visit(const set<int> & excludeSet, int segmentID, int startID, const std::vector<int> & targetPattern, const string & searchDirection, const int depth)
{
    if( depth < m_graph[startID].m_depthLooked )
       return SegmentPath();

    assert( excludeSet.find(startID) == excludeSet.end() );

    m_graph[startID].m_color = GREY;

    if( m_graph[startID].m_adjList.empty() )
        m_graph[startID].buildAdjList( m_lc );
 
    //If target node is adjacent to 'this' node, return immediately
    deque<int> tmp(m_graph[startID].m_adjList);
    random_shuffle( tmp.begin(), tmp.end(), RNG);

    for( const auto & t: tmp )
    {
        if( m_graph[t].m_color == WHITE && 
            matchPattern(t,targetPattern) && 
            excludeSet.find(t) == excludeSet.end() 
          )
        {
            m_graph[t].m_mark = segmentID;
            m_graph[t].m_ParentID = startID;

            SegmentPath seg;
            seg.push_back(t);

            m_graph[startID].m_color = WHITE;
            return seg;
        }
    }
    if( depth == 1)
    {//Target node can not be reached within 1 step and search depth is met.
        m_graph[startID].m_color = WHITE;
        return SegmentPath();
    }
    


    //Decide how to traverse the adjacent list: forward, backward or random 
    tmp.clear();
    tmp.assign(m_graph[startID].m_adjList.begin(), m_graph[startID].m_adjList.end() );
    if( searchDirection == "forward")
        ;
    else if( searchDirection == "backward")
        reverse( tmp.begin(), tmp.end() );
    else if( searchDirection == "random")
        random_shuffle( tmp.begin(), tmp.end(), RNG);
    else
    {
        cerr<<"Error: Unknown search order in DFS_Visit()\n"<<endl;
        exit(-1);
    }



    //Traverse nodes in adjacent list
    for( const auto & NodeID: tmp )
    {    
        if( excludeSet.find( NodeID ) == excludeSet.end() && m_graph[NodeID].m_color ==  WHITE)
        {
            m_graph[NodeID].m_ParentID = startID;

            SegmentPath result( DFS_Visit(excludeSet, segmentID, NodeID, targetPattern, searchDirection, depth-1) );
            if( !result.empty() )
            {
                m_graph[NodeID].m_mark = segmentID;
                result.push_back( NodeID );

                m_graph[startID].m_color = WHITE;
                return result;
            }
        }
    }
    m_graph[startID].m_color = WHITE;

    //If not returned from previous steps, it means that with the 'depth', no result could be found.
    //This node will be explored again next time only when the given 'depth' is greater than m_depthLooked
    m_graph[startID].m_depthLooked = depth;



    //When search detph is huge, memory will become an issue
    //Then program will try to use adjacent lists on-the-fly.
    if( m_graph.size() > 100000 && m_SearchDepth > 10 )
        m_graph[startID].clearAdjList(); 

    return SegmentPath();
}






const SegmentPath kStop::SIDDFS(const std::set<int> & excludeSet, const int segID, const int startID, const std::vector<int> & targetPattern, const std::string & searchDirection, const int depth)
{
    //No need to start with searchDetph=0, since state 'startID' can NOT match
    //the starting and ending dynamic patterns at the same time.
    int iteration = 1;
    SegmentPath result;
    
    while( iteration < depth)
    {
        for(auto & node: m_graph)
            node.m_depthLooked = 0;

        m_graph[startID].m_ParentID = startID;

        assert( result.empty() );
        result = DFS_Visit(excludeSet, segID, startID, targetPattern, searchDirection, iteration);

        if( !result.empty() )
        {
            m_graph[startID].m_mark = segID;
            result.push_back( startID );
            return result;
        }

        ++iteration;
    }


    for(auto & node: m_graph)
        node.m_depthLooked = 0;
    
    
    //Do it one last time when search depth equals to 'depth'
    m_graph[startID].m_ParentID = startID;
    assert( result.empty() );
    result = DFS_Visit(excludeSet, segID, startID, targetPattern, searchDirection, depth);

    if( !result.empty() )
    {
        m_graph[startID].m_mark = segID;
        result.push_back( startID );
        return result;
    }

    return SegmentPath();
}




//Given a vertex ID (as decimal)
//First, convert the ID into bit form
//Then, compare to 'pattern' to see if it is a match
const bool kStop::matchPattern(const int ID, const vector<int> & pattern)
{
    vector<int> bitForm( getBits(ID, m_gc.m_radixs) );
    assert( bitForm.size() == pattern.size() );

    for(size_t t=0; t<bitForm.size(); ++t)
    {
        if( pattern[t] != -1 && pattern[t] != bitForm[t])
                return false;
    }
    return true;
}

//Output final path on standard output
void kStop::showFinalPath() const 
{
    if( m_FinalPath.size() == 0)
    {
        cout<<"Final path is empty. No path to show"<<endl;
    }
    else
    {
        for(unsigned t=0; t<m_FinalPath.size(); t++)
        {
            m_FinalPath[t].showSegmentpath();
        }
    }

}

//Save final path as Trajectory file format version 2.1
void kStop::saveTra(const string outputfile) const 
{
    if( m_FinalPath.empty() )
    {
        cout<<"No result to be saved\n"<<endl;
        return;
    }

    ofstream ou(outputfile.c_str());

    ou<<"TRAJECTORY_VER2.1"<<endl;
    ou<<"1\t"<<m_lc.m_radixs.size()<<"\t"<<"0"<<endl;

    //Output radix of each node
    for(size_t i=0; i<m_lc.m_radixs.size(); i++)
        ou<<m_lc.m_radixs[i]<<"\t";
    ou<<endl;

    //Output name of each node
    for(size_t i=0; i<m_lc.m_nodes.size(); i++)
        ou<<m_lc.m_nodes[i]<<"\t";
    ou<<endl<<endl;

    ou<<endl; //Line of names of external nodes is empty


    //Ouput if a node can be parent
    for(const auto & i: m_lc.m_nodes )
    {
        if( find( m_lc.m_ExcludedParent.cbegin(), m_lc.m_ExcludedParent.cend(), i) != m_lc.m_ExcludedParent.cend() )
        {
            ou<<"No\t";
        }
        else
            ou<<"Yes\t";
    }
    ou<<endl;


    //Output if a node can have parent(s)
    for(const auto & i: m_lc.m_nodes )
    {
        if( find( m_lc.m_ExcludedChild.cbegin(), m_lc.m_ExcludedChild.cend(), i) != m_lc.m_ExcludedChild.cend() )
        {
            ou<<"No\t";
        }
        else
            ou<<"Yes\t";
    }
    ou<<endl<<endl;



    //Calculate how many rows in this trajectory
    int totalRows = 0;
    for(size_t i=0; i<m_FinalPath.size(); i++)
        totalRows += m_FinalPath[i].m_SegmentPath.size();
    totalRows = totalRows - (int)m_FinalPath.size() + 1;
    ou<<totalRows<<endl;


    //Output the actual trajectory
    for(size_t i=0; i<m_FinalPath.size();i++)
    {
        for(size_t j=m_FinalPath[i].m_SegmentPath.size()-1; j>0; j--)
        {
            vector<int> bitForm( getBits(m_FinalPath[i].m_SegmentPath[j], m_lc.m_radixs) );
            for(size_t k=0; k<bitForm.size(); k++)
                ou<<bitForm[k]<<"\t";
            ou<<endl;
        }
    }

    //Output last must-include vertex
    vector<int> bitForm( getBits(m_FinalPath.back().m_SegmentPath.front(), m_lc.m_radixs) );
    for(size_t k=0; k<bitForm.size(); k++)
        ou<<bitForm[k]<<"\t";
    ou<<endl;


    ou.close();
}



const TrajectoryCollection kStop::getTrajectoryCollection() const
{
    TrajectoryCollection result;

    result.setVersion("TRAJECTORY_VER2.1");
    result.setn_nodes( m_gc.m_nodes.size() );
    result.setBases( m_gc.m_radixs );
    result.setIntNodeNames( m_gc.m_nodes );
    result.setExtNodeNames( vector<string>() );
    
    vector<bool> beParent;
    for(const auto & i: m_lc.m_nodes )
    {
        if( find( m_lc.m_ExcludedParent.cbegin(), m_lc.m_ExcludedParent.cend(), i) != m_lc.m_ExcludedParent.cend() )
        {
            beParent.push_back(false);
        }
        else
            beParent.push_back(true);
    }
    result.setBeParent( beParent );
    
    
    vector<bool> haveParent;
    for(const auto & i: m_lc.m_nodes )
    {
        if( find( m_lc.m_ExcludedChild.cbegin(), m_lc.m_ExcludedChild.cend(), i) != m_lc.m_ExcludedChild.cend() )
        {
            haveParent.push_back(false);
        }
        else
            haveParent.push_back(true);
    }
    result.setHaveParent( haveParent );


    vector< vector<int> > tra;

    for(size_t i=0; i<m_FinalPath.size();i++)
    {
        for(size_t j=m_FinalPath[i].m_SegmentPath.size()-1; j>0; j--)
        {
            vector<int> bitForm( getBits(m_FinalPath[i].m_SegmentPath[j], m_lc.m_radixs) );
            tra.push_back( bitForm );
        }
    }
    vector<int> bitForm( getBits(m_FinalPath.back().m_SegmentPath.front(), m_lc.m_radixs) );
    tra.push_back( bitForm );

    result.setexcludeTrajectory( vector< vector<int> >() );
    result.setTrajTables(  vector< vector< vector<int> > >(1, tra) );

    return result;
}


#include <map>
#include <utility>
using std::map;
using std::pair;


bool kStop::checkMust_Include_Set_Pure() const
{
    map<int, int> histo;
    for( const auto & t: m_MustInclude )
    {
        if( histo.find( t ) == histo.end() )
            histo.insert( make_pair(t, 1) );
        else
            ++histo[t];
    }

    for( const auto & p: histo)
        if( p.second > 1 )
            return false;


    return true;
}





const GLNConfig kStop::getGLNConfig(const string & fileName, string t) const
{
    GLNConfig gc;
    gc.m_alpha = 1;
    gc.m_allowSelfCycle = true;

    vector<string> tmp;
    tmp.push_back( fileName + "-" + t + ".tra" );
    gc.m_listTrajFiles = tmp ;
    
    TrajectoryCollection trajCol;
    trajCol.scan( fileName + "-" + t + ".tra");
    
    gc.m_mustIncludeParents.setNodeNames( trajCol.getIntNodeNames() );
    gc.m_mustIncludeParents.scan(fileName + "LocalConstraintTopology.txt"); 
    gc.m_mustIncludeParents.setFile(fileName + "LocalConstraintTopology.txt");

    gc.m_BP = trajCol.getBeParent();
    gc.m_HP = trajCol.getHaveParent();

    return gc;
}
