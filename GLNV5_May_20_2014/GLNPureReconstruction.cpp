#include "GLNGlobals.h"

#include <cassert>
#include <algorithm>

#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;

#include <numeric>

#include <sstream>
using std::stringstream;

#include <set>
using std::set;

#include "EnvGLNRec.h"

#include "Transition.h"
#include "TransitionTable.h"
#include "TrajectoryCollection.h"


/*
Given a row from trajectory collection, a vector containing parents' index and another vector containing parents' corresponding radix,
this function returns the decimal form of parents.
*/
int getParentCombinationIndex(const vector<int> & v, const vector<int> & p, const vector<int> & base)
{
    vector<int> bits2Trans(p.size());

    for(size_t i=0; i<p.size(); i++)
    {
        assert( p[i] <= (int)v.size() );
        bits2Trans[i] = v[p[i]];
    }

    //return ditgits2Num(bits2Trans, base);
    return digitsToNumber(bits2Trans, base);
}

/*
GeneralizedLogicalNetwork ReconstructPureGLN( const TrajectoryCollection & traCol, const bool & selfControl)
{
    ;
}

GeneralizedLogicalNetwork ReconstructPureGLN( const TrajectoryCollection & traCol, const vector<int> & radixs, const vector<string> & nodeNames, const vector<bool> & beParent, const vector<bool> & haveParent, const vector< vector<int> > & excludeLists, const bool & selfControl)
*/

GeneralizedLogicalNetwork ReconstructPureGLN( const TrajectoryCollection & traCol, const bool & selfControl, const Topology & mustIncludeParents)
{cout<<"Start to recon...\n";
    vector< vector< vector<int> > > tra = traCol.getTrajTables();
    const vector<int> radixs( traCol.getBases() );
    vector<string> nodeNames( traCol.getIntNodeNames() );
    nodeNames.insert( nodeNames.end(), traCol.getExtNodeNames().cbegin(), traCol.getExtNodeNames().cend()  );
    const vector<bool> beParent( traCol.getBeParent() );
    const vector<bool> haveParent( traCol.getHaveParent() );
    const vector< vector<int> > excludeLists( traCol.getexcludeTrajectory() );

    if( tra.size() == 1 && tra[0].size() == 2)
    {
        cerr<<"ReconstructionPureGLN(): There are too few rows (only two rows) in the trajectory collection to be used to reconstruct the GLN.\nReconstructionPureGLN() requires at least three rows.\nMaybe try GLN::Reconstruction()\n"<<endl;
        exit(-1);
    }

    assert( tra[0][0].size() == radixs.size() );

    if( excludeLists.size() != 0 )
        assert( excludeLists.size() == tra.size() );

	vector< set<int> > parentLists;
	parentLists.resize( tra[0][0].size() );
    
    //Construct the parent set for each node
    for(size_t n = 0; n<tra.size(); ++n) //Check each trajectory
    {
        for(size_t time=2; time<tra[n].size(); ++time) //Check time points in each trajectory
        {
            set<int> possibleParentSet;   //Parent ID: 0-based
            for(int p=0; p<(int)tra[0][0].size(); ++p)
            {
                if( beParent[p] )
                {   
                    //If node p is exclude from trajectory n, then don't add p as parent
                    if( excludeLists.size() != 0 && find( excludeLists[n].begin(), excludeLists[n].end(),p) != excludeLists[n].end() ) 
                        continue;

                    //If p changed, add p to the possible parent set
                    if( tra[n][time-1][p] != tra[n][time-2][p]){
                        possibleParentSet.insert(p);
                    }
                }
            }
            
            //Add nodes in the possible parent set as parent of nodes who changed from time t-1 to t
            for(int child=0; child<(int)tra[0][0].size(); ++child)
            {
                if( ! haveParent[child] ) 
                {//child node can not have parent, set its parent set to empty
                    parentLists[child].clear();
                    continue;
                }

                //Child's value changed, merge possibleParentSet to this child's parent set
                if( tra[n][time][child] != tra[n][time-1][child]) {
                    parentLists[child].insert( possibleParentSet.begin(), possibleParentSet.end() );
                }
            }   
        }
    }


    /*Check to see if previous procedure chose wrong parent set.
    If wrong parent set is chosen, the transition from parents' combination to child may contain conflict.
    By looking at contingency table, we may found such conflict.
    If conflict is found, we need to add more nodes into the parent set.
    (Fictitious parent(s) will be eliminated later)
    */
    vector<TransitionTable> transitionTables;

    for(size_t t=0; t<tra[0][0].size(); ++t)
    {
        vector<int> tempParentSet( parentLists[t].begin(), parentLists[t].end() );
        
        for(size_t x=0; x<tempParentSet.size(); ++x)
            ++tempParentSet[x];

        transitionTables.push_back( TransitionTable((int)t,tempParentSet,vector<int>(tempParentSet.size(),-1), traCol)  );
    }

    /*
    for(size_t t=0; t<transitionTables.size(); ++t)
    {
        transitionTables[t].fill(traCol);
        
        while ( ! transitionTables[t].pure() )
        {
            for(size_t n=0; n<tra[0][0].size(); ++n)
            {
                if( parentLists[(int)t].find(n) == parentLists[(int)t].end() )
                {
                    if( beParent[n] )
                    {
                        parentLists[t].insert((int)n);
                        break;
                    }
                }
            }

            vector<int> tempParentSet( parentLists[t].begin(), parentLists[t].end() );
            for(size_t x=0; x<tempParentSet.size(); x++)
                ++tempParentSet[x];
            
            transitionTables[t] = TransitionTable((int)t,tempParentSet,vector<int>(tra[0][0].size(),-1), traCol);
            transitionTables[t].fill(traCol);
        }
    }
    */

    for(size_t t=0; t<transitionTables.size(); ++t)
    {
        transitionTables[t].fill(traCol);
        if( !transitionTables[t].pure() ){
            for(size_t n=0; n<tra[0][0].size(); ++n){
                if( parentLists[(int)t].find(n) == parentLists[(int)t].end() ){
                    if( beParent[n] ){
                        parentLists[t].insert((int)n);
                    
                        vector<int> tempParentSet( parentLists[t].begin(), parentLists[t].end() );
                        for(size_t x=0; x<tempParentSet.size(); ++x)
                            ++tempParentSet[x];
                    
                        transitionTables[t] = TransitionTable((int)t,tempParentSet,vector<int>(tra[0][0].size(),-1), traCol);
                        transitionTables[t].fill(traCol);
                        if( transitionTables[t].pure() ){
                            break;
                        }
                    }
                }
            }
        }
    }


    //If selfControl is false (default value), then a node can NOT be itself`s parent
    if( !selfControl )
    {
        for(int index=0; index<(int)radixs.size(); ++index)
        {//Check to see if one node is itself`s parent. If it is, delete it.
            if( parentLists[index].find(index) != parentLists[index].end() ){
                parentLists[index].erase(index);
            }
        }
    }

    


    //Added by Haizhou Wang, Feb 2, 2013
    //When given must-include parents for each node, Pure Reconstructon will include them in the final output
    if( mustIncludeParents.size() != 0){
        for(int node=0; node<(int)mustIncludeParents.get().size(); ++node){
            for(const auto & parentID: mustIncludeParents.get()[node][0] ){
                //If the must-include parent for this node in not in parent list, add it.
                if( parentLists[node].find( parentID - 1 ) == parentLists[node].end() ){
                    parentLists[node].insert( parentID - 1);
                }
            }
        }
    }







    //Once get the parent set for each node
    //Build truth value vector for each node
    vector<GeneralizedTruthTable> gtts;
    gtts.resize( tra[0][0].size() );

    for(size_t child =0; child<gtts.size(); ++child)
    {
        if( parentLists[child].size() == 0)
            continue;

        vector<int> parentIndex( parentLists[child].begin(), parentLists[child].end() );
        vector<int> parentBases;
        vector< vector<double> > truthValueProb;
        

        for(size_t r=0; r<parentIndex.size(); ++r)
            parentBases.push_back(radixs[ parentIndex[r] ]);


        int NumOfTruthValues = 1;
        for(size_t n=0; n<parentIndex.size(); ++n)
            NumOfTruthValues *= radixs[ parentIndex[n] ];

        vector<int> truthValue(NumOfTruthValues, -1); //-1 means that value is unknown.

        //Initialize truthValueProb, set all initial probabilities to 0
        truthValueProb.resize(NumOfTruthValues);

        for(int v=0; v<NumOfTruthValues; ++v)
        {    
            truthValueProb[v] = vector<double>(radixs[child], 0);
            //truthValueProb[v][0] = 1;
        }
        
        /*
        Since trajectory (trajectory collection) is usually shorter than the size of truth value vector (hundreds compare to millions), it's better to traverse the trajectory (trajectory collection) once and initialize the truth value vector other than traverse the truth value vector itself.
        If, trajectory is much longer than the size of truth value vector, should consider the other way around.
        */
        for(size_t n=0; n<tra.size(); ++n)
        {
            for(size_t row=0; row<tra[n].size()-1; ++row)
            {
                int childValue = tra[n][row+1][child];

                int index = getParentCombinationIndex( tra[n][row], parentIndex, parentBases);
                truthValue[index] = childValue;
                truthValueProb[index][childValue] = 1;
            }
        }

        for(const auto & row: truthValueProb )
            assert( std::accumulate( row.cbegin(), row.cend(), 0.0 ) <= 1.0 );


         vector<string> parentNames;
         for( const auto & index: parentIndex)
             parentNames.push_back( nodeNames[index] );

        //In the final gln file, node starts from 1. But in program, node starts from 0
        for(size_t t=0; t<parentIndex.size(); ++t)
            ++parentIndex[t];

        gtts[child] = GeneralizedTruthTable((int)child, radixs[child], parentIndex, parentBases, vector<int>(parentLists[child].size(),-1), truthValue, truthValueProb );
        gtts[child].setParentNames( parentNames );
    }


    GeneralizedLogicalNetwork gln;
    gln.setVersion("GENERALIZED_LOGICAL_NETWORK_VER2");
    gln.initialize( (int)radixs.size() );
    gln.setMaxMarkovOrder(-1);
    gln.setMinMarkovOrder(-1);

    for(size_t t=0; t<radixs.size(); t++)
    {
        if( haveParent[t] == true)
            gln.initializeNode((int)t, tra[0][0][t], radixs[t], nodeNames[t], 'i');
        else
            gln.initializeNode((int)t, tra[0][0][t], radixs[t], nodeNames[t], 'e');
    }

    gln.setGTTs( gtts );
    gln.eliminateFictitiousParents(mustIncludeParents);

    cout<<"Recon finished...\n";
    return gln;
}




GeneralizedLogicalNetwork ReconstructPureGLN( const GLNConfig & gc)
{
    EnvGLNRec envman;
    envman.initialize(gc);

    vector<string> names(envman.getTraCol().getIntNodeNames());
    
    /*
    Takes care of trajectory version 1, which doesn't have node names
    Will use node ID as node's name
    */
    size_t t=0;
    if( names.size() == 0) 
    {
        for(t=0; t<envman.getTraCol().getBases().size()-envman.getTraCol().getExtNodes().size(); ++t)
        {
            stringstream ou;
            ou<<(t+1);
            names.push_back( ou.str() );
        }
        envman.getTraCol().setIntNodeNames( names );
    }
    


    names.clear();
    names = envman.getTraCol().getExtNodeNames();
    if( names.size() == 0) 
    {
        for(; t<envman.getTraCol().getBases().size(); ++t)
        {
            stringstream ou;
            ou<<(t+1);
            names.push_back( ou.str() );
        }
        envman.getTraCol().setExtNodeNames( names );
    }
   


    //return ReconstructPureGLN( envman.getTraCol(), envman.getTraCol().getBases(), names, envman.getTraCol().getBeParent(), envman.getTraCol().getHaveParent(), envman.getTraCol().getexcludeTrajectory(),gc.m_allowSelfCycle );
    return ReconstructPureGLN( envman.getTraCol(), gc.m_allowSelfCycle, gc.m_mustIncludeParents);
}
