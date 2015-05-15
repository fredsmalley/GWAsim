#include "kStop_TestCase.h"


#include <algorithm>
#include <cassert>
#include <functional>
#include <numeric>

#include <map>
using std::map;

#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;
using std::to_string;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include "DS_LocalConstraint.h"
#include "DS_CompulsiveParents.h"
#include "EffCGLN_GlobalConstraint.h"
#include "GLN.h"
#include "GLNGlobals.h" //To use Pure Reconstruction
#include "kStop.h"
#include "TrajectoryCollection.h"
#include "Topology.h"



const DS_LocalConstraint kStop_TestCase::extractLocalConstraints(const GeneralizedLogicalNetwork & gln, const double percentage)
{
    DS_LocalConstraint lc;

    for(unsigned nameIndex=0; nameIndex<gln.size(); ++nameIndex){
        lc.m_nodes.push_back( gln.getNode(nameIndex).getName() );
    }

    for(const auto & node: gln.getGTTs() )
    {
        lc.m_radixs.push_back( node.getBase() );        
        

        vector<int> parentIDs( node.getParents() );  //Parent IDs, 1-based
        random_shuffle( parentIDs.begin(), parentIDs.end(), RNG); //RNG: Random Number Generator, defined in kStop.h
  

        //Select a random number of parents for this child
        int numOfParents = static_cast<int>( node.getParents().size() * percentage );

        DS_CompulsiveParents parentList;
        parentList.m_nodeList.assign( parentIDs.begin(), parentIDs.begin() + numOfParents );


        /*
        *In GenerlizedTruthTable class, only child's value is saved in m_truthValues
        *In order to extract the relationship between parent and child as pairwise local constraint file,
        *a vector of parent value must be rebuilt first.
        *Then based on the correlation between parent vector and child vector,
        *the control (positive/negative) from parent to child could be decided.
        */
        for(const auto & p : parentList.m_nodeList )
        {
            int index = -1;
            for(int t=0; t<(int)node.getParents().size(); ++t)
            {
                if( node.getParents()[t] == p )
                {
                    index = t;
                    break;
                }
            }

            if( index == -1 )
            {
                GLNExit( EXIT_FAILURE, "Can NOT locate parent's ID in child's parent list when extracting local constraints, something goes wrong!\n");
            }

            unsigned radixBehind;
            if( index == (int)node.getParents().size() )
                radixBehind = 1;
            else
                radixBehind = accumulate( node.getParentBases().begin()+index+1, node.getParentBases().end(), 1, std::multiplies<int>() );

            vector<int> parentValueList( node.getTruthValues().size() );
            for(size_t t=0, parentValue=0, iteration=0; t<parentValueList.size();)
            {
                while( iteration < radixBehind )
                {
                    parentValueList[t++] = parentValue;
                    ++iteration;
                }
                iteration = 0;
                parentValue = ( ++parentValue == (unsigned) node.getParentBases()[index]) ? 0 : parentValue;
            }

            double correlation = DS_LocalConstraint::correlation( parentValueList, node.getTruthValues() );
            parentList.m_polarityRelations.push_back( correlation );

        }//end of for() loop, proceed to next parent of node

        lc.m_compulsiveParents.push_back( parentList );

    }//end of outer for() loop, proceed to next node in GLN

    return lc;
}


const EffCGLN_GlobalConstraint kStop_TestCase::extractGlobalConstraints( const TrajectoryCollection & TraCol, const int numOfGlobalConstraints, const bool incompleteDynamicPattern, const bool sequential)
{
    EffCGLN_GlobalConstraint gc;
    gc.m_nodes.insert( gc.m_nodes.end(), TraCol.getIntNodeNames().begin(), TraCol.getIntNodeNames().end() );
    gc.m_nodes.insert( gc.m_nodes.end(), TraCol.getExtNodeNames().begin(), TraCol.getExtNodeNames().end() );
    gc.m_radixs = TraCol.getBases();


    //Check to see if the given trajectory contains loop/cycle
    vector< vector<int> > tra( TraCol.getTrajTables()[0] );
    int distance = tra.size();

    //Decide how many distinct rows are there in the trajectory.
    map< vector<int>, int> row2Index;
    for( int index=0; index<(int)tra.size(); ++index)
    {
        if( row2Index[tra[index]] == 0)
            ++row2Index[tra[index]];
        else
        {
            distance = index - 1;
            break;
        }
    }

    if( distance == 1)
    {
        GLNExit( EXIT_FAILURE, "kStop_TestCase::extractGlobalConstraints(): When extracting global constraints, the simulated trajectory of the original GLN has only one distinct row.\nThus, only one dynamic pattern can be extracted. Which is meaningless in the application of this program.");
    }
    else if( numOfGlobalConstraints > distance )
    {
        GLNExit( EXIT_FAILURE, "kStop_TestCase::extractGlobalConstraints(): number of distinct rows in simulated trajectory is less than the number of global constraints you want to extract.\n");
    }




    
    vector<unsigned> rows;
    for(int t=0; t<distance; ++t)
        rows.push_back(t);
    
    // sequential = true, extract consecutive rows from the begin of trajecotry.
    vector<unsigned> extraction;
    if( !sequential )
    {
        random_shuffle( rows.begin(), rows.end(), RNG); //RNG: Random Number Generator, defined in kStop.h
        extraction.assign( rows.begin(), rows.begin() + numOfGlobalConstraints );
        sort( extraction.begin(), extraction.end() );
    }
    else
    {
        extraction.assign( rows.begin(), rows.begin() + numOfGlobalConstraints );
    }




    EffCGLN_DynamicPattern dynamicPatterns;

    vector<unsigned> cols;
    for(size_t t=0; t<TraCol.getTrajTables()[0][0].size(); ++t)
        cols.push_back(t);

    for(const auto & i: extraction)
    {
        vector<int> pulledRow( TraCol.getTrajTables()[0][i]  );
        
        //Incomplete Dynamic Pattern: within each dynamic pattern, some bits are set to -1 (wild card)
        if( incompleteDynamicPattern ) 
        {
            random_shuffle( cols.begin(), cols.end(), RNG);

            // Randomly choose 0~50% bits to hide
            int numOfColsToUnmark = RNG(  (int)(TraCol.getTrajTables()[0][0].size()/2)  );
            for(auto t=0; t<numOfColsToUnmark; ++t)
            {
                pulledRow[ cols[t] ] = -1;
            }

        }

        dynamicPatterns.m_DynamicPattern.push_back( pulledRow );
    }

    gc.m_DynamicPatternList.push_back( dynamicPatterns );
    return gc;
}


//Amaong a list of GLNs, sort each edge by its occurance
const vector< pair< pair<int,int>, unsigned> > kStop_TestCase::RankEdges( const vector<GeneralizedLogicalNetwork> & GLNList )
{
    map< pair<int,int>, unsigned> occurance;

    for( const auto & GLN: GLNList){
        for( const auto & node: GLN.getGTTs() ){
            for( const auto & parent: node.getParents() ){
                //edge pair: <parent, child>
                pair<int, int> edge( parent, node.getChild()+1);  //edge<1-based, 1-based>
                ++occurance[ edge ] ;                
            }
        }
    }

    vector< pair<pair<int,int>, unsigned> > result( occurance.begin(), occurance.end() );
    //Sort descendly based on the occurance of each edge. Can't use map, since it's sorted by key
    sort( result.begin(), result.end(), kStop_TestCase::comp );

    return result;
}


inline const int kStop_TestCase::findEdgeCount(const vector< pair< pair<int,int>, unsigned> > & edgeCounts, const pair<int,int> & e)
{
    for(const auto & entry: edgeCounts){
        if( entry.first == e ){
            return entry.second;
        }
    }

    return -1;
}

//Among a list of GLNs, rank each GLN by the occurance of its median edge
const double kStop_TestCase::getGLNRankByEdgeMedian(const vector< pair< pair<int,int>, unsigned> > & edgeCounts, const GeneralizedLogicalNetwork & gln)
{
    vector< pair< pair<int,int>, unsigned > > edgeRanks;
    for(const auto & e: gln.getEdgeList() )
    {
        int count = kStop_TestCase::findEdgeCount( edgeCounts, e );
        if( count == -1 )
        {
            cerr<<"Error! kStop_TestCase::rankGLNByEdgeMedian(): edge<"<<e.first<<","<<e.second<<"> is NOT found in the edge list\n";
            exit(-1);
        }
        edgeRanks.push_back( make_pair(e,count) );
    }

    //Sort descendly based on the occurance of each edge.
    sort( edgeRanks.begin(), edgeRanks.end(), kStop_TestCase::comp );

    if( edgeRanks.size()%2 == 0 )
    {
        int upper_middle = edgeRanks.size()/2;
        int lower_middle = upper_middle - 1;
        assert( lower_middle >= 0 );

        return ( (edgeRanks[upper_middle].second + edgeRanks[lower_middle].second)/2.0);
    }
    else
    {
        int middle = (edgeRanks.size()-1)/2;
        return ( edgeRanks[middle].second );
    }
}



void kStop_TestCase::outputGLNsRankedByEdgeMedian( const vector< pair< pair<int,int>, unsigned> > & edgeCounts, const vector<GeneralizedLogicalNetwork> & GLNList, const vector<string> & GLNNames)
{
    assert( GLNList.size() == GLNNames.size() );

    vector< pair<int,double> > scoresByIndex;
    for( size_t i=0; i<GLNList.size(); ++i)
    {
        scoresByIndex.push_back( make_pair(i,  kStop_TestCase::getGLNRankByEdgeMedian( edgeCounts, GLNList[i] )) );
    }

    //Sort descendly based on the score.
    sort( scoresByIndex.begin(), scoresByIndex.end(), kStop_TestCase::comp2 );

    ofstream ou("GLNList_Ranked_By_Edge_Median.txt");
    ou<<"GLN Name\tMedian Edge Score"<<endl;

    for(const auto & r: scoresByIndex)
    {
        ou<<GLNNames[ r.first ]<<"\t"<<r.second/GLNList.size()*100<<"%"<<endl;
    }


    ou.close();
}


void kStop_TestCase::test(const string & inputGLNFileName, const bool incompleteDynamicPattern, const double percentage, const int numOfGlobalConstraints, const string SearchDirection, const int depth )
{
    /*
    * 1. Do simualation on the input GLN, get a trajectory
    * 2. Based on input variable 'percentage' (hidden percentage), extract local constraint
    * 3. Based on the input variable 'numOfGlobalConstraints', extract global constraint
    * 4. Feed local and global constraints to kStops program to get a simple path/trajectory if one exists
    * 5. Reconstruct a GLN out of the result of previous step
    * 6. Compare the reconstructed GLN with the original ipnut GLN
    */
    GeneralizedLogicalNetwork gln;
    gln.scan( inputGLNFileName + ".gln" );


/*********************************************************************************************/
    int timeSteps = 1;
    for(size_t t=0; t<gln.size(); ++t)
        timeSteps *= gln.getNode(t).getBase();

    if( timeSteps < 2){
        GLNExit( EXIT_FAILURE, "The original GLN is too simple.\nThere are too few time points to be simulated\n" );
    }

    //simulate(inittrjfile, numberoftrj, timesteps,levelNoise=0,houseModel="default")
    TrajectoryCollection traCol = gln.simulate("", 1, timeSteps);
/*********************************************************************************************/


/*********************************************************************************************/
   DS_LocalConstraint lc( kStop_TestCase::extractLocalConstraints( gln, percentage) );
   EffCGLN_GlobalConstraint gc( kStop_TestCase::extractGlobalConstraints( traCol, numOfGlobalConstraints, incompleteDynamicPattern) );

   //lc.saveLocalConstraint(inputGLNFileName + "ExtractedLocalConstraintFile.txt");
   //gc.saveGlobalConstraint(inputGLNFileName + "ExtractedGlobalConstraintFile.txt");

   gc.saveAsTraColFile( inputGLNFileName + "Tra.gln");
   lc.saveAsTopologyFile( inputGLNFileName + "LocalConstraintTopology.txt");

   //To get a truth edge list, for drawing PR and ROC curve later
   kStop_TestCase::extractLocalConstraints(gln, 1).saveConstraints(inputGLNFileName + "TruthLocalConstriants.txt");
   lc.saveLocalConstraintPair( inputGLNFileName + "LocalConstraints.txt");//To be used to remove extracted edges when plotting PR/ROC
/*********************************************************************************************/
   //DS_LocalConstraint lc("ExtractedLocalConstraintFile.txt");
   //EffCGLN_GlobalConstraint gc("ExtractedGlobalConstraintFile.txt");
    

/*********************************************************************************************/
//Direct reconstruction of the extracted global cosntraints without interpolation
    GLNConfig directRecon;
    directRecon.m_alpha = 1;
    directRecon.m_allowSelfCycle = true;

    vector<string> tmp;
    tmp.push_back( inputGLNFileName +  "Tra.gln" );
    directRecon.m_listTrajFiles = tmp ;
    
    TrajectoryCollection trajCol;
    trajCol.scan( inputGLNFileName + "Tra.gln" );
    directRecon.m_BP = trajCol.getBeParent();
    directRecon.m_HP = trajCol.getHaveParent();

    //Unconstrained Network Modeling
    GeneralizedLogicalNetwork directPureWithoutLocalGLN( ReconstructPureGLN( directRecon  ) ) ;


    //Topology-constrained Network Modeling
    directRecon.m_mustIncludeParents.setNodeNames( trajCol.getIntNodeNames() );
    directRecon.m_mustIncludeParents.scan( inputGLNFileName + "LocalConstraintTopology.txt"); 
    directRecon.m_mustIncludeParents.setFile( inputGLNFileName + "LocalConstraintTopology.txt");
    GeneralizedLogicalNetwork directPureWithLocalGLN( ReconstructPureGLN( directRecon  ) ) ;
/*********************************************************************************************/


/*********************************************************************************************/
//kStop, dynamic-constrained network modeling
    vector<GeneralizedLogicalNetwork> GLNList;
    
    //If there are too many satisfying candidates, only use the first 100.
    //If there are no satisfying candidates, stop after trying 100 times.
    int t = 0;
    while( t< 100 )
    {
        //kStop::kStop(local_constraint, global_constraint, traverse_order, deepth_limit_for_IDDFS)
        kStop testExample( lc, gc, SearchDirection, depth);
        
        if ( !testExample.FindPath() ){
            cerr<<"Test failed.\n No solution found by CGLN with search detph "<<depth<<".\n Considering increase the search detph if viable.\n"<<endl;
            return;
        }
        
        if( !testExample.getTrajectoryCollection().CheckTraDeterministic() ){
            GLNExit( EXIT_FAILURE, "kStop() produced an impure trajectory collection (It shouldn't).\n Program will be terminated.\n");
        }


        testExample.saveTra(inputGLNFileName + "-" + to_string(t) + ".tra");

        GeneralizedLogicalNetwork tmp( ReconstructPureGLN( testExample.getGLNConfig(inputGLNFileName, to_string(t) ) ) );
        
        if( find( GLNList.cbegin(), GLNList.cend(), tmp ) != GLNList.cend() )
            break; //Duplicated gln reconstructed, stop.
 
        GLNList.push_back( tmp );  
        tmp.save( inputGLNFileName + "-PureRecon-" + to_string(t) + ".gln");

        ++t;
    }
    assert( GLNList.size() != 0 );
    
/*********************************************************************************************/


/********************************************************************************************/

    //<parent,child>:occurrence
    vector< pair< pair<int,int>, unsigned> > trueEdgeRank( kStop_TestCase::RankEdges( vector<GeneralizedLogicalNetwork>(1, gln) ) );
    vector< pair< pair<int,int>, unsigned> > CGLNEdgeRank( kStop_TestCase::RankEdges( GLNList ) );
    vector< pair< pair<int,int>, unsigned> > directGLNWithLocalEdgeRank( kStop_TestCase::RankEdges( vector<GeneralizedLogicalNetwork>(1, directPureWithLocalGLN) ) );
    vector< pair< pair<int,int>, unsigned> > directGLNWithoutLocalEdgeRank( kStop_TestCase::RankEdges( vector<GeneralizedLogicalNetwork>(1, directPureWithoutLocalGLN) ) );


    //<parent,child> both 1-based
    set< pair<int, int> > trueEdgeList;
    for( const auto & e: trueEdgeRank)
        trueEdgeList.insert( e.first );

    vector< pair<int,int> > CGLNEdgeList;
    for( const auto & e: CGLNEdgeRank)
        CGLNEdgeList.push_back( e.first );

    vector< pair<int, int> > directGLNWithLocalEdgeList;
    for( const auto & e: directGLNWithLocalEdgeRank)
        directGLNWithLocalEdgeList.push_back( e.first );

    vector< pair<int, int> > directGLNWithoutLocalEdgeList;
    for( const auto & e: directGLNWithoutLocalEdgeRank)
        directGLNWithoutLocalEdgeList.push_back( e.first );

    int totalRows = gln.size() * gln.size() ;

    output_PR_ROC( totalRows, trueEdgeList, CGLNEdgeList, "pureCGLN" + inputGLNFileName);
    output_PR_ROC( totalRows, trueEdgeList, directGLNWithLocalEdgeList, "directPureWithLocalGLN" + inputGLNFileName);
    output_PR_ROC( totalRows, trueEdgeList, directGLNWithoutLocalEdgeList, "directPureWithoutLocalGLN" + inputGLNFileName);
//*********************************************************************************************
}



//vector< <parentID,childID> >
void kStop_TestCase::output_PR_ROC(const unsigned totalRows, const set< pair<int,int> > & truthList, const vector< pair<int,int> > & predictList, const string & fileName)
{
	ofstream ou(fileName+"-PR.txt");
	ofstream ou2(fileName+"-ROC.txt");
	
	ou<<"Precision\tRecall\n";
	ou2<<"TPR\tFPR\n";
	
	for(size_t row=0; row<predictList.size(); ++row)
	{
		unsigned tp=0, tn=0, fp=0, fn=0;
		
		for(size_t t=0; t<=row; ++t){
			if( truthList.find( predictList[t] ) != truthList.cend() ){
				++tp;
            }
			else{
				++fp;
            }
		}
		
        fn = truthList.size() - tp;
        tn = totalRows - (row+1) - fn;

        assert( fn>=0 && tn>=0 && fp>=0 && tp>=0);
        assert( tp+fp == row+1 );
		
		double precision = double(tp)/(tp+fp);
        double recall = double(tp)/truthList.size();
		
		double tpr = recall;
		
        double fpr;
        if( (fp+tn) == 0)
            fpr = 0;
        else
            fpr = double(fp)/(fp+tn);
		
		ou<<precision<<"\t"<<recall<<endl;
		ou2<<tpr<<"\t"<<fpr<<endl;
	}

    ou<<"1\t1\n";
    ou2<<"1\t1\n";

	ou.close();
	ou2.close();
}
