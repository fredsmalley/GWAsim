// EnumMethod.cpp
//
// Joe Song
// Created: November 30, 2008
// Last modified: MS November 21, 2011.  Added new functions enumerateParentSubsets(),
//	next_N_Elements_To_K_Sets(), nextAssignmentWithRepeat(), nextAssignmentNoRepeat().
//  Made older functions enumerateParents() and enumerateParents2() obsolete.

#include <cassert>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>

#include "EnumMethod.h"
#include "SetOps.h"

int EnumMethod::next_N_Elements_To_K_Sets(vector<unsigned long> & a, unsigned K)
{
    // Each element is assigned to a K-bit number 1 to 2^K-1 except 0
    //
    // Enumerate a size N integer vector, each element goes from 1 to 2^K-1
    //   and repeat is allowed
    //
    // Return a negative value when the last enumeration in the sequence
    //   has been passed and the assignment will go to the first enumeration
    //
    // Example: n=3
    // a= (  1,   1,   1)
    //    (  1,   1,   2)
    //       .
    //       .
    //       .
    //    (2^K-1, 2^K-1, 2^K-1)
    
    unsigned N = a.size();
    
    int i;
    
    for(i = N-1; i >= 0; i --) {
        if(a[i] < static_cast<unsigned>( (1 << K) - 1) ) {
            a[i] ++;
            break;
        } else {
            a[i] = 1;
        }
    }
    
    return i; // i<0 implies all combination has been enumerated once
}

int EnumMethod::nextAssignmentWithRepeat(vector<unsigned> & a, unsigned K)
{
    // Each of N elements is assigned exactly one number between 0 and K-1
    //
    // Enumerate a size N integer vector, each element goes from 0 to K-1
    //   and repeat is allowed
    //
    // Return a negative value when the last enumeration in the sequence
    //   has been passed and the assignment will go to the first enumeration
    //
    // Example: n=3
    // a= (  0,   0,   0)
    //    (  0,   0,   1)
    //       .
    //       .
    //       .
    //    (K-1, K-1,  K-1)
    
    unsigned N = a.size();
    
    int i;
    
    for(i=N-1; i>=0; i--) {
        if(a[i] < K-1) {
            a[i] ++;
            break;
        } else {
            a[i] = 0;
        }
    }
    
    return i; // i<0 implies all combination has been enumerated once
}

int EnumMethod::nextAssignmentWithRepeat(vector<int> & a, const vector<int> & K)
{
    // Each of N elements is assigned exactly one number between 0 and K-1
    //
    // Enumerate a size N integer vector, each element goes from 0 to K-1
    //   and repeat is allowed
    //
    // Return a negative value when the last enumeration in the sequence
    //   has been passed and the assignment will go to the first enumeration
    //
    // Example: n=3
    // a= (  0,   0,   0)
    //    (  0,   0,   1)
    //       .
    //       .
    //       .
    //    (K[0]-1, K[1]-1,  K[2]-1)
    
    unsigned N = a.size();
    
    int i;
    
    for(i=N-1; i>=0; i--) {
        if(a[i] < K[i]-1) {
            a[i] ++;
            break;
        } else {
            a[i] = 0;
        }
    }
    
    return i; // i<0 implies all combination has been enumerated once
}


int EnumMethod::nextAssignmentNoRepeat(vector<int> & a, int K)
{
    // Enumerate the next size m vector, each element goes from 1 to at most K
    //   and no repeat is allowed.
    //
    // Return a negative value when the last enumeration in the sequence
    //   has been passed and the assignment will go to the first enumeration
    //
    // Example: N=3
    // a=   1, 2, 3
    //      1, 2, 4
    //         .
    //         .
    //         .
    //    K-2, K-1, K
    
    int N = a.size();
    
    int position = N - 1;
    
    while(position >= 0){
        
        //Check to see if current position has reached it's max
        if(a[position] < static_cast<int>(K - (N - position - 1)) ){
            break;
        }//end if
        position--;
        
    }//end while
    
    if(position >= 0){
        
        a[position]++;
        
        if(position < N - 1) {
            
            //fill in the rest of the positions
            for(int i = position + 1; i < N; i++) {
                
                a[i] = a[i-1]+1;
                
            } //end for
        }//end if
    }
    
    return position;
}

bool EnumMethod::enumerateParentSubsets(const int child,
                                        const vector<int> & parentsBegin,
                                        const vector<int> & parentsEnd,
                                        const int nNodes,
                                        const GLNConfig & gc,
                                        const vector<int>  & parentCandidates,
                                        const vector<int> & mustIncludeParents)
{
    // enumerateParentSubsets() allows different parent sets to be used for each
    //   condition under comparison.  It applies to two or more experimental
    //   conditions.
    
    // INPUT:
    // 'parentsBegin' and 'parentsEnd' give the randge of possible parent set.
    //   This functoin will go through each combanation inbetween the range
    //   and find the best parent set as defined by the statistical test used.
    
    // 'nNodes' is the number of nodes who can be parent.
    
    // 'parentCandidates' records nodes' ID (1-based) who can be parent,
    //   it is the upper bound of parents.
    
    // 'mustIncludeParents' contains the lower bound of parents
    //   nodes in 'mustIncludeParents' must be included.
    
	bool improved = false;
    
    size_t nParents = parentsBegin.size();
    
	if( nParents <= 0) {
		return improved;
	}
    
    assert( parentsBegin.size() == parentsEnd.size() );
    
	vector<int> parentSuperset( parentsBegin );
    vector<int> mappedParentSuperset( nParents );
    
    int position;
    
	do {
        
        if(! parentCandidates.empty() ) {
            
            // Convert the parent id to those on the list of candidate parents
            for(size_t i=0; i<parentSuperset.size(); ++i) {
                mappedParentSuperset[i] = parentCandidates[parentSuperset[i] - 1];
            }
            
        } else {
            //If no candidate parent set is given, start with the
            // lower bound of possible parent sets.
            mappedParentSuperset = parentSuperset;
        }
        
        //*****************************************************************************************
        //Added by Haizhou Wang, Jan 23, 2013
        if( ! mustIncludeParents.empty() ) {
            
            // Check to see if all must-included parents are in 'mappedParentSuperset'.
            
            bool mustParentsContained = true;
            for(const auto & parent: mustIncludeParents)
            {
                if( find( mappedParentSuperset.cbegin(), mappedParentSuperset.cend(), parent)
                   ==
                   mappedParentSuperset.cend() )
                {
                    mustParentsContained = false;
                    break;
                }
            }
            
            if( ! mustParentsContained ) { //If not, enumerate next 'mappedParentSuperset'.
                position = nextAssignmentNoRepeat(parentSuperset, nNodes);
                continue;
            }
            
            // If yes, proceed to enumerateDelays().
        }
        //*****************************************************************************************
        
        
		if(! gc.allowSelfCycle()) { // Check for self cycle
            //If self cycle is NOT allowed
            if(mappedParentSuperset.end()
               != find(mappedParentSuperset.begin(), mappedParentSuperset.end(), child+1)) {
                
                position = nextAssignmentNoRepeat(parentSuperset, nNodes);
                continue;
            }
		}
        
        if((gc.m_listTrajFiles.size() > 1 && gc.m_acceptableTopologies == DIFFERENT) // 2)
           ||gc.m_candidateTopologys.size()>1||gc.m_pathways.size()>1) {
            // GLN Comparison different parents
            
            unsigned nConditions = gc.m_listTrajFiles.size();
            
            vector<unsigned long> assignment(nParents, 1); // Bit j of
            // element(parent) in assignment indicates whether that parent i
            // is included in the subset j or not
            
            do { // Enumerate all possible parent subsets, one for each condition,
                //   that together cover the parent superset
                
                vector<vector<int> > parentSubsets(nConditions, vector<int>(0));
                
                // The for-loop translate the assignment to actual parents
                for(size_t m=0; m<nParents; m++) {
                    unsigned long mask = 1 << (nConditions-1);
                    for(size_t k=0; k<nConditions; k++) {
                        if(assignment[m] & mask) {
                            parentSubsets[k].push_back(mappedParentSuperset[m]);
                        }
                        mask = mask >> 1;
                    }
                }
                
                if( withinRange(parentSubsets) ) {
                    // Only continue when sizes of parent subsets are within
                    //   required range specified by the user command line
                    if(child == -1) {
                        for (size_t k=0; k<parentSubsets.size(); k++) {
                            for (size_t j=0; j<parentSubsets[k].size(); j++) {
                                cout << parentSubsets[k][j] << ",";
                            }
                            cout << " vs ";
                        }
                        cout << endl; 
                    }
                    bool updated = enumerateDelays2(child, parentSubsets, gc);
                    improved = updated ? updated : improved;
                }
                
            } while (next_N_Elements_To_K_Sets(assignment, nConditions) >= 0);
            
        } else {
            // Covers:
            //   1. GLN Reconstruction
            //   2. GLN Comparison with same topology
            
            bool updated = enumerateDelays(child, mappedParentSuperset, gc);
            improved = updated ? updated : improved;
        }
        
        position = nextAssignmentNoRepeat(parentSuperset, nNodes);
        
	} while( ( position>=0 || parentSuperset[0] != nNodes - static_cast<int>(nParents - 1)  )
            && parentSuperset[position] <= parentsEnd[position]
            );
    
	return improved;
}


bool EnumMethod::enumerateParents2(const int child,
						 		   const vector<int> & parentsBegin,
								   const vector<int> & parentsEnd,
								   const int nNodes,
								   const GLNConfig & gc,
								   const vector<int>  & listCandidateParents)
{
    // enumerateParents2() allows different parent sets to be used for each condition under comparison.
    //   This code only applies when there are two experimental conditions.
    
	bool improved = false;
    
	if(parentsBegin.size()<=0)
	{
		return improved;
	}
    
	int nParents = (int) parentsBegin.size();
	vector<int> parents = parentsBegin;
    
	int position;
    
	do{
        
		int j = nParents;
        
		if(! gc.allowSelfCycle()) {
            
            vector<int> mappedParents(parents.size());
            
            if(! listCandidateParents.empty() ) {
                for(size_t i=0; i < mappedParents.size(); i++)
                {
                    mappedParents[i] = listCandidateParents[parents[i]-1];
                }
            }
            
			for(j=0; j < nParents ; j++) {
				//added by Yang Zhang 7.25.2010 to convert candidate parent ids back to GLN id
				if( listCandidateParents.size()!=0 )
				{
					if( child == mappedParents[j] - 1 ) {
						break;
					}
				}else
				{
					if( child == parents[j] - 1 ) {
						break;
					}
				}
			}
		}
        
		if( j == nParents ) {
			
			//added by Yang Zhang 2.16.2010
			vector<vector<int> > powerset =  powerSet(parents);
			//vector<vector<vector<int> > > enmuerateSet;
            
            vector<vector<int> > OneEnmueration(2); //need to change 2 to enable multiple condition comparison
            
            vector<vector<int> > mappedOneEnmueration(OneEnmueration.size());
            
			for(size_t i=0; i<powerset.size(); i++) {
                
				vector<vector<int> > powersetforsubset = powerSet(powerset[i]);
				for(size_t j = 0; j< powersetforsubset.size(); j++)
				{
					OneEnmueration[0] = powerset[i];
					vector<int> complementset = complementSet(powerset[i],parents);
					OneEnmueration[1] = unionSet(complementset,powersetforsubset[j]);
                    
					if(listCandidateParents.size()!=0) //using candidate parent file to enumerate parents, added by Yang Zhang 7.25.2010
					{
						for(size_t i=0;i<OneEnmueration.size();i++)
						{
							mappedOneEnmueration[i].resize(OneEnmueration[i].size());
							for(size_t j=0;j<mappedOneEnmueration[i].size();j++)
							{
								mappedOneEnmueration[i][j] = listCandidateParents[OneEnmueration[i][j]-1];
							}
						}
                        
						bool updated = enumerateDelays2(child, mappedOneEnmueration, gc);
						improved = updated ? updated : improved;
					}
					else
					{
						bool updated = enumerateDelays2(child, OneEnmueration, gc);
						improved = updated ? updated : improved;
					}
				}
			}
            
		}
        
		position = nParents - 1;
        
		while(position >= 0){
			//Check to see if current m_position has reached it's max
			if(parents[position] < (nNodes - (nParents - position - 1)) ){
				break;
			}//end if
			position--;
		}//end while
        
		if(position >= 0){
            
			parents[position]++;
            
			if(position < nParents-1) {
                
				//fill in the rest of the m_indexParents
				for(int i = position+1; i < nParents; i++) {
                    
					parents[i] = parents[i-1]+1;
                    
				} //end for
			}//end if
		}
	} while((position >=0 || parents[0] != nNodes - (nParents - 1))&&parents[position]<=parentsEnd[position] );
    
	return improved;
}


/*bool EnumMethod::enumerateDelays2(const int child,
 const vector<int> & parents1, const vector<int> & parents2,
 const GLNConfig & ep)
 {
 bool improved = false;
 
 vector<int> delays1(parents1.size());
 vector<int> delays2(parents2.size());
 
 bool first = true;		// bool to indicate to generate the first set of delays
 bool last1;	// bool to indicate last enumerate has reached.
 bool last2;	// bool to indicate last enumerate has reached.
 
 do{
 // Generate parent time offsets, until all possible
 // offsets have been gone through, last = true
 // when this has happened
 
 // generate_new = false;
 
 last1 = generateDelays(delays1, gc.m_min_Markov_order, gc.m_max_Markov_order, first);
 last2 = generateDelays(delays2, gc.m_min_Markov_order, gc.m_max_Markov_order, first);
 
 bool updated = processOneEnumeration(child, parents1, parents2, delays1, delays2);
 
 improved = updated ? updated : improved;
 
 first = false;
 
 } while( ! (last1 && last2) );
 
 return improved;
 }*/
//this function only enumerate same delays for parents,
//modified by Yang 2.23.2011
/*
 bool EnumMethod::enumerateDelays2(const int child,
 const vector<vector<int> > & parents,const GLNConfig & ep)
 {
 bool improved = false;
 
 vector<vector<int> > delays(parents.size());
 for(size_t i=0; i<parents.size(); i++)
 {
 delays[i].resize(parents[i].size());
 }
 
 bool first = true;		// bool to indicate to generate the first set of delays
 bool last = false;	// bool to indicate last enumerate has reached. changed from true to false by YZ 2010.12.6
 
 do{
 // Generate parent time offsets, until all possible
 // offsets have been gone through, last = true
 // when this has happened
 
 for(size_t i=0; i<delays.size(); i++)
 {
 //	last = last && (generateDelays(delays[i], gc.m_min_Markov_order, gc.m_max_Markov_order, first));
 //  modified by Yang Zhang, is last is false
 if(delays[i].size()==0)
 {
 continue;
 }
 last = generateDelays(delays[i], gc.m_min_Markov_order, gc.m_max_Markov_order, first);
 }
 
 bool updated = processOneEnumeration(child, parents, delays);
 
 improved = updated ? updated : improved;
 
 first = false;
 
 } while( ! last );
 
 return improved;
 }*/

bool EnumMethod::enumerateDelays2(const int child,
								  const vector<vector<int> > & parentSubsets,
                                  const GLNConfig & gc)
{
    // enumerateDelays2() allows different parent sets to be used for each condition under comparison
    
    /*
     if(child == 1) {
     
     cout << "[";
     for(int j=0; j<parentSubsets.size(); j++) {
     cout << "C" << j << ": ";
     for(int i=0; i<parentSubsets[j].size(); i++) {
     cout << parentSubsets[j][i];
     if(i < parentSubsets[j].size()-1) {
     cout << ',';
     }
     }
     if(j < parentSubsets.size()-1) {
     cout << "; ";
     }
     }
     cout << "]" << endl;
     
     }
     */
    
	bool improved = false;
    
	vector<vector<int> > delays(parentSubsets.size());
    
	vector<vector<int> > delaysFromPathway(parentSubsets.size());
	
	for(size_t i=0; i<parentSubsets.size(); i++)
	{
		delays[i].resize(parentSubsets[i].size());
		delaysFromPathway[i].resize(parentSubsets[i].size());
	}
    
    const vector<vector<vector<int> > > * pTopology = NULL;
    
	if(gc.m_pathways.size()==1)
	{
		pTopology = & (gc.m_pathways[0].get());
        
	} else if(gc.m_candidateTopologys.size()==1)
	{
		pTopology = & (gc.m_candidateTopologys[0].get());
	}
    
    for(size_t indexcond = 0; indexcond < parentSubsets.size(); indexcond ++) //index for conditions
    {
        
        if(gc.m_pathways.size()>1)
        {
            
            pTopology = & (gc.m_pathways[indexcond].get());
            
        }else if(gc.m_candidateTopologys.size()>1)
        {
            pTopology = & (gc.m_candidateTopologys[indexcond].get());
        }
        
        if(pTopology)
        {
            for(size_t indexparent = 0; indexparent < parentSubsets[indexcond].size(); indexparent++ )
            {
                for(size_t indexcand = 0; indexcand < (*pTopology)[child][0].size(); indexcand ++ )
                {
                    //index for candidates or pathways
                    if(parentSubsets[indexcond][indexparent] == (*pTopology)[child][0][indexcand])
                    {
                        //changed by Yang 2013.5.21
                        //before delay for both condition are all saved, now one topology file only save delay info for one condition
                        //e.g. before suppose parents are 1,2, delays by default are saved as 9,9,9,9
                        //now becomes to 9,9
                        //int delayindex = (parentSubsets.size())*indexcand + indexcond;
                        //delaysFromPathway[indexcond][indexparent] = (*pTopology)[child][1][delayindex];
                        delaysFromPathway[indexcond][indexparent] = (*pTopology)[child][1][indexcand];
                    }
                }
            }
        }
    }
    
	int totalNumberOfParents = 0;
    
	for(size_t i = 0; i < parentSubsets.size(); i++)
	{
		totalNumberOfParents += parentSubsets[i].size();
	}
    
    vector<int> delaysCombined(totalNumberOfParents);
    
	bool first = true;		// bool to indicate to generate the first set of delays
	bool last = false;	// bool to indicate last enumerate has reached. changed from true to false by YZ 2010.12.6
    
	do{
		// Generate parent time offsets, until all possible
		// offsets have been gone through, last = true
		// when this has happened
		last = generateDelays(delaysCombined, gc.m_min_Markov_order, gc.m_max_Markov_order, first);
        
		int tempIndex = 0;  //to assign combined delays to individual ones
        
		for(size_t i = 0; i < parentSubsets.size(); i++)
		{
			for(size_t j = 0; j < parentSubsets[i].size(); j++)
			{
				delays[i][j] = delaysCombined[ tempIndex + j ];
			}
			tempIndex += parentSubsets[i].size();
		}
        
		bool flag = true; // whether or not consider this enumeration of delays
        
		if(pTopology)
		{
			for(size_t indexcond = 0; indexcond < delays.size(); indexcond++)
			{
				for(size_t indexparent = 0; indexparent < delays[indexcond].size(); indexparent ++)
				{
					if(delays[indexcond][indexparent]!=delaysFromPathway[indexcond][indexparent]
                       && delaysFromPathway[indexcond][indexparent]<=0)
					{
						flag = false;
						break;
					}
				}
			}
		}
        
        //test whether enumerated parent subsets match individual topology file provided
        
        for(size_t i = 0; i < parentSubsets.size(); i++)
        {
            vector<int> candidates;
            vector<int> parentsForOneCond = parentSubsets[i];
            if(gc.m_pathways.size()>1)
            {
                //enumerated parent sets have to be exactly the same as defined by pathway
                candidates = gc.m_pathways[i].get()[child][0];
                if(intersectionSet(parentsForOneCond, candidates).size()!=candidates.size())
                {
                    flag = false;
                    break;
                }
                
            }else if(gc.m_candidateTopologys.size()>1)
            {
                //enumerated parent sets must be included in pathway definition
                candidates = gc.m_candidateTopologys[i].get()[child][0];
                
                if(unionSet(parentsForOneCond, candidates).size()!=candidates.size())
                {
                    flag = false;
                    break;
                }
            }
        }
        
		if(flag)
		{
			bool updated = processOneEnumeration(child, parentSubsets, delays);
            
			improved = updated ? updated : improved;
		}
		first = false;
        
	} while( ! last );
    
	return improved;
}

bool EnumMethod::enumerateParents(const int child,
                                  const vector<int> & parentsBegin,
                                  const vector<int> & parentsEnd,
                                  const int nNodes,
                                  const GLNConfig & gc,
                                  const vector<int>  & parentList)
{
	bool improved = false;
    
	if(parentsBegin.size()<=0)
	{
		return improved;
	}
    
	int nParents = (int) parentsBegin.size();
	vector<int> parents = parentsBegin;
	
	int position;
    
	do{
        
		int j = nParents;
		vector<int> mapedParents;
		if(parentList.size()!=0)
		{
			mapedParents.resize(parents.size());
			for(size_t i=0;i<mapedParents.size();i++)
			{
				mapedParents[i] = parentList[parents[i]-1];
			}
		}
		if(! gc.allowSelfCycle()) {
            
			for(j=0; j < nParents ; j++) {
				//added by Yang Zhang 7.25.2010 to convert candidate parent ids back to GLN id
				if(parentList.size()!=0)
				{
					if(child == mapedParents[j]-1) {
						break;
					}
				}else
				{
					if(child == parents[j]-1) {
						break;
					}
				}
			}
		}
        
		if( j == nParents ) {
            
			if(parentList.size()!=0) //using candidate parent file to enumerate parents, added by Yang Zhang 7.25.2010
			{
				bool updated = enumerateDelays(child, mapedParents, gc);
				improved = updated ? updated : improved;
			}
			else
			{
				bool updated = enumerateDelays(child, parents, gc);
				improved = updated ? updated : improved;
			}
		}
        
		position = nParents - 1;
        
		while(position >= 0){
			//Check to see if current m_position has reached it's max
			if(parents[position] < (nNodes - (nParents - position - 1)) ){
				break;
			}//end if
			/*moved to end
             if(parents[position]<parentsEnd[position]){
             break;
             }*/
            
			position--;
		}//end while
        
		if(position >= 0){
            
			parents[position]++;
            
			if(position < nParents-1) {
                
				//fill in the rest of the m_indexParents
				for(int i = position+1; i < nParents; i++) {
                    
					parents[i] = parents[i-1]+1;
                    
				} //end for
			}//end if
		}
	} while((position >=0 || parents[0] != nNodes - (nParents - 1))&&parents[position]<=parentsEnd[position] );
    
	return improved;
}

bool EnumMethod::enumerateDelays(const int child,
                                 const vector<int> & parents,
                                 const GLNConfig & gc)
{
    /*
     * With the given set of parents, this function will try ALL the possible delays,
     * i.e. from the minimum Markov order to the possible maximum Makrov order (both given by user)
     * (default of min_Markov_order = max_Markov_order = -1)
     * and set the one with largest p-value as the child's parent set.
     */
	bool improved = false;
    
	vector<int> delays(parents.size());
	vector<int> delaysFromPathway(parents.size());
    
	const vector<vector<vector<int> > > * pTopology = NULL;
    
    
	if(gc.m_pathways.size() == 1)
	{
		pTopology = & (gc.m_pathways[0].get());
        
	}
    else if(gc.m_candidateTopologys.size() == 1)
	{
		pTopology = & (gc.m_candidateTopologys[0].get());
	}
    
	if(pTopology)
	{
		for(size_t indexparent = 0; indexparent < parents.size(); indexparent ++) //index for parents
		{
			for(size_t indexcand = 0; indexcand < (*pTopology)[child][0].size(); indexcand ++ )
			{
				//index for candidates or pathways
				if(parents[indexparent] == (*pTopology)[child][0][indexcand])
				{
					delaysFromPathway[indexparent] = (*pTopology)[child][1][indexcand];
				}
			}
		}
	}
	bool first = true;	// bool to indicate to generate the first set of delays
	bool last;	        // bool to indicate last enumerate has reached.
    
	do{
		// Generate parent time offsets, until all possible
		// offsets have been gone through, last = true
		// when this has happened
        
		// generate_new = false;
        
		last = generateDelays(delays, gc.m_min_Markov_order, gc.m_max_Markov_order, first);
		
		bool flag = true; // whether or not consider this enumeration of delays
		
		if(pTopology) {
			for(size_t index = 0; index < delays.size(); index++)
			{
				if(delays[index] != delaysFromPathway[index]
                   && delaysFromPathway[index]<=0)
				{
					flag = false;
					break;
				}
			}
		}
        
		if(flag)
		{
			bool updated = processOneEnumeration(child, parents, delays);
			improved = updated ? updated : improved;
		}
		first = false;
        
	}while( ! last );
    
	return improved;
}

bool EnumMethod::generateDelays (vector<int> & delays, const int min_Markov_order, const int max_Markov_order, bool first)
{	// Function added by Joe Song on August 7, 2008
    
	bool last = false;
    
	if (first){
        
		for (size_t i = 0; i < delays.size(); ++i) {
			delays[i] = min_Markov_order; // -1; // 8/4/2008 Modified by Joe Song
		}
		// first = false;
        
	} else {
        
		/*
         int flag = 1;  //flag used for max offset reached at array position
         
         for (int i = ((int) delays.size()) - 1 ; i >= 0; i--){
         //if flag is one enter to adjust offset
         if(flag == 1) {
         
         //max offset reached for array position at i
         if(delays[i] <= max_Markov_order){ // 8/4/2008 Modified by Joe Song
         
         delays[i] = min_Markov_order; // -1;  // 8/4/2008 Modified by Joe Song
         flag = 1;
         
         } else { //not the max offset decrement one
         delays[i]--;
         flag = 0;
         break;		// 8/4/2008 Added by Joe Song
         }//end else
         }//end if
         
         }//end for
         */
        
		//this piece of code can only generate same delays, need to consider different delays.
		//e.g. need to change (0,0),(-1,-1) to (0,0),(0,-1),(-1,0),(-1,-1)
		/*
         for (int i = ((int) delays.size()) - 1 ; i >= 0; i--) {
         
         (delays[i] == max_Markov_order) ? delays[i] = min_Markov_order : delays[i] -- ;
         
         if ( delays[i] == min_Markov_order ) break;
         
         } //end for*/
		for (int i = ((int) delays.size()) - 1 ; i >= 0; i--) {
			
			if(delays[i] == max_Markov_order)
			{
				delays[i] = min_Markov_order;
				continue;
			}else
			{
				--delays[i];
				break;
			}
		}
        
        
	}//end else
	
	//check to see if all array positions are at the max delays
	/*
     for (size_t i = 0; i < delays.size(); i++){
     if (delays[i] == max_Markov_order) // 8/4/2008 Modified by Joe Song
     last = true;
     else{
     last = false;
     break;
     }//end else
     }//end for
     */
    
	size_t i;
	for ( i = 0; i < delays.size(); ++i ) {
		if (delays[i] != max_Markov_order) break;
	}
    
	last = ( i == delays.size() ) ? true : false;
    
	return last;
}


void findCommonParents(const vector<vector< int > > & parents,
                       const vector<vector< int > > & delays,
                       vector<int> & intersectparents,
                       vector<int> & intersectdelays)
{
    vector<int> intersectionSet(const vector<int> & firstSet, const vector<int> & secondSet);
    
    // ATTENTION: This portion of parent intersection code does not
    //   handle the same parent id with different delays.  Needs to
    //   be fixed.  10/6/2011.  Joe Song
    //
    // A function that check both parents and delays called findCommonParentsWithDelays()
    //   is created.  MS. 11/21/2011.
    
    intersectparents = parents[0];
    
    for(size_t k=1; k<parents.size(); k++) {
        intersectparents = intersectionSet(intersectparents, parents[k]);
    }
    
    intersectdelays.resize( intersectparents.size() );
    
    for(size_t i=0; i<intersectparents.size(); i++) {
        for(size_t j=0; j<parents[0].size(); j++) {
            if(intersectparents[i] == parents[0][j]) {
                intersectdelays[i] = delays[0][j]; // intersectdelays.push_back(delays[0][j]);
            }
        }
    }
}

vector<Parent> intersectParents(const vector<vector<Parent> > & sets)
{
    vector<Parent> shared = sets[0];
    
    for(size_t i=1; i < sets.size(); i++) {
        
        vector<Parent> result( min(shared.size(), sets[i].size()) );
        
        set<Parent> s1( shared.begin(), shared.end() ); // use "set" to sort the data, required by set_intersection()
        
        set<Parent> s2( sets[i].begin(), sets[i].end() ); // use "set" to sort the data, required by set_intersection()
        
        vector<Parent>::iterator itr = set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), result.begin());
        
        shared = vector<Parent>(result.begin(), itr);
    }
    
    return shared;
}

void findCommonParentsWithSameDelays(const vector<vector< int > > & parents,
                                     const vector<vector< int > > & delays,
                                     vector<int> & commonParents,
                                     vector<int> & commonDelays)
{
    // This function checks both parents and their delays.  Two parents are considered common
    //   if and only if they have the same id and delay.  MS. 11/21/2011.
    
    assert(parents.size() == delays.size());
    
    vector< vector<Parent> > sets( parents.size() ); 
    
    for (size_t i=0; i < sets.size(); i++) {
        sets[i].resize( parents[i].size() );
        for (size_t j=0; j < sets[i].size(); j++) {
            sets[i][j] = Parent( parents[i][j], delays[i][j] );
        }
    }
    
    vector<Parent> shared = intersectParents(sets);
    
    commonParents.resize( shared.size() );
    commonDelays.resize( shared.size() );
    
    for (size_t j=0; j < shared.size(); j++) {
        commonParents[j] = shared[j].getId();
        commonDelays[j] = shared[j].getDelay();
    }
}

void mergeParentsWithSameDelays(const vector<vector< int > > & parents,
                                const vector<vector< int > > & delays,
                                vector<int> & parentSuperset,
                                vector<int> & delaySuperset)
{
    // This function merges all parent sets into a superset. It will eliminate
    //   duplicate parents with same delays.  MS 6/2/2013.
    
    assert(parents.size() == delays.size());
    size_t total = 0;
    
    for (size_t i=0; i < parents.size(); i++) {
        total += parents[i].size();
    }
    
    vector<Parent> superList( total );
    
    for (size_t i=0, k=0; i < parents.size(); i++) {
        for (size_t j=0; j < parents[i].size(); j++, k++) {
            superList[k] = Parent( parents[i][j], delays[i][j] );
        }
    }
    
    set<Parent> uniqueSet;
    uniqueSet.insert(superList.begin(), superList.end());
    
    // unique(superSet.begin(), superSet.end());
    
    parentSuperset.resize( uniqueSet.size() );
    delaySuperset.resize( uniqueSet.size() );
    
    set<Parent>::iterator it = uniqueSet.begin();
    for (size_t j=0; j < uniqueSet.size(); ++j, ++it) {
        parentSuperset[j] = (*it).getId();
        delaySuperset[j] = (*it).getDelay();
    }
    assert(total >= uniqueSet.size());
}


bool operator < (const Parent & p1, const Parent & p2)
{
    return (p1.m_delay < p2.m_delay)
    || (p1.m_delay == p2.m_delay && p1.m_id < p2.m_id);
}

