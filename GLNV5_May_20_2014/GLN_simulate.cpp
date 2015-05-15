// GLN_simulate.cpp
//
// Joe Song
//
// Created: August 8, 2006
//
// Modified: 
//
//   September 7, 2008
// 
//   April 23, 2011.  Add a more precise initialization of srand() using
//    microseconds elapsed since Jan 1, 1970, but it only works for linux or MacOS, 
//    but not windows.  It addresses the problem of trajectory simulation with random noise,
//    when multiple runs of simulate() function is called too quickly next to each other,
//    resulting the same random seed if we use time(NULL) to get seconds elapsed since
//    1/1/1970. This problem is more likely to occur in the future when the computer 
//    runs at a very high clock frequency.

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::bad_alloc;

#include <fstream>
using std::ofstream;

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <algorithm>
#include <set>

using std::set;

#include <cmath>
using std::min;

#include "GLN.h"
#include "GLNGlobals.h"

list< vector<int> > GeneralizedLogicalNetwork::computeStateTrans()
{ // Works only for Markovian order of 1.
    
	long long unsigned numb_of_states=1;
    long long unsigned num_to_index;
    
	// int start;
	// int current;
	
	vector<int> edge(2);				//The array that will be pushed onto the list
	list< vector<int> > edgeList;			//The list edges

	vector<int> bases(size()); //Declaration for array of bases

	vector<int> old_values(size()); //Declaration for array of old current values

	vector<int> state(size()); //Declaration for array of current states

	/********************************************************************************
	******Find the number of total states and fill in the array of bases*************
	*Also keeps original values inside of the nodes before we modify them
	*********************************************************************************/
	for(size_t i=0; i < size(); i++){
		old_values[i] = m_nodes[i].getValue();		//Keep all the old values before simulation
		bases[i] = m_nodes[i].getBase();		//Get the base and store into array
		numb_of_states = numb_of_states * bases[i];	//Multiply all the bases for the total states
	}//end for

	cout << "Total possible states in this network: " << numb_of_states << endl;

	for(long long unsigned i=0; i < numb_of_states; i++){

		/******************************************************
		******Convert i to a trajectory state in an array******
		*******************************************************/
		
        LinearIndexToArrayIndex(bases, i, state);

        /*
        current = i;
		for(int j = (int)(size()) -1; j >=0; j--){
			state[j] = current % bases[j];
			current = current / bases[j];
		}//end for
        */
        
		/**************************************************
		******Change the current values of the nodes ******
		***************************************************/
		for(size_t j = 0; j < size(); j++){
			m_nodes[j].setValue(state[j]);
		}//end for

		/*************************************************************
		******Do predictor on each array index and replace value******
		**************************************************************/

        /*
         for(size_t j = 0; j < size(); j++){
         state[j] = m_nodes[j].predict(trajectory);
         }//end for
         */
        
        vector< vector<int> > trajectory(1, state);        
		state = simulate(trajectory, 1);

		/*********************************************************
		******Convert the state(array) into a regular number******
		**********************************************************/
        /*
        num_to_index = 0;
        for(int j=0; j < ((int)size())-1; j++){
			start = j+1;
			num_to_index = (num_to_index + state[j]) * bases[start];		 
		}//end for

		num_to_index = num_to_index + state[size()-1];	//add the last digit in the array
        */
        
        num_to_index = ArrayIndexToLinearIndex(bases, state);
        
		/***********************************************************
		******Store the edges and push the pair onto the list*******
		************************************************************/
		edge[0]= static_cast<int>(i);
		edge[1]= static_cast<int>(num_to_index);
		edgeList.push_back(edge);

	}//end for

	/**************************************************
	************Restore old original values************
	***************************************************/
	for(size_t i = 0; i < size(); i++){
		m_nodes[i].setValue(old_values[i]);
	}//end for

	return edgeList;
}//end computeStateTrans

void GeneralizedLogicalNetwork::generateRandInitTrajCol(TrajectoryCollection & trjcollection, int numberoftrj, int timesteps) const
{
	// NOTE: (MS Nov 24, 2011) 
    //       Random values are no longer generated for an external node at all time points
	//       except the initial time.  The random value of a external node set at t=0 will 
    //       be repeated for t>0 for that node.
    
    int nNodes = size();
    
	trjcollection.allocate(numberoftrj, timesteps, nNodes, -1);

	// Maximum markovian order,to make sure we have enough timepoints to start simulate
	// We use negative number.
	// We must generate at least one timestep for each trajectory.
	// int maxoffset = min(-1, getMaxMarkovOrder());

	for(int i=0; i<numberoftrj; i++)
	{
		for(int j=0; j<timesteps; j++)
		{
			for(int k=0; k<nNodes; k++)
			{
				// if((m_nodes[k].getType()=='i' && j<-maxoffset) || m_nodes[k].getType()=='e')
				// {
				//		trjcollection.set( i, j, (int) k, rand() % m_nodes[k].getBase() );
				// }
				// Random values are generated for all nodes
                
                switch(m_nodes[k].getType()) {
                        
                    case 'i': 
                        trjcollection.set( i, j, k, rand() % m_nodes[k].getBase() );
                        break;
                        
                    case 'e': // MS Nov 24, 2011 added this case to handle external node 
                        if (j == 0) { // randomly set the initial value at t=0 for the external node
                            trjcollection.set( i, j, k, rand() % m_nodes[k].getBase() );                 
                        } else { // repeat the initial value for the external node
                            trjcollection.set( i, j, k, 
                                              trjcollection.get( i, 0, k ) );
                        }
                        break;
                        
                    default:
                        GLNExit(EXIT_FAILURE, "ERROR: Invalid node type!");
                        break;
                }
			}
		}
	}
}

vector<int> GeneralizedLogicalNetwork::simulate(const vector< vector<int> > & traj, int t) const
{
	vector<int> state_t(size());

	vector<int> parentValues;

	for(int i=0; i<(int)size(); i++) {

		if(m_nodes[i].getType() == 'e') {
		
			if((int) traj.size() > t) {

				// *** THIS REQUIRES ALL EXTERNAL NODES VALUES BE PROVIDED
				state_t[i] = traj[t][i];

			} else {

                state_t[i] = traj[t-1][i]; // MS Nov 27, 2011

                /* MS Nov 27, 2011
				cerr << "ERROR: Trajectory not initialized properly for simulation!" << endl;
				GLNExit(EXIT_FAILURE);
                */
			}

		} else if(m_nodes[i].getType() == 'i') {

			// get parent values at the right offsets
			const vector<int> & parents = m_gtts[i].getParents();
			const vector<int> & delays = m_gtts[i].getDelays();
			parentValues.resize(parents.size());

			size_t k;

			for(k=0; k < parents.size(); k++) {

				if(t + delays[k] >= 0) {
					// *** THIS REQUIRES ALL PARENT NODES VALUES WITH THE DELAYS BE PROVIDED
					parentValues[k] = traj.at(t + delays[k])[ parents[k] - 1 ];
				} else {
					break;
				}
			}

			if(k == parents.size()&& k!=0) {
			
				// Could inconsistent when some parents have zero order.
				state_t[i] = m_gtts[i].lookup(parentValues, traj[t-1][i]);

			} else {

				if(traj.size() > (size_t) t) {

					state_t[i] = traj[t][i];

				} else {

                    state_t[i] = traj[t-1][i]; // MS Nov 27, 2011

                    /*
					cerr << "ERROR: Trajectory not initialized properly for simulation!" << endl;
					GLNExit(EXIT_FAILURE);
                    */
				}
			}
		}
	}

	return state_t;
}

void GeneralizedLogicalNetwork::simulate(TrajectoryCollection & trjCol) const
{
	for(int i = 0; i < (int) trjCol.getNumTrajectories(); i++)	{

		for(int t = 0; t<(trjCol.getTrajectoryLength((int)i)); t++) {

			vector<int> state = simulate(trjCol.getTrajectory((int) i), t);
			trjCol.setState(i, t, state);

		}
	}
}

void GeneralizedLogicalNetwork::simulateExhaustively(TrajectoryCollection & trjCol) const
{  // Works only for Markovian order of 1.
    
    long long unsigned nStates=1;
    
    vector< vector< vector< int > > > trjs;
    
	vector<int> bases(size()); // bases of each node
	vector<int> ai(size()); // array index
    
	for(size_t i=0; i < size(); i++){
		bases[i] = m_nodes[i].getBase();    // Get the base and store into array
		nStates = nStates * bases[i];       // Multiply all the bases for the total states
	}//end for

    set<long long unsigned> visited;
    
	for(long long unsigned li=0; li < nStates; li++){

        // Generate next initial state not previously visited
        
        if(visited.find(li) != visited.end()) { // check li
            continue;
        }
        
        LinearIndexToArrayIndex(bases, li, ai);

        vector< vector<int> > trajectory(1, ai);
        
        long long unsigned currLi;
        
        do {
            // Simulate from the initial state
            //   until a cycle or steady state has been reached    
            
            ai = simulate(trajectory, trajectory.size());
            
            trajectory.push_back(ai);
            
            currLi = ArrayIndexToLinearIndex(bases, ai);

            // check currLi
            if(currLi <= li || visited.find(currLi) != visited.end()) {
                break;
            } else {
                // insert to set if not already visited
                visited.insert(currLi);
            }

        } while (true);
        
        trjs.push_back(trajectory);

	} //end for

    trjCol.setTrajTables(trjs);
}

//added by Yang Zhang 7/22/2012, add a new -C parameter to specify the new house noise model
//the new model can achieve pure random noise when noise level is 1 and do not need to separate quantization level 2 from others
void GeneralizedLogicalNetwork::simulate(const string & fileInitTrj, 
                                         const string & fileOutputTrj, 
										 int nTrajs, int length, double levelNoise, const string & houseModel) const
{
	TrajectoryCollection trjCol;

	if(fileInitTrj.empty())	{
        
        trjCol.setVersion("TRAJECTORY_VER2");
        trjCol.setBases(getBases());
        trjCol.setn_nodes(size());

        trjCol.setIntNodeNames(getNodeNames('i'));
        trjCol.setExtNodeNames(getNodeNames('e'));
                
	}

    if( (nTrajs > 0 && length > 0 ) || ! fileInitTrj.empty()) { 

        if(fileInitTrj.empty()) {
            // Initialize the TrajectoryCollection object randomly if no 
            //   initial trajectory collection file is provided
            generateRandInitTrajCol(trjCol, nTrajs, length);
            
        } else {
            // read initial trajectory file
            trjCol.scan(fileInitTrj);
        }
        
        simulate(trjCol);
        
    } else { // perform exaustive simulation if nothing is specified
        
        simulateExhaustively(trjCol);

    }
    
	if(levelNoise > 0) trjCol.addNoise(levelNoise, houseModel);

	trjCol.save(fileOutputTrj.c_str());
}



//Haizhou Wang, Oct 11, 2012
const TrajectoryCollection GeneralizedLogicalNetwork::simulate(const string & fileInitTrj, 
										 int nTrajs, int length, double levelNoise, const string & houseModel) const
{
	TrajectoryCollection trjCol;

	if(fileInitTrj.empty())	{
        
        trjCol.setVersion("TRAJECTORY_VER2");
        trjCol.setBases(getBases());
        trjCol.setn_nodes(size());

        trjCol.setIntNodeNames(getNodeNames('i'));
        trjCol.setExtNodeNames(getNodeNames('e'));
                
	}

    if( (nTrajs > 0 && length > 0 ) || ! fileInitTrj.empty()) { 

        if(fileInitTrj.empty()) {
            // Initialize the TrajectoryCollection object randomly if no 
            //   initial trajectory collection file is provided
            generateRandInitTrajCol(trjCol, nTrajs, length);
            
        } else {
            // read initial trajectory file
            trjCol.scan(fileInitTrj);
        }
        
        simulate(trjCol);
        
    } else { // perform exaustive simulation if nothing is specified
        
        simulateExhaustively(trjCol);

    }
    
	if(levelNoise > 0) trjCol.addNoise(levelNoise, houseModel);

	return trjCol;
}
