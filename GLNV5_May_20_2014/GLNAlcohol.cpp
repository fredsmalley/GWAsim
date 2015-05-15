//
//  GLNAlcohol.cpp -- Extracted from GLN_simulate.cpp on 11/27/2011
//  gln
//
//  Created by Joe Song on 11/27/11.
//  Copyright (c) 2011 New Mexico State University. All rights reserved.
//

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

#include "GLN.h"

using std::cout;
using std::endl;

void outputStatesDotFormat(const list< vector<int> > & lst, const string & filename);

void GeneralizedLogicalNetwork::detectCycleAcohol(vector< vector<int> > & trajectory, 
												  const vector<int> & bases, 
												  list< vector< int > > & edgeList,
												  vector<int> & vecStates, int prev_li)
{
	cout << "Start traversing ... "; 
    
	int li;
	// int prev_li;
    
	while(true) { // Simulation with alcohol on:
        
		vector<int> edge(2);    // The array that will be pushed onto the list
		vector<int> ai(size()); // array index of current state
        
		// simulate one step
		/*
         for(size_t j = 0; j < size(); j++){
         ai[j] = m_nodes[j].predict(trajectory);
         if(j > 0) {
         cout << ai[j];
         }
         } //end for
         cout << endl;
         */
		ai = simulate(trajectory, (int) trajectory.size());
        
		ai[0] = 0; // quick hack -- external node should be set according to the input
        
		// Convert array index to linear index
		/*********************************************************
         ******Convert the array index of state into a linear index of a state
         **********************************************************/
		li = ArrayIndexToLinearIndex(bases, ai);
        
		/***********************************************************
         ******Store the edges and push the pair onto the list*******
         ************************************************************/
		edge[0]= prev_li;
		edge[1]= li;
		edgeList.push_back(edge);
		
		cout << edge[0] << " -> " << edge[1] << endl;
        
		// mark edge in first color
		// ...
        
		// check if a cycle has resulted with previous states
        
		vector<int>::iterator itr = find(vecStates.begin(), vecStates.end(), li);
		if(itr == vecStates.end()) {
			vecStates.push_back(li);
			/**************************************************
             ******Change the current values of the nodes ******
             ***************************************************/
			for(size_t j = 0; j < size(); j++){
				m_nodes[j].setValue(ai[j]);
			}//end for
            
			prev_li = li;
            
		} else {
			// if so stop
			break;
		}
	}
}

void GeneralizedLogicalNetwork::simulateAcohol(const char * fileAcoNet, const char * fileDOT) 
{
	list< vector< int > > edgeList;
	list< string > listColors;
	vector<int> vecStates;
    
	vector< vector<int> > trajectory;
    
	vector<int> bases(getBases());				//Declaration for array of bases
	vector<int> old_values;			//Declaration for array of old current values
    
	scan(fileAcoNet);
    
	for(size_t i=0; i < size(); i++){
		old_values.push_back(m_nodes[i].getValue());		//Keep all the old values before simulation
		// bases.push_back(m_nodes[i].getBase());		//Get the base and store into array
		if(i > 0) {
			cout << old_values[i];
		}
	} //end for
	cout << endl;
    
	old_values[0]=0; // The array index of the states does NOT include the alcohol node!
	int li = ArrayIndexToLinearIndex(bases, old_values);
	vecStates.push_back(li);
    
	int prev_li = li;	// linear index of current state
    
	m_nodes[0].setValue(1); // alcohol: 1 0 0 0 ...
    
    // #include "DetectCycle.cpp"
	detectCycleAcohol(trajectory, bases, edgeList, vecStates, prev_li);
    
	// restore the initial state
    
	/**************************************************
     ************Restore old original values************
     ***************************************************/
	for(size_t i = 0; i < size(); i++){
		m_nodes[i].setValue(old_values[i]);
		if(i > 0) 
			cout << old_values[i];
	}//end for
	cout << endl;
    
	// change the first node (alcohol) initial value to 0, which was 1 originally
	old_values[0] = 0;
	prev_li = ArrayIndexToLinearIndex(bases, old_values);
    
	m_nodes[0].setValue(0); // alcohol: 0 0 0 0 ...
    
	// vecStates.clear();
    
    // #include "DetectCycle.cpp"
	detectCycleAcohol(trajectory, bases, edgeList, vecStates, prev_li);
    
	// External input 2:
	// simulate one step
	// mark edge in second color
	// check if a cycle has resulted with previous state
	// if so stop
    
	// output the state transition diagram
	outputStatesDotFormat(edgeList, fileDOT);
}
