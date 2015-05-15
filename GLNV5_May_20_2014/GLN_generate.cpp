// GLN_generate.cpp
//
// Curtis Luce
// Created: 2007
// Modified: September 7, 2008.  Joe Song
// Last modified: September 7, 2008.  Joe Song. Named changed from generateRandomNetwork.cpp

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::bad_alloc;

#include <sstream>
using std::stringstream;

#include <fstream>
using std::ofstream;

#include <vector>
using std::vector;

#include <string>
using std::to_string;

#include <ctime>
#include <cstdlib>
#include "GLN.h"

int randomNumberGenerator(int minNumber, int maxNumber)
{
	int random;

	int N = maxNumber-minNumber+1;

	random = (int) (N * (rand() / ((double)RAND_MAX+1.0)));

	return random+minNumber;
} //end randomNumberGenerator

//Haizhou Wang, Feb 6, 2013
//GeneralizedLogicalNetwork::generateRandomNetwork() may generate fictitious parent for a node
void GeneralizedLogicalNetwork::generateRandomNetwork
(int networkSize, int minBase, int maxBase,int minParents, int maxParents, 
 int min_Markov_order, int max_Markov_order, string init_traj_file, 
 int numb_of_init)
{
	ofstream out(init_traj_file.c_str());
	vector<int> base_lookup(networkSize);			//Base table
	int initial_val;			//Holds initial value
	int numb_of_parents;			//Holds the number of parents
	int tt_size;				//Truth table size
	int temp, temp2;			//Temporary variables
	int base, tt_val;
	vector<int> parents;			//Holds all the parents
	vector<int> parentBases;	//holds parent bases
	vector<int> truth_table;
	// Predictor pred;
	int node_id = 1;			//start node id at 1
	vector<int> parent_offsets;		//holds parent time offsets
	int offset;				//offset read in

    vector<string> parentNames;
	
	// std::srand( (unsigned)time( NULL ) );	//Seed for srand
	initialize(networkSize);		//Initialze network size	

	if((numb_of_init != -1) && (init_traj_file.size() != 0)){
		//ofstream out(init_traj_file.c_str());
		out << "TRAJECTORY_VER2" << endl;
		out << numb_of_init << " " << networkSize << " " << 0 << endl;
	}
	
	cout << "Network size: " << networkSize;
	cout << ", (min, max) Base: " << minBase << "," << maxBase;
	cout << ", Max #Parents: " << maxParents << endl;
	//add maximum parent offset
	cout << "Bases of each node: ";

	//Generate the entire base 
	for(int i=0; i < networkSize; i++){
		base_lookup[i] = randomNumberGenerator(minBase, maxBase);		//no base 0! Base 1 and up
		cout << base_lookup[i] << " ";
		if((numb_of_init != -1) && (init_traj_file.size() != 0))
			out << base_lookup[i] << " ";
	}//end for
	cout << endl;
	//print out default node names
	if((numb_of_init != -1) && (init_traj_file.size() != 0)){
		out << endl;
		for(int k = 0; k < networkSize; k++)
			out << "Node" << k+1 << " ";
		out << "\n\n";
	}//end if
	
	for(int i=0; i < networkSize; i++){

		//Reset variables
		tt_size=1;
		temp=1;
		temp2=1;

		initial_val = randomNumberGenerator(0, base_lookup[i]-1);		//Get a initial value between 0 and base-1
		cout << "[Node " << i+1 << "] >> Initial value: " << initial_val << ", ";
		cout << "Base: " << base_lookup[i] << ", ";
		//changed minParents to 1 by YangZhang 12052008
		//changed 1 to minParents by Yang Zhang 11/9/2009
		numb_of_parents = randomNumberGenerator(minParents, maxParents);		//Get the number of parents. From 0 to maxParents
		// cout << "#parents: " << numb_of_parents << ", ";

		//Get random parents based on number of parents, find from 0 to networkSize-1
		cout << "Parents IDs: ";

		/****************************************************
		********Generate Parents and Truth Table size********
		*****************************************************/

		for(int j=0; j < numb_of_parents; j++){

			//Needs to go back and see if any other nodes have this parent
			//parent starts from 1 modified by YangZhang12052008
			//temp = randomNumberGenerator(0, networkSize-1);
			temp = randomNumberGenerator(1, networkSize);
			for(int k=0; k < j; k++){	
				if(temp == parents[k]){
					//We already have this parent -- find another parent
					k = -1;			//Reset k
					//temp = randomNumberGenerator(0, networkSize-1);
					temp = randomNumberGenerator(1, networkSize);
				}//end if

			}//end for
			cout << temp+1 << " ";

			parents.push_back(temp);
			//base starts from 0 modified by YangZhang12052008
			//temp2 = base_lookup[temp];

            parentNames.push_back( "Node" + to_string( temp ) );

			temp2 = base_lookup[temp - 1];
			//cout << "[" << temp2 <<"] ";
			parentBases.push_back(temp2);				//get the base of the parents also
			tt_size = tt_size * temp2;				//figure out the size of truth table
		}//end for


		cout << endl;
		
		//no parents then push zero as the time offsets
		if (numb_of_parents == 0)
			for (int u = 0; u < numb_of_parents; u++)
				parent_offsets.push_back(0);
				
		else{
			cout << "Parent time offsets:  ";
			for (int p = 0; p < numb_of_parents; p++){
				//generate offsets for each parent that has been generated
				offset = randomNumberGenerator(min_Markov_order, max_Markov_order);
				cout << offset << " ";
				parent_offsets.push_back(offset);
			}//end for
			cout << endl;
		}//end else
		
		cout << "Truth Table size: " << tt_size << endl;
		//Generate numbers for the Truth Table
		cout << "Truth Table: ";

		/*************************************
		********Truth Table generation********
		**************************************/
		if(numb_of_parents == 0){
			truth_table.push_back(initial_val);	//If no parents, truth table with just initial value
		}//end if	

		else{
			for(int j=0; j < tt_size; j++){
				base = base_lookup[i]-1;
				tt_val = randomNumberGenerator(0, base);
				cout << tt_val << ",";
				truth_table.push_back(tt_val);

			}//end if
		}//end else
		//create a default name for the node, 	
		
		//char node_name[50];
		//sprintf(node_name, "Node%d", node_id);
		
		stringstream intStream;
		intStream << node_id;
		string nodeName = "Node" + intStream.str();

		GeneralizedTruthTable gtt(i, base_lookup[i], parents, parentBases, parent_offsets, truth_table); 

		m_gtts[i] = gtt;
        m_gtts[i].setParentNames( parentNames );

		initializeNode(i, initial_val, base_lookup[i], node_id, nodeName.c_str(), 'i');

		cout << endl;

		node_id++;
		truth_table.clear();
		parentBases.clear();
		parents.clear();
        parentNames.clear();


	}//end for

	//prints out a initial trajectory file
	if((numb_of_init != -1) && (init_traj_file.size() != 0)){
		int traj = 0;  //number of trajectories
		for (int i = max_Markov_order; i <= min_Markov_order; i++)
			traj++;
		for (int k = 0; k < numb_of_init; k++){
			out << traj << endl;  //print out number of trajectories
			for(int i = max_Markov_order; i <= min_Markov_order; i++){
				for (int l = 0; l < networkSize; l++){
					//output bases for the parents
					out << randomNumberGenerator (0, m_nodes[l].getBase()-1) << " ";
				}//end for
			out << endl;
			}//end for
			out << "\n\n";
		}//end for
	}//end if
	
	out.close();

}//end generateRandomNetwork

void GeneralizedLogicalNetwork::initialStates(string infile, string outfile, int numTraj, int min_Markov_order, int max_Markov_order){
	ofstream out(outfile.c_str()); //temp output
	ifstream in(infile.c_str());	//input stream
	string buffer;
	int base;
	int numNodes;
	char property;
	vector<int> parent_bases;  //holds parent bases
	vector<string> nodeNames;  //holds node names
	vector<char> nodeProp;   //holds node properties
	int extCount = 0;  //external node count
	int intCount = 0;  //internal node count
	if(!in){
		cout << "Cannot open file, make sure the file is in the same directory as the program \n";
	}//end if
	
	//check for logical network header, if not the right file then do nothing
	in >> buffer;
	if ((buffer != "LOGICAL_NETWORK_VER1")  && (buffer != "LOGICAL_NETWORK_VER2"))
		cout << "This is not a logical network file.  Choose a different file."  << endl;
		
	else{
		out << "TRAJECTORY_VER2" << endl;
	
		out << numTraj << " ";	
	
		in >> numNodes;
		cout << "numNodes: " << numNodes << endl;
		for (int i = 0; i < numNodes; i++){
			in >> base;	//node id
	
			in >> buffer;	//node name
			nodeNames.push_back(buffer);
		
			in >> base;	//node base		
			parent_bases.push_back(base);
		
			in >> base;  	//initial value
		
			in >> property;	//property

			if (property == 'e')
				extCount++;  //external node count
			else
				intCount++;  //internal node count
			nodeProp.push_back(property);
		
			in >> base;	//number of parents
		
			//no parents do nothing.
			if (base == 0)
				in >> property;		//increment past end of node
			
			if (base > 0){
				//read in parents
				for (int i = 0; i < base; i++)
					in >> base;
				//read in time offsets
				for (int i = 0; i < base; i++)
					in >> base;
				//read in parents bases
				for (int i = 0; i < base; i++)
					in >> base;
			
				while (property != 42){
					in >> property;
				}//end while
			}//end if
		}//end for
		
		out << intCount << " " << extCount << endl;
		
		for (int i = 0; i < numNodes; i++)
			out << parent_bases[i] << " ";
		out << endl;
		
		for (int i = 0; i < numNodes; i++)
			if(nodeProp[i] == 'i')
				out << nodeNames[i] << " ";
			
		out << endl;
		
		for (int i = 0; i < numNodes; i++)
			if(nodeProp[i] == 'e')
				out << nodeNames[i] << " ";
			
		out << endl; //spacing
		out << endl; //spacing
		int traj = 0;  //number of initial trajectories
		//number of initial trajectories are decided by the max offset
		for (int i = max_Markov_order; i <= min_Markov_order; i++)
			traj++;
		for (int k = 0; k < numTraj; k++){
			out << traj << endl;
			for(int i = max_Markov_order; i <= min_Markov_order; i++){
				for (int l = 0; l < numNodes; l++){
					//print out parent bases for parents
					out << randomNumberGenerator (0, parent_bases[l]-1) << " ";
				}  //end for
			out << endl;
			}  //end for
			out << "\n\n";
		}  //end for
	}  //end else
}//end generate_random_states
