// GLN_IO.cpp
//
// Curtis Luce, Joe Song
//
// Created: 2007
// Modified: September 7, 2008.  Joe Song
// Last modified: September 27, 2008.  Joe Song. Name changed form discrete_network_IO.cpp

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::ofstream;

#include <cstdlib>

#include "GLNGlobals.h"
#include "GLN.h"

void GeneralizedLogicalNetwork::print() const
{
    
    
}


void GeneralizedLogicalNetwork::save(const string & filename)
/***************************************************************
 **Outputs data of the network into a file selected by the user**
 **The format of the output file will be:
 
 *Total number of nodes in the network
 *Node id (NOTE: Although printed as "Node 1", it is actually Node 0 when inputed into the network)
 *Node name
 *The base of the node
 *The initial value of the Node
 *Node properties
 *Number of parents
 *List of parents
 *Time offset of parent nodes
 *Base of parents
 *Converted Truth Table
 ****************************************************************/
{
    if(filename.empty()) return;
    
	//network needs to be initialized first -- check somehow
	ofstream out(filename.c_str()); //temp output

	// The header is required for scanning later
	// if(ver == "TRAJECTORY_VER1")
	// out << "GENERALIZED_LOGICAL_NETWORK_VER1" << endl;
	out << "GENERALIZED_LOGICAL_NETWORK_VER2" << endl;
    
    //Comment(Haizhou): It's better to use:  out<<m_version<<endl; here

	out << size() << "\n";	//Print out size of network

	for(size_t i=0; i < size(); i++){

		out << i+1 << endl;	//Print Node id
		
		out << m_nodes[i].getName() << endl;	//print out node name
		out << m_nodes[i].getBase() << endl;	//Get the base 
		out << m_nodes[i].getValue() << endl;	//Get initial value
		out << m_nodes[i].getType() <<  endl;	//print out node properties e or i

		int nParents = (int) m_gtts[i].getParents().size();	//get the number of parents for node

		if(nParents == 0){

			out << "0" << endl;	//no parents

		} else {

			out << nParents << endl;		//print out the number of parents

			for(int j=0; j < nParents; j++){
				out << m_gtts[i].getParents()[j] << " ";  // removed +1
			}//end for
			out << endl;
				
			vector <int> offset = m_gtts[i].getDelays();
			for(int j=0; j < nParents; j++){
				out << offset[j]<< " ";
			}//end for
			out << endl;

			//print the base of it's parents
			for(int j=0; j < nParents; j++){
				out << m_gtts[i].getParentBases()[j] << " ";
			}//end for

			out << endl;

			//output the truth table on one line (Joe Song 10/2/2011)
			for(size_t j=0; j < m_gtts[i].getTruthValues().size(); j++){
				out << m_gtts[i].getTruthValues()[j] << ' ';
			}//end for
            out << endl;
            
			// Output the transition probability table
			const vector< vector<double> > & ptt = m_gtts[i].getTransProbTable();

			for(size_t j=0; j < ptt.size(); j++) {
				for(size_t c=0; c < ptt[j].size(); c++) {
					out << ptt[j][c] << " ";
				}
				out << endl;
			}//end for

		}//end else

		out << "*" << endl;	//denotes end of node
	}//for
	out << "done" << endl;	//printed out to show end of file

	out.close();
}

void GeneralizedLogicalNetwork::scan(const string & filename)
// Reads in a GLN format file to initialize a network
{
	setMaxMarkovOrder(0);
	string strBuffer;
	int count = 0;

	//Variables that temporarily hold in the read data, reused for nodes
	int total_nodes;			//Holds the total nodes in the network

	int base, par_base;			//Holds base, holds base of current parent
	int number_of_parents;			//The number of parents the current node has
	int init_val;				//Initial value of the node
	int value;				//Temporarily holds truth table value before inserting into truth table
	int tt_rows;				//Number of Truth Table rows
	int time_offset;			//holds time of parent in time offset
	unsigned short parent;			//Temporarily holds parent value before inserting into parent_vector
	// Predictor pred;				//Predictor function that gets manipulated
	vector<int> parent_vector;		//Holds all the parents for each node
	vector<int> radixes;		//Hold's base of all the parents

	vector<int> truth_table;	//Truth table - one dimensional
	
	vector< vector<double> > prob_truth_table; 

	vector<int> parent_time;	//holds the time offset of the parents
	unsigned char node_type;
	string node_name;
	int node_id;

    vector<string> nodeNameList; //Store the names of each node
    vector< vector< int> > nodeParentsIDs; //Store the names of parents for each node
	
	ifstream in(filename.c_str());	//input stream

	if(!in.is_open()){
		cerr << "ERROR: Cannot open file \"" << filename << "\"!" << endl;
		GLNExit(EXIT_FAILURE);
	}//end if

	/**********************************File Format**************************************
	The file has no words in the file and the file format must be followed with no
	extra words or characters in the file, no error checking can be done except for the 
	header at the top the code depends on the file being in the right format
	**********************************************************************************/

	in >> strBuffer;
	setVersion(strBuffer);

	/** Give option of inputting another file **/
	if (strBuffer != "GENERALIZED_LOGICAL_NETWORK_VER1" 
		&& strBuffer != "GENERALIZED_LOGICAL_NETWORK_VER2" 
		&& strBuffer != "LOGICAL_NETWORK_VER1" 
		&& strBuffer != "LOGICAL_NETWORK_VER2" ){
		cerr << "ERROR: The GLN file version \"" << strBuffer 
			<< "\" is not supported by this program!" << endl;
		GLNExit(EXIT_FAILURE);
	}//end if
	
	else{
		
		in >> total_nodes;		//Read in the total nodes for the network
        nodeParentsIDs.resize( total_nodes );
		initialize(total_nodes);	//send total nodes
		while(count < total_nodes){
			in >> node_id;  /***increment past node id****/
			
			in >> node_name;  //takes in node name, increment past
			nodeNameList.push_back( node_name );
			//in base case
			in >> base;

			//initial value case
			in >> init_val;
			
			in >> node_type;	//takes in node property, increment past	
			
			//number of parents case
			in >> number_of_parents;

			//list of parents case
			if(number_of_parents == 0){
				//Do nothing since there is no parent to read in
			}//end if
	
			else{
				//loop to insert into vector
				for(int i=0; i < number_of_parents; i++){

					in >> parent;
					parent_vector.push_back(parent);  // removed -1
				}//end for
			}//end else

			if(number_of_parents == 0){
				//do nothing no time offsets for parents
			}//end if
			
			else{
				for(int i = 0; i < number_of_parents; i++){
					in >> time_offset;
					if(time_offset < getMaxMarkovOrder())
						setMaxMarkovOrder(time_offset);
					parent_time.push_back(time_offset);
				}//end for
			}//end else
				
			
			//base case
			if(number_of_parents == 0){
				radixes.push_back(base);	//If no parents, push it's own base on there
			}//end if

			else{
				for(int i=0; i <number_of_parents; i++){
						
					in >> par_base;
					//cout << "parent base:  " << par_base << endl;
					radixes.push_back(par_base);
				}//end for
			}//end else

			//truth table case
			//0 parent case
			if(number_of_parents == 0){
				truth_table.push_back(init_val);	//Make a truth table size 1 with it's initial value
				/** does this work ? **/
			}//end if

			else{
				tt_rows = 1;

				for(int i = 0; i < number_of_parents; i++){
					tt_rows = radixes[i]* tt_rows;
				}//end for

				//Take in truth table
				for(int i = 0; i < tt_rows; i++){
					in >> value;
					//cout << "truth table:  " << value << endl;
					truth_table.push_back(value);
				}//end for

				prob_truth_table.resize(tt_rows);

				if( getVersion() == "GENERALIZED_LOGICAL_NETWORK_VER2" ) {

					for(int i = 0; i < tt_rows; i++) {

						prob_truth_table[i].resize(base);

						for(int j=0; j < base; j ++) {
							in >> prob_truth_table[i][j];
						}
					}

				} else if( getVersion() == "GENERALIZED_LOGICAL_NETWORK_VER1" ) {

					for(int i = 0; i < tt_rows; i++) {

						prob_truth_table[i].resize(base);

						for(int j=0; j < base; j ++) {
							prob_truth_table[i][j] = 0;
						}

						prob_truth_table[i][ truth_table[i] ] = 1.0;

					}
				}
			}//end else

			GeneralizedTruthTable gtt(count, base, parent_vector, radixes, 
                                      parent_time, truth_table, prob_truth_table); 

			m_gtts[count] = gtt;
            nodeParentsIDs[count] = parent_vector;

			initializeNode(count, init_val, base, node_id, node_name, node_type);

			//clear vectors
			truth_table.clear();
			parent_vector.clear();
			radixes.clear();
			parent_time.clear();
			prob_truth_table.clear(); //added by Yang Zhang 2011.1.17 if not clear, then even when truth table is empty, the probability table can still be not empty
			++count;
	
			//check for end of node, if not their print out error message
			char check;
			in >> check;
			if (check != '*'){
				GLNExit(EXIT_FAILURE, "ERROR: GLN format invalid, missing * after node!"); 
			}
		}//end while

	}//end else

	in.close();

    vector<string> parentNames;
    for( size_t node=0; node<nodeParentsIDs.size(); ++node )
    {
        for( const auto & parentID: nodeParentsIDs[node])
            parentNames.push_back( nodeNameList[parentID-1] );

        m_gtts[node].setParentNames( parentNames );
        parentNames.clear();
    }



	//require internal nodes be put before external nodes
	//added by Yang Zhang 2011.10.6
	bool externalFlag = false; // check whether external nodes have appeared
	for(size_t index = 0; index < m_nodes.size(); index++)
	{
		if(m_nodes[index].getType()=='i')
		{
			if(externalFlag)
			{
				cerr << "ERROR: internal nodes must be put before external nodes!" << endl;
				exit(EXIT_FAILURE);
			}
		}else if(m_nodes[index].getType()=='e')
		{
			externalFlag = true;
		}
	}

}//end scanFile

/*
void GeneralizedLogicalNetwork::saveNetworkDotFormat(const string & filename){

	if(m_ready != true || size() <= 0){
		cout << "The network hasn't been initialized with nodes yet, cannot create the file for the network in dot format" << endl;
	}//end if

	else{
		outputNetworkDotFormat(m_nodes, filename);
	}//end else

}//end outputNetworkDotFormat
*/

void GeneralizedLogicalNetwork::saveDot(const string & filename)
// Output the GLN network into dot format
{
	if(m_ready != true || size() <= 0){

		cout << "The network hasn't been initialized with nodes yet, cannot create the file for the network in dot format" << endl;

		return;

	}//end if

	int color=0;
	ofstream out(filename.c_str()); 	//The file the output is going to

	//An array of strings of random colors????
	//string colors[6] = {"purple","blue","green","yellow","orange","red"};

	//Option to change the options below??
	out << "digraph GRN {" << endl
		<< "ratio=auto;" << endl
		<< "margin=0;" <<endl
		<< "node[fontname=\"Helvetica\"];" << endl
		<< "node[fontcolor=\"black\"];" << endl;

	vector<bool> isSingleton(size(), true);

	// Go through the nodes and put out the connection to it's parents
	for(unsigned int i=0; i < size(); i++){

		int dest = i+1;
		string nodeNameDest = m_nodes[i].getName();
		
		// Joe Song, August 2, 2007
		vector<int> delays = m_gtts[i].getDelays();
		
		int nParents = (int) m_gtts[i].getParents().size();

		if(nParents > 0) {
		  
			isSingleton[i] = false;
		  
		}

		for(int j=0; j < nParents; j++){

			int src = m_gtts[i].getParents()[j];

			string nodeNameSrc = m_nodes[ src-1 ].getName();

			isSingleton[src-1] = false;

			if(nodeNameSrc.empty()) {
			}
			
			if(nodeNameSrc.empty() && nodeNameDest.empty()) {

				out	<< "\"" << src << "\" " << "-> " << "\"" << dest << "\"";

			} else if(nodeNameSrc.empty() && ! nodeNameDest.empty()) {

				out	<< "\"" << src << "\" " << "-> " << "\"" << nodeNameDest << "\"";

			} else if(! nodeNameSrc.empty() && nodeNameDest.empty()) {

				out	<< "\"" << nodeNameSrc << "\" " << "-> " << "\"" << dest << "\"";

			} else {
				out	<< "\"" << nodeNameSrc << "\" " << "-> " << "\"" << nodeNameDest << "\"";
			}

			// Joe Song, August 2, 2007
			out << " [label=" << delays[j] << ","   
				<< "arrowhead=normal,color=black,style=solid, ];" << endl;

			// out	<< " [arrowhead=normal,color=black,style=solid, ];" << endl;

		}//end for

		/*
		if(nodeNameDest.empty()) {

			out << "\"" << dest << "\" ";

		} else {

			out << "\"" << nodeNameDest << "\" ";

		}

		switch(m_nodes[i].getType()) {
			case 'e':
				out << "[shape=invtriangle,style=filled,fillcolor=tan1];";
				break;
			case 'i':
			default:
				out << "[shape=ellipse,style=filled,fillcolor=orangered3];";
		}
		*/
	    out << endl << endl;

		color++;
		if(color >= 6)
		{ color=0; }		//Reset color

	}//end for

	for(unsigned int i=0; i < m_nodes.size(); i++){

		if(! isSingleton[i]) {
			int dest = i+1;
			string nodeNameDest = m_nodes[i].getName();

			if(nodeNameDest.empty()) {

				out << "\"" << dest << "\" ";

			} else {

				out << "\"" << nodeNameDest << "\" ";

			}

			switch(m_nodes[i].getType()) {
			case 'e':
				out << "[shape=invtriangle,style=filled,fillcolor=tan1];";
				break;
			case 'i':
			default:
				out << "[shape=ellipse,style=filled,fillcolor=wheat];";
			}
			out << endl;

		}
	}	

	out << endl << "}";
	out.close();
}

void GeneralizedLogicalNetwork::printGTTs() const
{
    cout << "Significant generalized truth tables:" << endl;
    
	for(size_t i=0; i < size(); i++){

        cout << i+1 << '(' << getNodeName(i) << ')' << ": ";

		if(m_gtts[i].size() == 0){

			cout << "None.";

		} else {

			for(int j=0; j < (int) m_gtts[i].getTruthValues().size(); j++){
				cout << m_gtts[i].getTruthValues()[j] << ", ";
			}//end for
		}

		cout<< endl;

	}//end for
}

vector<int> ScanExternalNodes(const string & externalNodesFile)
// Read in the external nodes of this network
{
	vector<int> extNodeIndices(0);
    
	ifstream extfile(externalNodesFile.c_str());		//Read in the file with external nodes
	if(!extfile){
		cerr << "ERROR: file \"" << externalNodesFile << "\" does not exist!" << endl;
		GLNExit(EXIT_FAILURE);
	}//end if
    
	//Read the files into a node
	//--Files in the external node folder must equal to or greater than 1, a 0 node is not accepted
	cout << "External nodes in the network: ";
    
	while( !extfile.eof() ){
		int node_index;
		extfile >> node_index;
		if(node_index == 0){
			cerr << "You cannot enter a 0 node in the external nodes file. The node must be greater than 0" << endl;
			GLNExit(EXIT_FAILURE);
		}//end if
		extNodeIndices.push_back(node_index);
		cout << node_index << ", ";
        
	}//end while
	cout << endl;
	//print the external nodes?
    
	return extNodeIndices;
}

vector<string> ScanListFile(const char * listfile) 
{
	vector<string> vec;
    
	ifstream in(listfile);
    
	if(!in){
		cerr << "ERROR: Cannot open list file \"" << listfile << "\"" << endl;
        GLNExit(EXIT_FAILURE);
	}
    
	while( !in.eof() ){
		string str;
		in >> str; 
		vec.push_back(str);
	}
    
	in.close();
    
	return vec;
}//end scanListFile

void outputStatesDotFormat(const list< vector<int> > & lst, const string & filename){
    
	int i=0;
	// filename = filename + ".DOT";	//User must not send an extention or must handle it incase they do
	ofstream out(filename.c_str()); 	//The file the output is going to
    list< vector<int> >::const_iterator p;
    
	//An array of strings of random colors????
	//string colors[6] = {"purple","blue","green","yellow","orange","red"};
    
	//Option to change the options below??
	out << "digraph GRN {" << endl
    << "ratio=auto;" << endl
    << "margin=0;" <<endl
    << "node[fontname=\"Helvetica\"];" << endl;
    
	p = lst.begin();
	while(p != lst.end()){
		out	<< "\"" << (*p)[0] << "\" " << "-> " << "\"" << (*p)[1] << "\""
        << " [arrowhead=normal,color=black,style=solid, ];" << endl
        << "\"" << (*p)[0] << "\" " 
        << "[shape=box,style=filled,fillcolor=lightsalmon1];" << endl << endl;
        
		p++;			//increment iterator
		i++;			//increment for colors
		if(i >= 6)
		{ i=0;	}		//Reset i
	}//end while
    
    
	out << "\n}";
	out.close();
}

void GeneralizedLogicalNetwork::saveEdgesDotFormat(const string & filename){
    
	if(m_ready != true || size() <= 0){
		cout << "The network hasn't been initialized with nodes yet, cannot create the file for edges dot format" << endl;
	}//end if
    
	else{
		list< vector<int> > edges;
        
		edges = computeStateTrans();
		outputStatesDotFormat(edges, filename);
	}//end else
    
}//end 


