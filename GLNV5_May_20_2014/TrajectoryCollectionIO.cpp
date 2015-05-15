// TrajectoryCollectionIO.cpp
//
// Joe Song
// Created: September 27, 2008
// Last updated: October 5, 2011. Replaced save() by using C++ iostream, instead of C file I/O.

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
using std::ios;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <cassert>
#include <cstdlib>

#include "TrajectoryCollection.h"
#include "GLNGlobals.h"

void readTrajForCompare(string filename, int& offset, vector< vector<int> >& traj_table, 
						ofstream& out, int & num_traj, int &numb_of_init, vector<int> &numbTraj, 
						vector<unsigned short> &maxsteps, vector< vector<unsigned short> > &externalNodes, 
						int &intNodes, int &extNodes)
{
	ifstream in(filename.c_str());
	if(!in){
		cerr << "Cannot open the initial trajectory file, make sure the file is in the same directory as the program \n";
		GLNExit(EXIT_FAILURE);
	}//end if
    
	string ver("");
	int numTraj;  //number of trajectories
	int intNode;  //internal node count
	int extNode;  //external node count
	int base;    //node base
	int numNodes;  //total node count
	int rows, columns = 0;
	int value;  //number of initial trajectories
	int count = 0;  //offset counter variable
	string nodeName("");
	vector<int> temp;
	vector<unsigned short> tempExt;
    
	in >> ver;
	if(ver != "TRAJECTORY_VER2" && ver != "TRAJECTORY_VER1"){
		cerr << "ERROR: Neither TRAJECTORY_VER1 nor TRAJECTORY_VER2 trajecotry file!" << endl;
		GLNExit(EXIT_FAILURE);
	}
	//header
	out << ver << endl;
    
	//number of trajectories
	in >> numTraj;
	out << numTraj << " ";
    
	num_traj = numTraj;
    
	//internal node count
	in >> intNode;
	out << intNode << " ";
	intNodes = intNode;
	//external node count
	if(ver == "TRAJECTORY_VER2"){
		in >> extNode;
		out << extNode << endl;
	}
	else
		extNode = 0;
	extNodes = extNode;
	numNodes = intNode + extNode;
	//print out bases to trajectory file
	for (int i = 0; i < numNodes; i++){
		in >> base;
		out << base << " ";
	}//end for
    
	out << endl;
	if(ver == "TRAJECTORY_VER2"){
		//read in internal node names, print to trajectory file
		for (int i = 0; i < intNode; i++){
			in >> nodeName;
			out << nodeName << " ";
		}//end for
        
		out << endl;
		//read in external node names, print to trajectory file
		for (int i = 0; i < extNode; i++){
			in >> nodeName;
			out << nodeName << " ";
		}//end for
        
		out << endl;
	}
	//int count = 0;
	for (int i = offset; i < 0; i++)
		count++;
    
	offset = count;  //set offset to positive offset
    
    
	//int rows, columns = 0;
	//int value;
    
	columns = numNodes;
	//vector<unsigned short> temp;
    
	//read in initial trajectories and store them in a trajectory table
	for (int k = 0; k < numTraj; k++){
		in >> value;
        
		maxsteps.push_back(value);
		// cout << "max steps value is:  " << maxsteps[k] << endl;
		numb_of_init = value;
		rows = offset;
		numbTraj.push_back(offset);
		for(int i=0; i < rows && i < numb_of_init; i++){ // Joe Song, 12/8/2007: added "&& i < numb_of_init"
			for(int j=0; j < columns; j++){
				//if(j <= intNode){
				in >> value;
				temp.push_back(value);
				//}
				//else{
				// in >> value;
				// tempExt.push_back(value);
				//}
			}//end for
			traj_table.push_back(temp);
			temp.clear();
			//externalNodes.push_back(tempExt);
			//tempExt.clear();
		}//end for
		// cout << "external nodes:  " ;
		for (int i = 0; i < numb_of_init-offset; i++){
			for(int j = 0; j < columns; j++){
				if(j < intNode){
					in >> value;
				}
				else{
					in >> value;
					tempExt.push_back(value);
					// cout << value <<"\t";
				}
			}
			externalNodes.push_back(tempExt);
			tempExt.clear();
			// cout << "\n";
		}
	}//end for
	in.close();
	cout << "finished reading trajectory file"  << endl;
}

vector< vector<int> > TrajectoryCollection::readTraj(ifstream &in, int iteration,int &columns, 
											   		int &rows, int &numb_of_traj)
{
	int value;
	vector< vector<int> > traj_table;

	// int rows_total=0;
	int intcount;
	int extcount;
	string nodeName;
	string keyword;
	string temp;

	/** Possible modifications/additions:
	*Check to see if columns match in consecutive files
	*
	*
	***/

	/***************************************************************************
	****************************************************************************
	Check to see if the string is the valid keyword
	*Done on **first** traj of file IF there are multiple, also reads in:
	*# of trajectories in file
	*# of nodes (samea as columns)
	*All of the bases
	****************************************************************************
	****************************************************************************/
	if(iteration == 0){
		in >> keyword;			//First word should say TRAJECTORY_VER1
		setVersion(keyword);

		setVersion(keyword);
		if((keyword != "TRAJECTORY_VER1") && (keyword != "TRAJECTORY_VER2")){
			cerr << "ERROR: Invalid trajectory file format!" << endl;
			GLNExit(EXIT_FAILURE);
		}//end if
		if(keyword == "TRAJECTORY_VER1"){
			in >> numb_of_traj;	//Read in the number of trajectories
			in >> columns;		//Columns are also the number of nodes
		}//end if
		//case for trajectory_ver2 different inputs
		else{
			in >> numb_of_traj;
			in >> intcount;  //internal node count
			in >> extcount;  //external node count
			columns = intcount + extcount;  //number of columns in trajectory
			for(int i = 1; i <= extcount; i++)
				m_extNodes.push_back(intcount + i);
		}//end else
		//Only want 1 base table
		
		if(m_bases.size() == 0){
			/************************************************
			*****Create Base Table and read in the bases*****
			*************************************************/		
			m_bases.resize(columns);

			for(int j=0; j < columns; j++){
				in >> m_bases[j];
			}//end for
		}//end if
		//read in node names internal then external nodes
		if(keyword == "TRAJECTORY_VER2"){
			for (int i = 0; i < intcount; i++){
				in >> nodeName;
				m_intNodeNames.push_back(nodeName);
			}//end for
			for (int i = 0; i < extcount; i++){
				in >> nodeName;
				m_extNodeNames.push_back(nodeName);
			}//end for
		}//end if


	}//end if


	/***********************************************************************
	***********Can read in multiple trajectories in the same file***********
	************************************************************************/

	/***Read in the rows and columns***/
	in >> rows;				//read in the rows
	
	// traj_table = createTable(rows, columns);	//Create traj table
	traj_table.resize(rows);
	for(int i=0; i<rows; i++) { traj_table[i].resize(columns); }

	//reads the data that goes into traj_table
	for(int i=0; i < rows; i++){

		//Goes through the columns
		for(int j=0; j < columns; j++){
			in >> value;
			traj_table[i][j]=value;
			//cout << traj_table[i][j] << " ";
		}//end for
		//cout << "\n";
	}//end for

	return traj_table;
}//end readfile


vector< vector<int> > TrajectoryCollection::readTraj(ifstream &in, int iteration,int &columns,int &rows, int &numb_of_traj, vector <bool> &beParent,vector <bool> &haveParent)
//added by YangZhang 4/24/2009
//add two attributes saying whether or not one node can have parent and be parent
{
	int value;
	vector< vector<int> > traj_table;
	// int rows_total=0;
	int intcount;
	int extcount;
	string nodeName;
	string keyword;
	string temp;

	string be_parent("");
	string have_parent("");
	string excludeSign("");
	vector<int> excludenodes;

	if(iteration == 0){
		in >> keyword;			//First word should say TRAJECTORY_VER
		setVersion(keyword);
		if((keyword != "TRAJECTORY_VER1") && (keyword != "TRAJECTORY_VER2") 
			&& (keyword != "TRAJECTORY_VER2.1") && (keyword != "TRAJECTORY_VER2.2")){
			cerr << "ERROR: Invalid trajectory file format!" << endl;
			GLNExit(EXIT_FAILURE);
		}//end if
		if(keyword == "TRAJECTORY_VER1"){
			in >> numb_of_traj;	//Read in the number of trajectories
			in >> columns;		//Columns are also the number of nodes
		}//end if
		//case for trajectory_ver2 different inputs
		else{
			in >> numb_of_traj;
			in >> intcount;  //internal node count
			in >> extcount;  //external node count
			columns = intcount + extcount;  //number of columns in trajectory
			for(int i = 1; i <= extcount; i++)
				m_extNodes.push_back(intcount + i);
		}//end else
		//Only want 1 base table
		
		if(m_bases.size() == 0){
			/************************************************
			*****Create Base Table and read in the bases*****
			*************************************************/		
			m_bases.resize(columns);

			for(int j=0; j < columns; j++){
				in >> m_bases[j];
			}//end for
		}//end if
		//read in node names internal then external nodes
		if(keyword == "TRAJECTORY_VER2"|| keyword == "TRAJECTORY_VER2.1" || keyword == "TRAJECTORY_VER2.2"){
			for (int i = 0; i < intcount; i++){
				in >> nodeName;
				m_intNodeNames.push_back(nodeName);
			}//end for
			for (int i = 0; i < extcount; i++){
				in >> nodeName;
				m_extNodeNames.push_back(nodeName);
			}//end for
		}//end if
		
		if(keyword == "TRAJECTORY_VER2.1"|| keyword == "TRAJECTORY_VER2.2"){
			for (int i = 0; i < (intcount + extcount); i++){
				in >> be_parent;
				beParent.push_back((be_parent=="Yes"||be_parent=="YES")?true:false);
			}
			for (int i = 0; i < (intcount + extcount); i++){
				in >> have_parent;
				haveParent.push_back((have_parent=="Yes"||have_parent=="YES")?true:false);	
			}	
		}
        /*
        //added by Haizhou Wang, June 27, 2012
        //To handle the reconstruction of TRAJECTORY_VER1
        //But there is a "vector out of range" error

        else if( keyword == "TRAJECTORY_VER1") 
        { 
            beParent.assign(columns, true);
            haveParent.assign(columns, true);
        }
        */
		else{
			beParent.assign(intcount+extcount, true);
			haveParent.assign(intcount+extcount, false);
			for(int inindex=0; inindex<intcount; inindex++)
			{
				haveParent[inindex] = true;
			}
		}
	}//end if


	/***********************************************************************
	***********Can read in multiple trajectories in the same file***********
	************************************************************************/

	/*read in which nodes are excluded in this trajectory*/
	if(getVersion() == "TRAJECTORY_VER2.2")
	{
		in >> value;
		in >> excludeSign;
		if(excludeSign.compare("exclude")&&excludeSign.compare("EXCLUDE")
			&&excludeSign.compare("excludes")&&excludeSign.compare("EXCLUDES"))
		{
			cerr << "ERROR: Invalid trajectory file format!" << endl;
			GLNExit(EXIT_FAILURE);
		}
		if(value>0)
		{
			excludenodes.resize(value);
			for(int i=0; i<value;i++)
			{
				in>>excludenodes[i];
				excludenodes[i]--;//node id starts from 1,but index should starts from 0
			}
			m_excludeTrajectory.push_back(excludenodes);
		}
		else if(value==0)
		{
			excludenodes.resize(1);
			excludenodes[0] = -1;
			m_excludeTrajectory.push_back(excludenodes);
		}
	}
	/***Read in the rows and columns***/
	in >> rows;				//read in the rows
	
	// traj_table = createTable(rows, columns);	//Create traj table
	traj_table.resize(rows);
	for(int i=0; i<rows; i++) { traj_table[i].resize(columns); }

	//reads the data that goes into traj_table
	for(int i=0; i < rows; i++){

		//Goes through the columns
		for(int j=0; j < columns; j++){
			in >> value;
			traj_table[i][j]=value;
			//cout << traj_table[i][j] << " ";
		}//end for
		//cout << "\n";
	}//end for

	return traj_table;
}//end readfile




void readInitTraj(string filename, int& offset,  vector< vector<unsigned short> >&  traj_table, 
				  ofstream& out, int& num_traj,int &numb_of_init, vector<int> &numbTraj, 
				  vector< vector<int> > &externalNodes, int &intNodes,
				  int &extNodes, vector<int> & max_steps) 
{
	ifstream in(filename.c_str());
	if(!in){
		cerr << "Cannot open the initial trajectory file, make sure the file is in the same directory as the program \n";
		GLNExit(EXIT_FAILURE);
	}//end if

	string ver("");
	int numTraj;  //number of trajectories
	int intNode;  //internal node count
	int extNode;  //external node count
	int base;    //node base
	int numNodes;  //total node count
	int rows, columns = 0;
	int value;  //number of initial trajectories
	int count = 0;  //offset counter variable
	string nodeName("");
	vector<unsigned short> temp;
	vector<int> tempExt;
	in >> ver;
	if(ver != "TRAJECTORY_VER2"){
		cerr << "This is not a 'TRAJECTORY_VER2' file" << endl;
		GLNExit(EXIT_FAILURE);
	}
	//header
	out << ver << endl;

	//number of trajectories
	in >> numTraj;
	out << numTraj << " ";

	num_traj = numTraj;

	//internal node count
	in >> intNode;
	out << intNode << " ";
	intNodes = intNode;

	//external node count
	in >> extNode;
	out << extNode << endl;
	extNodes = extNode;

	numNodes = intNode + extNode;
	//print out bases to trajectory file
	for (int i = 0; i < numNodes; i++){
		in >> base;
		out << base << " ";
	}//end for

	out << endl;
	//read in internal node names, print to trajectory file
	for (int i = 0; i < intNode; i++){
		in >> nodeName;
		out << nodeName << " ";
	}//end for

	out << endl;
	//read in external node names, print to trajectory file
	for (int i = 0; i < extNode; i++){
		in >> nodeName;
		out << nodeName << " ";
	}//end for

	out << endl;
	//int count = 0;
	for (int i = offset; i < 0; i++)
		count++;

	offset = count;  //set offset to positive offset


	//int rows, columns = 0;
	//int value;

	columns = numNodes;
	//vector<unsigned short> temp;
	rows = offset;
	//read in initial trajectories and store them in a trajectory table
	for (int k = 0; k < numTraj; k++){
		in >> value;
		max_steps.push_back(value);
		numb_of_init = value;
		//rows = value;
		//numbTraj.push_back(value);
		numbTraj.push_back(offset);
		for(int i=0; i < rows; i++){
			for(int j=0; j < columns; j++){
				in >> value;
				temp.push_back(value);
			}//end for
			traj_table.push_back(temp);
			temp.clear();
		}//end for
		//}//end for
		//cout << "external node values:  " << endl;
		for (int i = 0; i < numb_of_init-offset; i++){
			for(int j = 0; j < columns; j++){
				if(j < intNode){
					in >> value;
				}
				else{
					in >> value;
					tempExt.push_back(value);
					// cout << value <<"\t";
				}
			}
			externalNodes.push_back(tempExt);
			tempExt.clear();
			// cout << "\n";
		}
	}//end for
	in.close();
}//end readInitTraj

size_t TrajectoryCollection::scan(const vector<string> & files)
{
	size_t i;

	for(i=0; i<files.size(); i++) {
		scan( files[i] );
	}

	return i; // total number of files read
}

bool TrajectoryCollection::scan(const string & file)
{
	m_files.push_back(file);

	ifstream in(file.c_str());

	if(!in.is_open()) {
        cerr << "ERROR: openning file \"" << file << "\"" << endl;
		GLNExit(EXIT_FAILURE);
	}

	//Checks for multiple trajectories in each file -- numb_of_traj will update upon funct call
	int num_of_traj=1;	//

	for(int j=0; j < num_of_traj; j++){
		//Read in the trajectory
		int lenTraj = 0;		// the length of the trajectory to be read
		vector< vector<int> > p_traj_table;  // pointer to the trajectory to be read

		//p_traj_table = readTraj(in, j, m_nNodes, lenTraj, num_of_traj); 
		//modified by YangZhang 4/24/2009
		p_traj_table = readTraj(in, j, m_nNodes, lenTraj, num_of_traj, m_beParent, m_haveParent); 

		if(p_traj_table.size()>0){
			m_trajTables.push_back(p_traj_table);	//Put the table into the node
			// m_trajLengths.push_back(lenTraj);	Removed MS 4/28/2013 //The row from readTraj in prev line
		}else{
			cerr << "ERROR: Invalid trajectory collection file \"" << file << "\"!" << endl;
            GLNExit(EXIT_FAILURE);
		} 

	}//end for

	in.close();

    // Check if the trajectories has any missing values
    checkMissingValues();
    
	return true; // one file read
}

void TrajectoryCollection::save(const char * file) const
// add by yangzhang 11.12.2008,save trajectory collection to .trj file
// Joe Song Oct 5, 2011: changed from C file IO to C++ file IO
{
	ofstream ofs;

	if(file==NULL) {
		ofs.open("outputTrj.txt", ios::out);
	} else {
		ofs.open(file, ios::out);
	}

	if(! ofs.is_open() ) {
		cerr << "ERROR: cannot open file \"" << file << "\" for writing!" << endl;
		GLNExit(EXIT_FAILURE);
	}

	if(m_version.find("2") != string::npos) { //version2 or version 2.1

		ofs << m_version << endl;

		ofs << m_trajTables.size() << ' ' << m_intNodeNames.size() << ' ' 
			<< m_extNodeNames.size() << endl;

		for(size_t i=0; i < m_bases.size(); i++) {
			ofs << m_bases[i] << '\t';
		}

		ofs << endl;

		for(size_t i=0; i < m_intNodeNames.size(); i++) {
			ofs << m_intNodeNames[i].c_str() << '\t';
		}

		ofs << endl;

		for(size_t i=0; i < m_extNodeNames.size(); i++) {
			ofs << m_extNodeNames[i].c_str() << '\t';
		}

		if(m_extNodeNames.size()>0)	{
			ofs << endl;
		}

        if( m_version.find("1") != string::npos ) //Take care of version 2.1
        {
            ofs<<endl;
            for( const auto & haveP: m_haveParent )
            {
                if( haveP == true )
                    ofs << "Yes\t";
                else
                    ofs << "No\t";
            }
            ofs<<endl;
            for( const auto & beP: m_beParent )
            {
                if( beP == true )
                    ofs << "Yes\t";
                else
                    ofs << "No\t";
            }
            ofs<<endl<<endl;
        }

		for(size_t i=0; i<m_trajTables.size(); i++)	{

			ofs << m_trajTables[i].size() << endl;

			for(size_t j=0; j < m_trajTables[i].size(); j++) {
				for(size_t k=0; k < m_trajTables[i][j].size(); k++) {
					ofs << m_trajTables[i][j][k] << '\t';
				}
				ofs << endl;
			}
		}

	}
}

void TrajectoryCollection::print() const
// printing the content of a trajectory collection
{
	cout << "The number of trajectories is: " << m_trajTables.size() << endl;

	// Print out all the trajectory tables and then bases
	for(size_t i=0; i < m_trajTables.size(); i++){
		// for(int j=0; j < m_trajLengths[i]; j++) {  MS 4/28/2013
        for(size_t j=0; j < m_trajTables[i].size(); j++) {
			for(int k=0; k < nNodes(); k++) { 
				cout << m_trajTables[i][j][k] << " ";
			}
			cout << endl;
		}//end for
		cout << endl;
	}//end for

	cout << "bases: " << endl;

	// Print out the base table
	for(int i=0; i < nNodes(); i++){
		cout << base(i) << " ";
	}//end for
	cout << endl;

	cout << "_____________________external node info_______________" << endl;
	if(nNodes()){
		cout << "There are " << m_extNodes.size() << " external nodes:" << endl;
		for (size_t i = 0; i < m_extNodes.size(); i++)
			cout << m_extNodes[i] << ", ";
        cout << endl;
	}
}
