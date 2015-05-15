// TrajectoryCollection.cpp
//
// Joe Song
// Created: September 27, 2008
// Last modified: MS November 27, 2011

#include <iostream>
#include <ctime>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <map>
#include <set>
#include <utility>

using std::multimap;
using std::map;
using std::pair;
using std::set;

#include "TrajectoryCollection.h"

const vector<string> TrajectoryCollection::getNodeNames(const vector<int> & ids) const
{
    vector<string> names(ids.size());
    
    for(size_t j=0; j<ids.size(); j++) {
        
        if (ids[j] <= (int) m_intNodeNames.size()) {
            names[j] = m_intNodeNames[ids[j]-1];
        } else { // Assuming the external nodes always come after the internal nodes
            names[j] = m_extNodeNames[ids[j]- 1 - m_intNodeNames.size()];
        }
    }
    return names;    
}

void TrajectoryCollection::setTrajTables(const vector< vector< vector<int> > > & trjTables) 
{ // Revised by MS Nov 27, 2011. Added by yangzhang 11.15.2008

    m_trajTables = trjTables;
    
    /* Removed MS 4/28/2013
    if( m_trajLengths.size() != trjTables.size() ) {
        m_trajLengths.resize(trjTables.size());
    }

    for (size_t i=0; i<trjTables.size(); i++) {
        m_trajLengths[i] = trjTables[i].size();
    }
    */
    
    checkMissingValues();
}

bool TrajectoryCollection::checkMissingValues()
{
    m_hasMissingValues = false;
    
    int value = m_missingValue;
    
    for(size_t j=0; j<m_trajTables.size(); j++) {
        
        for(size_t t=0; t<m_trajTables[j].size(); t++) {
            
            if (m_trajTables[j][t].end() != find(m_trajTables[j][t].begin(),
                                                 m_trajTables[j][t].end(),
                                                 value) ) {
                m_hasMissingValues = true;
                break;
            }
        }
    }
    
    return m_hasMissingValues;
}

int compareTraj(const string & input1, const string & input2)
//modified to apply for version 2 by yangzhang 11.19.2008
{
	string temp;
	int rows, rows2;
	int intnodes1, intnodes2,extnodes1,extnodes2;
	int numboftrj1, numboftrj2;
	int value, value2;
	int different=0;			//Inc everytime it's different
	// bool status=true;

	ifstream file1(input1.c_str());
	ifstream file2(input2.c_str());

	//Eat up header
	file1 >> temp;
	file2 >> temp;

	//Eat up # of trajectories
	file1 >> numboftrj1;
	file2 >> numboftrj2;

	if(numboftrj1 != numboftrj2){
		cerr << "The number of nodes are not the same, therefore the trajectories can't be the same";
		return false;
	}//end if

	//Take in the number of internal nodes
	file1 >> intnodes1;
	file2 >> intnodes2;

	//take in number of external nodes
	file1 >> extnodes1;
	file2 >> extnodes2;

	if((intnodes1+extnodes1) != (intnodes2+extnodes2)){
		cerr << "The number of nodes are not the same, therefore the trajectories can't be the same";
		return false;
	}//end if

	//Eat up the bases
	for(int i=0; i < intnodes1+extnodes1; i++){
		file1 >> temp;
		file2 >> temp;
	}//end for

	//Eat up the node names
	for(int i=0; i < intnodes1+extnodes1; i++){
		file1 >> temp;
		file2 >> temp;
	}//end for

	for(int i=0;i<numboftrj1;i++)
	{
		//Take in the rows
		file1 >> rows;
		file2 >> rows2;

		if(rows != rows2){
			cerr << "The number of rows are not the same, therefore the trajectories can't be the same";
			return false;
		}//end if

		for(int i=0; i < rows; i++){

			for(int j=0; j < intnodes1+extnodes1; j++){
				file1 >> value;
				file2 >> value2;

				if(value != value2){
					different++;
					// cout << "Different at << row: " << i+1 << ">> and <<column: " << j+1 << ">>" << endl;
					// status=false;
				}//end if		
			}//end for

			//cout << endl;
		}//end for
	}
	cout << "Total different states: " << different << endl;
	return different;
}//end compareTraj

vector< vector<int> > compareTrajectoriesValues(string traj1, string traj2, string difffile, int offset)
{
	vector<vector< int> > trajectorydifferences;
	cout << "in compare traj values"  << endl;
	ifstream in1(traj1.c_str());
	if(!in1){
		cerr << "Cannot open the trajectory file, make sure the file is in the same directory as the program \n";
		exit(EXIT_FAILURE);	
	}//end if

	ifstream in2(traj2.c_str());
	if(!in2){
		cerr << "Cannot open the trajectory file, make sure the file is in the same directory as the program \n";
		exit(EXIT_FAILURE);
	}//end if

	ofstream out (difffile.c_str());

	string ver("");
	// int offset;
	int numTraj;  //number of trajectories
	int intNode;  //internal node count
	int extNode;  //external node count
	int base;    //node base
	int numNodes;  //total node count
	int columns = 0;
	int value;
	int value1;  
	int value2;
	int count = 0;  //offset counter variable
	string nodeName("");
	//vector<unsigned short> temp;
	//vector<unsigned short> tempExt;

	in1 >> ver;
	in2 >> ver;
	if(ver != "TRAJECTORY_VER2" && ver != "TRAJECTORY_VER1"){
		cerr << "ERROR: Invalid 'TRAJECTORY_VER2' file \"" << traj2 << "\"" << endl;
		exit(EXIT_FAILURE);
	}
	//header
	out << ver << endl;

	//number of trajectories
	in1 >> numTraj;
	in2 >> numTraj;
	out << numTraj << " ";

	//num_traj = numTraj;

	//internal node count
	in1 >> intNode;
	in2 >> intNode;
	out << intNode << " ";
	//intNodes = intNode;
	//external node count
	if(ver == "TRAJECTORY_VER2"){
		cout << "in external nodes:  " << endl;
		in1 >> extNode;
		in2 >> extNode;
		out << extNode << endl;
	}
	else
		extNode = 0;
	//out << extNode << endl;
	//extNodes = extNode;
	numNodes = intNode + extNode;
	cout << "numNodes" << numNodes << endl;
	//print out bases to trajectory file
	for (int i = 0; i < numNodes; i++){
		in1 >> base;
		in2 >> base;
		out << base << " ";
	}//end for

	out << endl;
	if(ver == "TRAJECTORY_VER2"){
		//read in internal node names, print to trajectory file
		for (int i = 0; i < intNode; i++){
			in1 >> nodeName;
			in2 >> nodeName;
			out << nodeName << " ";
		}//end for

		out << endl;
		//read in external node names, print to trajectory file
		for (int i = 0; i < extNode; i++){
			in1 >> nodeName;
			in2 >> nodeName;
			out << nodeName << " ";
		}//end for
	}
	out << endl;
	//int count = 0;
	for (int i = offset; i < 0; i++)
		count++;

	offset = count;  //set offset to positive offset


	//int rows, columns = 0;
	//int value;

	columns = numNodes;
	cout << "columns are:  " << columns << endl;
	//vector<unsigned short> temp;
	out << endl;
	//read in initial trajectories and store them in a trajectory table
	for (int k = 0; k < numTraj; k++){
		in1 >> value;
		in2 >> value;
		out << value << endl;
		//cout << "value1 value2:  " << value1 <<";" << value2 << endl;
		//if(value1 == value2){
		vector<int> temp;
		for(int i=0; i < value; i++){
			for(int j=0; j < columns; j++){
				in1 >> value1;
				in2 >> value2;
				//if(value1 >= value2){
				out << abs(value1 - value2) << "\t";
				//cout <<  abs(value1 - value2) << "\t";
				//}
				// cout << value1 - value2 << " ";
				temp.push_back(value1 - value2);
				//else{
				//out << value2 - value1 << "\t";
				//cout << value2 - value1 << "\t";
				//}
			}
			out << endl;
			// cout << endl;
			trajectorydifferences.push_back(temp);
			temp.clear();
		}
		out << "\n\n";
		// cout << "\n\n";
	}//end for

	in1.close();
	in2.close();
	out.close();
	// cout << "created diff table"  << endl;
	for (size_t i = 0; i < trajectorydifferences.size(); i++){
		for (int j = 0; j < numNodes; j++){
			// cout << trajectorydifferences[i][j] << " ";

		}
		// cout << endl;
	}

	return trajectorydifferences;
}

void summary(vector < vector< int> > &diff, int numNodes, string summaryfile)
{
	cout << "in summary" << endl;
	ofstream out (summaryfile.c_str());
	cout << "size:  " << diff.size() << "num:  " << numNodes << endl;
	vector<int> min(numNodes);
	vector<int> max(numNodes);
	vector<float> mean(numNodes);
	double overallMean = 0.0;
	double overallVariance = 0.0;
	vector<double> variance(numNodes);
	vector<double> varianceresult(numNodes);
	int overallMin = 20;
	int overallMax = -20;

	for (int i = 0; i < numNodes; i++){
		min [i] = 20;
		max [i] = -20;
		mean [i] = 0;
		variance[i] = 0;
		varianceresult [i] = 0.0;
	}
	cout << "calculations" << endl;
	for (size_t i = 0; i < diff.size(); i++){
		// cout << "i: " << i << endl;
		for (int j = 0; j < numNodes; j++){

			// cout << diff[i][j] << " ";
			if(min[j] > diff[i][j])
				min [j] = diff[i][j];
			if(max[j] < diff[i][j])
				max [j] = diff[i][j];
			mean[j] +=diff[i][j];

		}
		// cout << endl;
	}

	for(int i = 0; i < numNodes; i++)
		mean [i] = mean[i] / numNodes;

	for (size_t i = 0; i < diff.size(); i++){
		for (int j = 0; j < numNodes; j++){
			// cout << diff[i][j] << ":" << mean [j] << "; ";
			variance[j] += diff[i][j] - mean[j]; 	
			//cout << "variance at j" << variance[j] << endl;
		}
		// cout << endl;
	}	

	for (int j = 0; j < numNodes; j++){
		// cout << variance[j] <<endl;
		variance[j] = pow(variance[j],2.0); 	
		// cout << "variance at j" << variance[j] << endl;
	}	

	for (int j = 0; j < numNodes; j++){
		varianceresult[j] = variance[j]/numNodes;	
	}


	for (size_t i = 0; i < diff.size(); i++){
		//cout << "i: " << i << endl;
		for (int j = 0; j < numNodes; j++){
			overallMean += diff[i][j] * i;

		}
		//cout << endl;
	}
	cout << "overall mean is:  " << overallMean << endl;
	//overallMean = 0.0;
	overallMean = (1.0/(numNodes*(diff.size()-1))) * overallMean;
	//overallMean = 1/overallMean;
	cout << "final overal mean:  " << overallMean << endl;

	for (size_t i = 0; i < diff.size(); i++){
		for (int j = 0; j < numNodes; j++){
			//cout << diff[i][j] << ":" << mean [j] << "; ";
			overallVariance += diff[i][j] - overallMean; 	
			//cout << "variance at j" << variance[j] << endl;
		}
		//cout << endl;
	}	

	overallVariance = pow(overallVariance,2.0); 	

	overallVariance = overallVariance * (1.0/(numNodes*(diff.size()-1)));	

	for (int i = 0; i < numNodes; i++)
		if(min[i] < overallMin)
			overallMin = min[i];

	for (int i = 0; i < numNodes; i++)
		if(max[i] > overallMax)
			overallMax = max[i];

	//out << "minimum error for each node" << endl;
	for (int i = 0; i < numNodes; i++)
		out << min[i] << "\t";

	out << overallMin;
	out << endl;

	//out << "maximum error for each node" << endl;
	for (int i = 0; i < numNodes; i++)
		out << max[i] << "\t";

	out << overallMax;
	out << endl;

	//out << "mean error for each node" << endl;
	for (int i = 0; i < numNodes; i++)
		out << mean[i] << "\t";

	out << overallMean;
	out << endl;

	//out << "variance for each node" << endl;
	for (int i = 0; i < numNodes; i++)
		out << varianceresult[i] << "\t";

	out << overallVariance;
	out << endl;

	//out <<"overall minimum"<<endl;
	//out << overallMin << endl;

	//out <<"overall maximum" << endl;
	//out << overallMax << endl;

	//out << "overall mean" << endl;
	//out << overallMean << endl;

	//out << "overall variance" << endl;
	//out << overallVariance << endl;

	cout << "summary finished"  << endl;

	out.close();
}

TrajectoryCollection::TrajectoryCollection (const vector<string> & listTrajFiles)
: m_hasMissingValues(true)
{
	scan(listTrajFiles);
}

void TrajectoryCollection::allocate(int nTrajs, int lenEachTraj, int nVariables, int initValue)
{
	m_trajTables.resize(nTrajs);

	for(int i=0; i<nTrajs; i++) {

		m_trajTables[i].resize(lenEachTraj);

		for(int j=0; j<lenEachTraj; j++)			{

			m_trajTables[i][j].resize( nVariables );

			for(int k=0; k < nVariables; k++)
			{
				m_trajTables[i][j][k] = initValue; //to denote if needs to be simulated
			}
		}
	}
	// m_trajLengths = vector<int>(nTrajs, lenEachTraj); Removed MS 4/28/2013
    m_nNodes = nVariables;
}

/* change into house noise model for pure random noise
//added by yangzhang 2/8/2009 add noise to one node in trjectory table
//nodeid starts from 0
void TrajectoryCollection::addNoiseToNode(int base,double noiselevel,int nodeid)
{
	vector< vector<double> >  prvector;   //probability vector
	vector< vector<double> >  valuevector;   //if prvector is 0.2,0.5,0.1,0.2 then valuevector is 0.2,0.7,0.8,1
	prvector.resize(base);	
	vector<double> tempsum;
	tempsum.resize(base);
	for(int m=0;m<base;m++)
	{
		for(int n=0;n<base;n++)
		{
				tempsum[m] += (n-m>0)?(n-m):(m-n); 
		}
	}
	for(int m=0;m<base;m++)
	{
		prvector[m].resize(base);
		for(int n=0;n<base;n++)
		{
			if(m==n)
			{
				prvector[m][n] = 1-noiselevel;
			}
			else
			{
				if (base == 2)
				{
					prvector[m][n] = noiselevel;
				}
				else if (base > 2)
				{
					prvector[m][n] = (1-((n-m>0)?(n-m):(m-n))/tempsum[m])*noiselevel/(base-2);
				}
			}
		}
	}
	//transform prvector to valuevector
	double sum = 0; //temp variable
	valuevector.resize(prvector.size());
	for(size_t i=0;i<prvector.size();i++)
	{
		valuevector[i].resize(prvector[i].size());
		sum = 0;
		for(size_t j=0;j<prvector[i].size();j++)
		{
			sum = sum + prvector[i][j];
			valuevector[i][j] = sum;
		}
		//because 0.3 is stored as 0.29999999 in C, and last value in valuevector may not be 1.
		valuevector[i][prvector[i].size()-1] = 1;		
	}


	//srand( (unsigned)time(NULL) );
	for(size_t i = 0; i<m_trajTables.size(); i++)
	{
		for(size_t j = 0; j<(m_trajTables[i]).size();j++)
		{
				double randomnumber = static_cast<double> (rand()) / RAND_MAX;
				int k = 0;
				for(;k<(int) valuevector[m_trajTables[i][j][nodeid]].size();k++) 
				{
					if (randomnumber<= valuevector[m_trajTables[i][j][nodeid]][k])
					{
						break;
					}
					else
					{
						continue;
					}
				}
				m_trajTables[i][j][nodeid] = k;
		}
	}
}*/

//modified by Yang Zhang 7.19.2012 house noise model with more random noise
void TrajectoryCollection::addNoiseToNode(int base,double noiselevel,int nodeid, const string & HouseNoiseModel)
{
	vector< vector<double> >  prvector;   //probability vector
	vector< vector<double> >  valuevector;   //if prvector is 0.2,0.5,0.1,0.2 then valuevector is 0.2,0.7,0.8,1
	prvector.resize(base);	
	vector<double> tempsum;
	tempsum.resize(base);
	for(int m=0;m<base;m++)
	{
		for(int n=0;n<base;n++)
		{
			//tempsum[m] += (n-m>0)?(n-m+base):(m-n+base);  //to make the equation simpler
			tempsum[m] += (n-m>0)?(n-m):(m-n); 
		}
	}

	if(HouseNoiseModel == "house") //house noise model, has simpler form and achieve pure random noise when noise level is 1
	{
		for(int m=0;m<base;m++)
		{
			prvector[m].resize(base);
			for(int n=0;n<base;n++)
			{
				//prvector[m][n] = (1-((n-m>0)?(n-m+base):(m-n+base))/tempsum[m])*noiselevel/(base-1);
				prvector[m][n] = (1-((n-m>0)?(n-m):(m-n))/tempsum[m])*noiselevel/(base-1);
						//    prvector[m][n] = 1/base;  //purely random
			}
			prvector[m][m] = prvector[m][m]+1-noiselevel;
	                     
		}
		//apply house noise model to obtain purely random noise
		for(int m=0;m<base;m++)
		{
			for(int n=0;n<base;n++)
			{
				prvector[m][n] = prvector[m][n]*(1-noiselevel)+noiselevel/base;
			}

		}
	}else //original noise model, not random even when noise level is 1.
	{
		for(int m=0;m<base;m++)
		{
			prvector[m].resize(base);
			for(int n=0;n<base;n++)
			{
				if(m==n)
				{
					prvector[m][n] = 1-noiselevel;
				}
				else
				{
					if (base == 2)
					{
						prvector[m][n] = noiselevel;
					}
					else if (base > 2)
					{
						prvector[m][n] = (1-((n-m>0)?(n-m):(m-n))/tempsum[m])*noiselevel/(base-2);
					}
				}
			}
		}
	}


        //transform prvector to valuevector
	double sum = 0; //temp variable
	valuevector.resize(prvector.size());
	for(size_t i=0;i<prvector.size();i++)
	{
		valuevector[i].resize(prvector[i].size());
		sum = 0;
		for(size_t j=0;j<prvector[i].size();j++)
		{
			sum = sum + prvector[i][j];
			valuevector[i][j] = sum;
		}
		//because 0.3 is stored as 0.29999999 in C, and last value in valuevector may not be 1.
		valuevector[i][prvector[i].size()-1] = 1;		
	}


	//srand( (unsigned)time(NULL) );
	for(size_t i = 0; i<m_trajTables.size(); i++)
	{
		for(size_t j = 0; j<(m_trajTables[i]).size();j++)
		{
				double randomnumber = static_cast<double> (rand()) / RAND_MAX;
				int k = 0;
				for(;k<(int) valuevector[m_trajTables[i][j][nodeid]].size();k++) 
				{
					if (randomnumber<= valuevector[m_trajTables[i][j][nodeid]][k])
					{
						break;
					}
					else
					{
						continue;
					}
				}
				m_trajTables[i][j][nodeid] = k;
		}
	}
}

//added by yangzhang 2/8/2009 add noise to whole trjectory table
void TrajectoryCollection::addNoise(double noiselevel, const string & housemodel)
{
	//srand( (unsigned)time(NULL) );

	for(int i = 0;i< m_nNodes;i++)
	{
		addNoiseToNode(m_bases[i],noiselevel,i, housemodel);
	}
}

void TrajectoryCollection::permuteTrjtableSwapingColumn()
{
	for (size_t k = 0; k < m_trajTables[0][0].size(); ++k) {
		int m = rand() % (k + 1);
		//do the actual permute job 
		for(size_t trjIndex = 0; trjIndex < m_trajTables.size(); trjIndex++)
		{
			for(size_t rowIndex = 0; rowIndex < m_trajTables[trjIndex].size(); rowIndex++)
			{
				int temp = m_trajTables[trjIndex][rowIndex][k];
				m_trajTables[trjIndex][rowIndex][k] = m_trajTables[trjIndex][rowIndex][m];
				m_trajTables[trjIndex][rowIndex][m] = temp;
			}
		}
	}
	/* the following code only permute gene names, but won't work for the child working zone
	vector<string> NodeNames = m_trajCols[0].getIntNodeNames();

	for(size_t j=0; j<gc.m_listTrajFiles.size(); j++) {
		for (int k = 0; k < NodeNames.size(); ++k) {
			int m = rand() % (k + 1);
			string tempName = NodeNames[k];
			NodeNames[k] = NodeNames[m];
			NodeNames[m] = tempName;
		}
		m_trajCols[j].setIntNodeNames(NodeNames);
		}
	*/

	//gc.m_pathway.scan(filePathway.c_str(), m_trajCols[0].getIntNodeNames());

	//gc.m_candidateTopology.scan(fileCandidate.c_str(),m_trajCols[0].getIntNodeNames());

}

void TrajectoryCollection::permuteTrjtableOnTime()
{
	//convert <vector<vector<int>> data structure into vector<vector<int>> for easier permutation
	vector<vector<int> > trjdata;
	for(size_t trjIndex = 0; trjIndex < m_trajTables.size(); trjIndex ++)
	{
		for(size_t rowIndex = 0; rowIndex<m_trajTables[trjIndex].size(); rowIndex++)
		{
			trjdata.push_back(m_trajTables[trjIndex][rowIndex]);
		}
	}

	//permute two dimentional trjdata
	if(trjdata.size() == 0)
	{
		cerr << "ERROR: trajectory file is empty!" << endl;
		exit(EXIT_FAILURE);
	}
	for(size_t k=0; k<trjdata[0].size(); k++)
	{
		for(size_t rowIndex = 0; rowIndex < trjdata.size(); rowIndex++)	
		{
				int m = rand() % (rowIndex + 1);
				int temp = trjdata[rowIndex][k];
				trjdata[rowIndex][k] = trjdata[m][k];
				trjdata[m][k] = temp;
			}
	}

	//put permuted trj data back
	int rowIndexInPermutedData = 0;
	for(size_t trjIndex = 0; trjIndex < m_trajTables.size(); trjIndex ++)
	{
		rowIndexInPermutedData = (trjIndex>0)?(rowIndexInPermutedData + m_trajTables[trjIndex-1].size()):rowIndexInPermutedData;
		for(size_t rowIndex = 0; rowIndex < m_trajTables[trjIndex].size(); rowIndex++)
		{
			for(size_t colIndex = 0; colIndex < m_trajTables[trjIndex][rowIndex].size(); colIndex++)
			{
				m_trajTables[trjIndex][rowIndex][colIndex] = trjdata[rowIndexInPermutedData+rowIndex][colIndex];
			}
		}
	}

}

void TrajectoryCollection::permuteTrjtableOnColumn()
{
	for (size_t k = 0; k < getIntNodeNames().size(); ++k) {
		//do the actual permute job 
		for(size_t trjIndex = 0; trjIndex < m_trajTables.size(); trjIndex++)
		{
			for(size_t rowIndex = 0; rowIndex < m_trajTables[trjIndex].size(); rowIndex++)
			{
				int m = rand() % (k + 1);
				int temp = m_trajTables[trjIndex][rowIndex][k];
				m_trajTables[trjIndex][rowIndex][k] = m_trajTables[trjIndex][rowIndex][m];
				m_trajTables[trjIndex][rowIndex][m] = temp;
			}
		}
	}

}

//this function will append the second trj file to the first one by row
vector<vector<vector<int> > > concatenateTrjs(vector<vector<vector<int> > > & firstTrj, vector<vector<vector<int> > > & secondTrj)
{
	vector<vector<vector<int> > > trjTables;

	if(firstTrj[0][0].size() != secondTrj[0][0].size())
	{
		cerr << "ERROR: Incompatible trajectory file size!" << endl;
		exit(EXIT_FAILURE);
	}
	for (size_t i=0; i<firstTrj.size(); i++)
	{
		trjTables.push_back(firstTrj[i]);
	}
	for (size_t i=0; i<secondTrj.size(); i++)
	{
		trjTables.push_back(secondTrj[i]);
	}
	return trjTables;
}

//this function will append second trj to the first one by column
vector<vector<vector<int> > > mergeTrjs(vector<vector<vector<int> > > & firstTrj, vector<vector<vector<int> > > & secondTrj)
{
	vector<vector<vector<int> > > trjTables;

	if(firstTrj.size() != secondTrj.size() 
		|| firstTrj[0].size() != secondTrj[0].size()
		|| firstTrj[0][0].size() != secondTrj[0][0].size())
	{
		cerr << "ERROR: Incompatible trajectory file size!" << endl;
		exit(EXIT_FAILURE);
	}
	trjTables.resize(firstTrj.size());
	for(size_t trjIndex = 0; trjIndex<firstTrj.size(); trjIndex++)
	{
		trjTables[trjIndex].resize(firstTrj[trjIndex].size());
	}
	for(size_t trjIndex = 0; trjIndex<firstTrj.size(); trjIndex++)
	{
		for(size_t rowIndex = 0; rowIndex < firstTrj[trjIndex].size(); rowIndex++)
		{
			trjTables[trjIndex][rowIndex].resize(2*firstTrj[trjIndex][rowIndex].size());
		}
	}
	for(size_t trjIndex = 0; trjIndex<trjTables.size(); trjIndex++)
	{
		for(size_t rowIndex = 0; rowIndex < trjTables[trjIndex].size(); rowIndex++)
		{
			for(size_t colIndex = 0; colIndex < trjTables[trjIndex][rowIndex].size(); colIndex++)
			{
				if(colIndex<firstTrj[trjIndex][rowIndex].size())
				{
					trjTables[trjIndex][rowIndex][colIndex] = firstTrj[trjIndex][rowIndex][colIndex];
				}else{
					trjTables[trjIndex][rowIndex][colIndex] = secondTrj[trjIndex][rowIndex][colIndex-firstTrj[trjIndex][rowIndex].size()];
				}
			}
		}
	}
	return trjTables;
}


//this function will splitTrjs by rows
vector<vector<vector<int> > > splitTrjsByRow(vector<vector<vector<int> > > & trjTable, int begin, int end)
{
	vector<vector<vector<int> > > splittrjTables;
	if(begin>end)
	{
		cerr << "ERROR: Begin index must be smaller than end index!" << endl;
		exit(EXIT_FAILURE);
	}
	for(int i=begin; i<end; i++)
	{
		splittrjTables.push_back(trjTable[i]);
	}
	return splittrjTables;
}


//this function will splitTrjs by column(gene name)
vector<vector<vector<int> > > splitTrjsByCol(vector<vector<vector<int> > > & trjTable, int begin, int end)
{
	vector<vector<vector<int> > > splittrjTables;
	splittrjTables.resize(trjTable.size());

	if(begin>end)
	{
		cerr << "ERROR: Begin index must be smaller than end index!" << endl;
		exit(EXIT_FAILURE);
	}
	for(size_t trjIndex = 0; trjIndex<splittrjTables.size(); trjIndex++)
	{
		splittrjTables[trjIndex].resize(trjTable[trjIndex].size());
	}
	for(size_t trjIndex = 0; trjIndex<splittrjTables.size(); trjIndex++)
	{
		for(size_t rowIndex = 0; rowIndex < splittrjTables[trjIndex].size(); rowIndex++)
		{
			splittrjTables[trjIndex][rowIndex].resize(end - begin);
		}
	}
	for(size_t trjIndex = 0; trjIndex<splittrjTables.size(); trjIndex++)
	{
		for(size_t rowIndex = 0; rowIndex < splittrjTables[trjIndex].size(); rowIndex++)
		{
			for(size_t colIndex = 0; colIndex < splittrjTables[trjIndex][rowIndex].size(); colIndex++)
			{
				splittrjTables[trjIndex][rowIndex][colIndex] = trjTable[trjIndex][rowIndex][colIndex + begin];
			}
		}
	}

	return splittrjTables;
}



bool TrajectoryCollection::CheckTraDeterministic() const
{
    //First, make histogram of each row in the trajectory collection
    //Then, check to see is rows with same values have same next value

    //multimap< state, pair< m_trajTables number, row number> >
    multimap< vector<int>, pair<unsigned,unsigned> > indexMap;

    //map< state, count >
    map< vector<int>, int> countMap;

    for( unsigned it = 0; it<m_trajTables.size(); ++it) //Traverse each trajectory
    {
        for( unsigned row = 0; row < m_trajTables[it].size(); ++row) //Traverse each row in the trajectory
        {
            indexMap.insert( pair< vector<int>,pair<unsigned,unsigned> > ( m_trajTables[it][row], pair<unsigned,unsigned>(it, row) ) );

            if( countMap.find( m_trajTables[it][row] ) == countMap.end() ) //row m_trajTables[it][row] hasn't shown up before, initialize its count to 1
                countMap[ m_trajTables[it][row] ] = 1;
            else //row m_trajTables[it][row] showed up before, increase its count by 1
                countMap[ m_trajTables[it][row] ]++;
        }
    }



    for(map< vector<int>, int>::iterator it = countMap.begin(); it!=countMap.end(); ++it)
    {
        if( (*it).second > 1)//Check to see is a row appears more than one time
        {
            vector<int> nextRowValue;
            multimap< vector<int>, pair<unsigned,unsigned> >::iterator i;
            for( i = indexMap.equal_range( (*it).first ).first; i!=indexMap.equal_range( (*it).first ).second; ++i)
            {
                if( (*i).second.second < m_trajTables[(*i).second.first ].size()-1) //Make sure it's not the last row of the trajectory
                {
                    nextRowValue =  m_trajTables[ (*i).second.first][ (*i).second.second + 1];
                    break;
                }

            }

            for( i = indexMap.equal_range( (*it).first ).first; i!=indexMap.equal_range( (*it).first ).second; ++i)
            {
                if( (*i).second.second < m_trajTables[(*i).second.first ].size()-1)
                {
                    if ( m_trajTables[ (*i).second.first][ (*i).second.second + 1] != nextRowValue )
                        return false;
                }
            }
        }
    }

    return true;  //Trajectory collection is pure
}