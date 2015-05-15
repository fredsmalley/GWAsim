// GLN_evaluate.cpp
//
// Curtis Luce
// Created: 2007
// Modified: September 7, 2008.  Joe Song
// Last Modified: September 27, 2008.  Joe Song. 
//   File named changed from lognet_evaluat to GLN_evaluate.cpp.

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ofstream;

#include <vector>
using std::vector;

#include "GLN.h"

void summary(vector< vector<int> > &diff, int numNodes, string summaryfile);
vector< vector<int> > compareTrajectoriesValues(string traj1, string traj2, string difffile, int offset);

void readTrajForCompare(string filename, int& offset, vector< vector<int> >& traj_table, ofstream& out, int & num_traj,
						int &numb_of_init, vector<int> &numbTraj, vector<unsigned short> &maxsteps, vector< vector<unsigned short> > &externalNodes, int &intNodes, int &extNodes);

void GeneralizedLogicalNetwork::evaluate(int mode, string trajectoryfile, string outputfile, int offset, string differencefile, string summaryfile)
{
	/*******can simulate more then one trajectory and can do k-th order********/
	//Mode 1 is for standard output for the trajectory
	//Mode 2 is for file output for the trajectory
	ofstream out(outputfile.c_str()); 	//temp output
	vector<int> temp;	//this was moved out of the loop -- should be ok
	vector< vector<int> > traj_list;
	
	vector< vector<int> > trajectory;

	int numTraj = 0;
	// static int flag = 1;
	static int traj_used = 0;
	vector<int> num_of_init_traj;
	vector<unsigned short> maxsteps;
	vector< vector<unsigned short> > externalNodes;
	vector< vector<int> > trajdifferences;
	int intNodes = 0;
	int extNodes = 0;
	int extRow = 0;
	readTrajForCompare(trajectoryfile, offset, traj_list, out, numTraj, traj_used, num_of_init_traj, maxsteps, externalNodes,
		intNodes, extNodes);
	//cout << "after readTraj" << endl;
	// cout << "number of trajectories is:  " << numTraj << endl;
	traj_used = offset;
	// cout << "traj_used is:  " << traj_used << endl;
	vector<unsigned short> check;
	static int row = 0;  //rows of table to copy to traj table
	static int currow = 0;  //current row in traj table
	// cout << "starting the first for" << endl;

	for (int a = 0; a < numTraj; a++){
		// cout << "maxsteps at:  " << a << ";" << maxsteps[a] << endl;
		//adjust the row if offset is less than the number of initial trajectories given
		// cout << "initial row is: "  << row << endl;
		row = num_of_init_traj[a]+ row - offset;
		if(offset < traj_used){
			cout << "WARNING: For trajectory " << a+1 << " only " << offset << " initial states are being used out of " << traj_used << endl;
			//adjust the max number of steps, there will be less to simulate since the unused
			//trajectories still get printed out, case for offset < numb of init traj
			if(0 == a){
				maxsteps[a] -= (traj_used-offset);
			}//end if
		}//end if

		//case for initial trajectory file, file = 1 and no file = 0
		//if(flag == 1){
		if(mode==1){
			//will print out corrected number of trajectories
			if (offset < traj_used)
				cout << maxsteps[a] + (traj_used - offset) << endl;
			else
				cout << maxsteps[a] << endl;
		}//end if
		else{
			//will print out corrected number of trajectories
			if (offset < traj_used)
				out << maxsteps[a] + (traj_used - offset) << endl;
			else
				out << maxsteps[a]  << endl;
		}//end else

		//copy initial trajectories out of table that are needed for the current
		//trajectory simulation
		// cout << "traj_used is:  " << traj_used << endl;
		// for (int i = 0; i < traj_used; i++){
		for (int i = 0; i < traj_used && i < maxsteps[a]; i++){  // Joe Song 12/8/2009
			//currow = row;
			// cout << "current row is:  " << currow << endl;
			vector<int>::iterator iter = traj_list[currow].begin();
			//print out all initial trajectories to trajectory file
			for (size_t j = 0; j < size(); j++){
				// cout << "num of nodes:  " << numNodes << endl;
				if(mode == 1) {cout << *iter << "\t";}
				else {out << *iter << "\t";}
				iter++;
			}//end for
			if(mode==1){cout << endl;}
			else{out << endl;}
			//case for initial trajectory is needed for simulation, row is the
			//rows needed for simulation
			if(currow == row){
				trajectory.push_back(traj_list[row]);
				row++;
			}//end if
			currow++;
			//if (mode == 1){cout << "\n";}
		}//end for

		//}//end if

		// for(int t=offset-1;t < maxsteps[a]-offset; t++) {  Joe Song, 12/8/2007: changed to the following:
		for(int t=offset; t < maxsteps[a]; t++) {
			vector<int> temp;
			//for(int i = 0; i < numNodes; i++) {

			temp = simulate(trajectory, t);

			/*
			for(int i = 0; i < intNodes; i++){
				// cout << << "int nodes:  " << intNodes << "i:" << i << endl;
				temp.push_back(m_nodes[i].predict(trajectory));

				if(mode == 1){ cout << temp.back() << "\t"; }	//Standard output  (mode 1)
				else{ out << temp.back() << "\t"; }		//output to a file (mode 2)
			}*/

			for(int i = 0; i < extNodes; i++){
				temp.push_back(externalNodes[extRow][i]);
				if(mode == 1){cout << externalNodes[extRow][i] << "\t";}
				else{out << externalNodes[extRow][i] << "\t";}
			}//end for

			trajectory.push_back(temp);
			if(mode==1){ cout << "\n"; }
			else{ out << "\n"; }

			for (size_t i=0; i < size(); i++){
				m_nodes[i].setValue(temp[i]);
			}//end for
			extRow++;
			//******Doesn't temp need to be cleared??
		}//end outside for
		out << endl;
		m_current_time++;
		trajectory.clear();
		//case for offset < number of trajectories, increment row to next set of initial trajectories
		// cout << << "traj used is:  " << traj_used << endl;
		while (row < traj_used)
			row++;
		// cout << << "row is:  " << row << endl;
	}//end for

	out.close();

	trajdifferences = compareTrajectoriesValues(trajectoryfile, outputfile, differencefile, offset);
	// cout << << "original differences" << endl;
	for (size_t i = 0; i < trajdifferences.size(); i++){
		for (size_t j = 0; j < size(); j++){
			// cout << << trajdifferences[i][j] << " ";

		}
		// cout << << endl;
	}
	summary(trajdifferences, (int)size(), summaryfile);
}
