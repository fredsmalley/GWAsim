//
//  GLN_manpage.cpp
//  gln
//
//  Created by Joe Song on 10/1/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "GLN_copyright.h"

void manpage(char *command)
{
	cout << "NAME" << endl;
	cout << "\t" << command << " - Generalized logical nework reconstruction, comparison, simulation,\n\tand pathway adaptation analysis" << endl;
    
	cout << "SYNOPSIS" << endl;
    
	cout << "\t" << command << " -? \t (print help information)" << endl;
	cout << "\t" << command << " -h \t (print help information)" << endl;
	cout << "\t" << command << " -V \t (print copyright and version information)" << endl;
    
	cout << "\t" << command << endl
	<< "\t -M <mode> \t (estimation, comparison, simulation)" << endl
	<< "\t [-g <min_number_of_parents>] \t (must >=1, default = 1)" << endl
	<< "\t -p <max_number_of_parents>" << endl
	<< "\t [-J <min_Markovian_order>] \t (default = -1)" << endl
	<< "\t -K <max_Markovian_order>" << endl
	<< "\t [-P <p_value_mode>] \t (1,2,3(default),4,5)" << endl
	<< "\t -A <alpha_level> \t default=0.05" << endl
	<< "\t -R <trajectory_collection_file>" << endl
	<< "\t -O <output_GLN_file>" << endl
	<< "\t [-D <output_DOT_file>]" << endl
	<< "\t [-U <permutation_statistics_table>]" << endl
	<< "\t [-Z <0,1>] \t (allow self cycle, default: not allowed)" << endl
	<< "\t [-X <result_directory>]" << endl
	<< "\t [-1 <trajectory_collection_file>]" << endl
	<< "\t ..." << endl
	<< "\t [-9 <trajectory_collection_file>]" << endl
	<< "\t [-T <file_of_trajectory_collection_file_names>]" << endl
	<< "\t [-x <number_of_trajectories>]" << endl
	<< "\t [-S <number_of_time_points_on_each_trajectory>]" << endl
	<< "\t [-f <noise_level>] \t (default=0)" << endl
	<< "\t [-i <initial_trajectory_collection_file>]" << endl
    
	<< "\t [-Q DREAM3] \t (special option for DREAM3 Challenge 4)" << endl
	
	<< endl;
	
	cout << "OVERVIEW" << endl;
	cout 
	<< "\tGeneralized logical network (GLN) modeling software provides\n" 
	<< "\tdata-driven computational modeling for interactions in discrete\n" 
	<< "\ttime and discrete value dynamical systems." << endl;
	cout << endl;
	
	// cout << "DESCRIPTION" << endl;
	
	cout << "COPYRIGHT" << endl;
    
	cout << "\t" << GLNVersion << endl;
	cout << "\t" << Copyright << endl;
	cout << "\t" << Affiliation << endl;
    
}

void manpage_old(char *command) 
{
    cout << command << endl
    << "\t -M <mode>" << endl
    << "\t [-g <min# parents (>=1), default = 1>]" << endl
    << "\t -p <max# parents>" << endl
    << "\t [-J <min Markovian order, default = -1>]" << endl
    << "\t -K <max Markovian order>" << endl
    << "\t [-P <pval mode: 1,2,3(default),4,5>]" << endl
    << "\t -A <alpha level, default=0.05>" << endl
    << "\t -R <trajectory file>" << endl
    << "\t -O <output GLN file>" << endl
    << "\t [-D <output DOT file>]" << endl
    << "\t [-U <permutation statistics table>]" << endl
    << "\t [-Z <0,1>] (allow self cycle, default: not allowed)" << endl
    << "\t [-X <result directory>]" << endl
    << "\t [-Q DREAM3] (for DREAM3 Challenges 4)" << endl
    << "\t [-V] (print version information)" << endl
    << "\t [-1] ... [-9]: Trajectory collection files" << endl
    << "\t [-T] ... " << endl
    << "\t [-x <number of trajectories>]" << endl
    << "\t [-S <number of time points on each trajectory>]" << endl
    << "\t [-f <level of noise=0>]" << endl
    << "\t [-i <initial trajectory file>]" << endl
    << endl;
    
	/********* options I will need ***********
     *For estimation mode:
     *number of simulations to run
     *pvalue mode -- 1 or 2
     *pvalue's alpha (default at .05)
     *Restriction on the maximum parents to find?
     *max number of parents to try out (IF YES)
     *A file or list of trajectories it will need to read in
     *Must self calculate how many files were inputted
     *Output type -- print on screen(1) or to a file (2)
     *(IF output type =2) need the name of the outputfile
     *(ADD) -- print out the network? If yes, to what file?
     *-K need a max offset to estimate a trajectory file
     *-X will be the directory where Result files will be created
     
     *For simulation mode:
     *number of simulations to run
     *Which file that contains a network do you want to read in
     *Output type -- print on screen(1) or to a file (2)
     *(IF output type =2) need the name of the outputfile
     *-i will need a initial trajectory file for offsets greater than -1
     *also need a initial trajectory file if user wants to simulate
     *more then one trajectory
     
     *For generate(network) mode:
     *Size of the network
     *The minimum base for the nodes of the network
     *The maximum base for the nodes of the network
     *The maximum number of parents a node can have
     *-K max number of offsets to generate up to
     *-t the max number of trajectories to produce
     *-i a trajectory file to print initial trajectories to
     
     *For randomstates mode:
     *-R A logical network file name
     *-i A trajectory file name to print to
     *-K will need the max offset
     *-t number of trajectories to produce
     
     */
    
    
}
