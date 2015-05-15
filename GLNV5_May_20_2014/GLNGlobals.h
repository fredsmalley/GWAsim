//
//  GLNGlobals.h
//  gln
//
//  Created by Joe Song on 10/2/11.
//  Copyright 2011 New Mexico State University. All rights reserved.
//

#ifndef gln_GLNGlobals_h
#define gln_GLNGlobals_h

#include "GLN.h"
#include "GLNConfig.h"
#include "Topology.h"
#include "EnvGLNCmp.h"

GeneralizedLogicalNetwork ReconstructGLN(const GLNConfig & gc);

//Added by Haizhou Wang, June 18, 2012
/*
*ReconstructionPureGLN applies when the trajectory to be reconstructed is noise free.
In other word, the trajectory is deterministic.
Deterministism is defined as if given any row in the trajectory, the successor row of all its occurance will be the same.
*/
GeneralizedLogicalNetwork ReconstructPureGLN( const GLNConfig & gc );
GeneralizedLogicalNetwork ReconstructPureGLN( const TrajectoryCollection & traCol, const bool & selfControl, const Topology & mustIncludeParents);

//GeneralizedLogicalNetwork ReconstructPureGLN(const TrajectoryCollection & traCol, const vector<int> & radixs, const vector<string> & nodeNames, const vector<bool> & beParent, const vector<bool> & haveParent, const vector< vector<int> > & excludeLists, const bool & selfControl);


bool CheckTraDeterministic(const GLNConfig & gc);

bool CheckTraDeterministic(const vector< vector< vector<int> > > & tra);
//End of Haizhou's Code, June 18, 2012

//row corresponds to timesteps, col is genes.
enum ResampleType {ON_ROW=0, ON_COL=1, ON_COND=2, BOOTSTRAP=3, BOOTSTRAP_PLUS_ONE=4};

vector<TrajectoryCollection>
resample(const vector<TrajectoryCollection> & trajCols, ResampleType strategy);

EnvGLNCmp CompareGLNs(const GLNConfig & gc);
void CompareGLNs(EnvGLNCmp & envman);

EnvGLNCmp BuildTopologyThenCompareGLNs(const GLNConfig & gc);

EnvGLNCmp CompareGLNsResample(const GLNConfig & gc, const string & option="BY_MEM");

void ComparePathways(GLNConfig & gc, int permuteStrategy, const vector<string> & input_files);

void GLNTest();
void manpage(char *command);

void GLNExit(int status, const char * message = NULL);

// General-purpose functions

int randomNumberGenerator(int minNumber, int maxNumber);

void print_adjusted(const GeneralizedLogicalNetwork & net, const vector<TransitionTable> & tts, double alpha); 

#endif
