// JobResultRec.h -- The class for reconstruction job
// 
// Curtis Luce, Joe Song
// Created: June 2007
// Modified: September 28, 2008. Joe Song.  Removed dependency on JobEnvironmentGLN.h
// Last modified: Joe Song. November 30, 2008.  Class name changed from JobResultGLN to JobResultRec
#pragma once

#include <vector>
using std::vector;

#include "JobResult.h"
#include "JobEnvironment.h"
#include "TransitionTable.h"

class JobResultRec: public JobResult, public TransitionTable {

public:

	virtual void merge(JobEnvironment &env) const;

	virtual void extract(const Job & job, const JobEnvironment & env);

	virtual int pack(unsigned char* & buffer) const;

	virtual void unpack(unsigned char * buffer, int length);

};

