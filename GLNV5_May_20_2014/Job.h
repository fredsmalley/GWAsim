// Job.h
//
// Curtis Luce, Joe Song
// Created: 2007
// Modified: October 1, 2008
// Last modified: November 28, 2008. Joe Song. Removed dependence on JobResult by changing
//   the prototype of process().

#pragma once

#include "JobEnvironment.h"

class Job {

public:

	virtual void first(const JobEnvironment &env) = 0;

	virtual void next(const JobEnvironment &env) = 0;

	// virtual void transfer(JobEnvironment &env) = 0;

	virtual bool process(JobEnvironment &env) const = 0;

	virtual bool isTerminal() const = 0;

	virtual int pack(unsigned char* & jobBuffer) const = 0; 
	virtual void unpack(unsigned char * jobBuffer, int length)=0;

	virtual void print() const = 0;
	virtual void brief() const = 0;

};
