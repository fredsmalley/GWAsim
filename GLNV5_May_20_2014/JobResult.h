// JobResult.h
// 
// Curtis Luce, Joe Song
// Created: June 2007
// Last modified: November 28, 2008. Joe Song.
#pragma once

#include "JobEnvironment.h"
#include "Job.h"

 class JobResult{
	 
 public:

	virtual void merge(JobEnvironment &env) const = 0;

	virtual void extract(const Job & job, const JobEnvironment & env) = 0;

	virtual int pack(unsigned char* & jobResultBuffer) const = 0;

	virtual void unpack(unsigned char * jobResultBuffer, int length) = 0;
 
 };
