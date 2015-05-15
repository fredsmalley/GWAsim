// JobResultCmp.h -- The job result class for comparison based modeling
// 
// Joe Song
// Created: November 24, 2008
// Last modified: October 28, 2009.  Added working zone statistics.

#pragma once

#include <vector>
using std::vector;

#include "JobResult.h"
#include "JobEnvironment.h"
#include "TransitionTable.h"

class JobResultCmp:public JobResult {

private:

	TransitionTable m_diffTransTable;
	TransitionTable m_pooledTransTable;
	
	/*keep the best individual transition tables,since heterogeneity test can have multiple input trajectory files.*/
	vector<TransitionTable> m_TransTables;  

	// Total interaction strength statistics
	double m_chisq_t;
	int m_df_t;
	double m_pchisq_t;

	int m_type;

	// Working zone change statistics
	double m_chisq_z;
	int m_df_z;
	double m_pchisq_z;

public:

	virtual void merge(JobEnvironment &env) const;

	virtual void extract(const Job & job, const JobEnvironment & env);

	virtual int pack(unsigned char* & buffer) const;

	virtual void unpack(unsigned char * buffer, int length);

};

